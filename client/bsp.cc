#include "bsp.hh"
#include "utils.hh"
#include "shaders.hh"
#include <fstream>
#include <glm/gtc/type_ptr.hpp>

const float feps = 1e-4f, surface_clip_eps = 0.125f;

bsp_plane::bsp_plane()
  : normal({ 0.f, 0.f, 0.f })
  , dist(0.0f) {
}

bsp_vertex bsp_vertex::operator+(const bsp_vertex &v) const {
  bsp_vertex res;
  res.position = position + v.position;
  res.decal = decal + v.decal;
  res.lightmap = lightmap + v.lightmap;
  res.normal = normal + v.normal;
  return res;
}

bsp_vertex bsp_vertex::operator*(float factor) const {
  bsp_vertex res;
  res.position = position * factor;
  res.decal = decal * factor;
  res.lightmap = lightmap * factor;
  res.normal = normal * factor;
  return res;
}

enum class lump {
  entities = 0,
  textures,
  planes,
  nodes,
  leaves,
  leaffaces,
  leafbrushes,
  models,
  brushes,
  brushsides,
  vertices,
  meshverts,
  effects,
  faces,
  lightmaps,
  lightvols,
  visdata
};

enum class face {
  polygon = 1,
  patch,
  mesh,
  billboard
};

struct bsp_direntry {
  int offset;
  int length;
};

struct bsp_header {
  char magic[4];
  int version;
  bsp_direntry direntries[17];
};

template <typename T>
void load_lump(std::ifstream &ifs, const bsp_header *header
    , lump lump_type, std::vector<T> &container) {
  int num_elements = header->direntries[(int)lump_type].length / sizeof(T);
  container.reserve(num_elements);
  ifs.seekg(header->direntries[(int)lump_type].offset);
  for (int i = 0; i < num_elements; ++i) {
    T element;
    ifs.read((char*)&element, sizeof(T));
    container.push_back(element);
  }
}

void bsp::_load_file(const char *filename, float world_scale
    , int tesselation_level) {
  std::ifstream ifs(filename, std::ios::binary);

  bsp_header header;
  ifs.read((char*)&header, sizeof(bsp_header));
  std::string magic = std::string(header.magic, 4);
  if (magic != "IBSP") {
    ifs.close();
    die("\"%s\": invalid magic \"%s\"", filename, magic.c_str());
  }
  if (header.version != 0x2e) {
    ifs.close();
    die("\"%s\": invalid version %d = 0x%x", filename, header.version
        , header.version);
  }

  load_lump(ifs, &header, lump::textures, _textures);
  texture_ids.resize(_textures.size());

  load_lump(ifs, &header, lump::planes, _planes);
  for (bsp_plane &p : _planes) {
    std::swap(p.normal.y, p.normal.z);
    p.normal.z = -p.normal.z;
    p.dist /= world_scale;
  }

  load_lump(ifs, &header, lump::nodes, _nodes);
  for (bsp_node &n : _nodes) {
    n.mins /= world_scale;
    std::swap(n.mins.y, n.mins.z);
    n.mins.z = -n.mins.z;
    n.maxs /= world_scale;
    std::swap(n.maxs.y, n.maxs.z);
    n.maxs.z = -n.maxs.z;
  }

  load_lump(ifs, &header, lump::leaves, _leaves);
  for (bsp_leaf &l : _leaves) {
    l.mins /= world_scale;
    std::swap(l.mins.y, l.mins.z);
    l.mins.z = -l.mins.z;
    l.maxs /= world_scale;
    std::swap(l.maxs.y, l.maxs.z);
    l.maxs.z = -l.maxs.z;
    if (l.mins.y > l.maxs.y)
      std::swap(l.mins.y, l.maxs.y);
    if (l.mins.z > l.maxs.z)
      std::swap(l.mins.z, l.maxs.z);
  }

  load_lump(ifs, &header, lump::leaffaces, _leaffaces);

  load_lump(ifs, &header, lump::leafbrushes, _leafbrushes);

  load_lump(ifs, &header, lump::brushes, _brushes);

  load_lump(ifs, &header, lump::brushsides, _brushsides);

  load_lump(ifs, &header, lump::vertices, _vertices);
  for (bsp_vertex &v : _vertices) {
    v.position /= world_scale;
    std::swap(v.position.y, v.position.z);
    v.position.z = -v.position.z;
  }

  load_lump(ifs, &header, lump::meshverts, _meshverts);

  load_lump(ifs, &header, lump::faces, _faces);
  _visible_faces.resize(_faces.size(), 0);

  int patch_count = 0
    , patch_size = (tesselation_level + 1) * (tesselation_level + 1)
    , patch_index_size = tesselation_level * tesselation_level * 6;
  for (const bsp_face &f : _faces)
    if (f.type == (int)face::patch)
      patch_count += ((f.size[0] - 1) / 2) * ((f.size[1] - 1) / 2);

  size_t vertex_count = _vertices.size(), meshverts_count = _meshverts.size();
  _vertices.resize(_vertices.size() + patch_count * patch_size);
  _meshverts.resize(_meshverts.size() + patch_count * patch_index_size);
  for (size_t i = 0, voff = vertex_count, eoff = meshverts_count; i < _faces.size(); ++i)
    if (_faces[i].type == (int)face::patch) {
      int dim_x = (_faces[i].size[0] - 1) / 2, dim_y = (_faces[i].size[1] - 1) / 2;
      _faces[i].meshvert = eoff;
      for (int x = 0, n = 0; n < dim_x; x = 2 * (++n))
        for (int y = 0, m = 0; m < dim_y; y = 2 * (++m)) {
          _tesselate(tesselation_level
              , _faces[i].vertex + x + _faces[i].size[0] * y, _faces[i].size[0]
              , voff, eoff);
          voff += patch_size;
          eoff += patch_index_size;
        }
      _faces[i].n_meshverts = eoff - _faces[i].meshvert;
    } else
      for (int j = 0; j < _faces[i].n_meshverts; ++j)
        _meshverts[_faces[i].meshvert + j].offset += _faces[i].vertex;

  load_lump(ifs, &header, lump::lightmaps, _lightmaps);
  _lightmap_texture_ids.resize(_lightmaps.size());

  ifs.seekg(header.direntries[(int)lump::visdata].offset);
  ifs.read((char*)&_visdata.n_vecs, sizeof(int));
  ifs.read((char*)&_visdata.sz_vecs, sizeof(int));
  int size = _visdata.n_vecs * _visdata.sz_vecs;
  _visdata.vecs.resize(size);
  ifs.read((char*)&_visdata.vecs[0], size * sizeof(unsigned char));

  ifs.close();
}

void bsp::_tesselate(int tesselation_level, int control_offset
    , int control_width, int voff, int eoff) {
  bsp_vertex controls[9];
  for (int c = 0, c_idx = 0; c < 3; c++) {
    int pos = c * control_width;
    controls[c_idx++] = _vertices[control_offset + pos];
    controls[c_idx++] = _vertices[control_offset + pos + 1];
    controls[c_idx++] = _vertices[control_offset + pos + 2];
  }

  int L1 = tesselation_level + 1;

  for (int j = 0; j <= tesselation_level; ++j) {
    float a = (float)j / tesselation_level, b = 1.f - a;
    _vertices[voff + j] = controls[0] * b * b + controls[3] * 2 * b * a
      + controls[6] * a * a;
  }

  for (int i = 1; i <= tesselation_level; ++i) {
    float a = (float)i / tesselation_level, b = 1.f - a;
    bsp_vertex temp[3];
    for (int j = 0; j < 3; ++j) {
      int k = 3 * j;
      temp[j] = controls[k] * b * b + controls[k + 1] * 2 * b * a + controls[k + 2] * a * a;
    }
    for (int j = 0; j <= tesselation_level; ++j) {
      float n_a = (float)j / tesselation_level, n_b = 1.f - n_a;
      _vertices[voff + i * L1 + j] = temp[0] * n_b * n_b
        + temp[1] * 2 * n_b * n_a + temp[2] * n_a * n_a;
    }
  }

  for (int i = 0; i <= tesselation_level; ++i)
    for (int j = 0; j <= tesselation_level; ++j) {
      int offset = eoff + (i * tesselation_level + j) * 6;
      _meshverts[offset + 0].offset = (i    ) * L1 + (j    ) + voff;
      _meshverts[offset + 1].offset = (i    ) * L1 + (j + 1) + voff;
      _meshverts[offset + 2].offset = (i + 1) * L1 + (j + 1) + voff;
      _meshverts[offset + 3].offset = (i + 1) * L1 + (j + 1) + voff;
      _meshverts[offset + 4].offset = (i + 1) * L1 + (j    ) + voff;
      _meshverts[offset + 5].offset = (i    ) * L1 + (j    ) + voff;
    }
}

int bsp::_find_leaf(glm::vec3 position) {
  int index = 0;
  while (index >= 0) {
    bsp_node *node = &_nodes[index];
    bsp_plane *plane = &_planes[node->plane];
    if (plane->normal.x * position.x + plane->normal.y * position.y
        + plane->normal.z * position.z > plane->dist)
      index = node->front;
    else
      index = node->back;
  }
  return -(index + 1); // leaf index
}

int bsp::_cluster_visible(int vis_cluster, int test_cluster) {
  if (vis_cluster < 0 || _visdata.vecs.size() == 0)
    return 1;
  int i = vis_cluster * _visdata.sz_vecs + (test_cluster >> 3);
  unsigned char vis_set = _visdata.vecs[i];
  return vis_set & (1 << (test_cluster & 7));
}

void bsp::_set_visible_faces(glm::vec3 camera_pos) {
  int leaf_index = _find_leaf(camera_pos);
  std::fill(_visible_faces.begin(), _visible_faces.end(), 0);
  for (bsp_leaf &l : _leaves)
    if (_cluster_visible(_leaves[leaf_index].cluster, l.cluster))
      for (int j = 0; j < l.n_leaffaces; j++)
        _visible_faces[_leaffaces[l.leafface + j].face] = 1;
}

bsp::bsp(const char *filename, float world_scale, int tesselation_level)
  : sp(shaders::map_vert, shaders::map_frag) {
  _load_file(filename, world_scale, tesselation_level);

  sp.use_this_prog();
  _vertex_pos_attr = sp.bind_attrib("vertex_pos");
  _texture_coord_attr = sp.bind_attrib("texture_coord");
  _lightmap_coord_attr = sp.bind_attrib("lightmap_coord");
  _mvp_mat_unif = sp.bind_uniform("mvp");
  glUniform1i(sp.bind_uniform("texture_sampler"), 0);
  glUniform1i(sp.bind_uniform("lightmap_sampler"), 1);

  vbo.bind();
  vbo.upload(sizeof(_vertices[0]) * _vertices.size(), &_vertices[0]);

  ebo.bind();
  ebo.upload(sizeof(_meshverts[0]) * _meshverts.size(), &_meshverts[0]);

  glEnableVertexAttribArray(_vertex_pos_attr);
  glEnableVertexAttribArray(_texture_coord_attr);
  glEnableVertexAttribArray(_lightmap_coord_attr);

  glActiveTexture(GL_TEXTURE0);
  for (size_t i = 0; i < _textures.size(); ++i)
    texture_ids[i] = 0;

  glActiveTexture(GL_TEXTURE1);
  glGenTextures(_lightmaps.size(), &_lightmap_texture_ids[0]);
  for (size_t i = 0; i < _lightmaps.size(); ++i) {
    glBindTexture(GL_TEXTURE_2D, _lightmap_texture_ids[i]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 128, 128, 0, GL_RGB
        , GL_UNSIGNED_BYTE, _lightmaps[i].map);
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }
  sp.dont_use_this_prog();
}

bsp::~bsp() {
  glDeleteTextures(_textures.size(), texture_ids.data());
  glDeleteTextures(_lightmaps.size(), _lightmap_texture_ids.data());
}

void bsp::render(glm::vec3 position, const glm::mat4 &mvp) {
  _set_visible_faces(position);

  sp.use_this_prog();

  glUniformMatrix4fv(_mvp_mat_unif, 1, GL_FALSE, glm::value_ptr(mvp));

  glVertexAttribPointer(_vertex_pos_attr, 3, GL_FLOAT, GL_FALSE
      , sizeof(bsp_vertex), (void*)(long)offsetof(bsp_vertex, position));
  glVertexAttribPointer(_texture_coord_attr, 2, GL_FLOAT, GL_FALSE
      , sizeof(bsp_vertex), (void*)(long)offsetof(bsp_vertex, decal));
  glVertexAttribPointer(_lightmap_coord_attr, 2, GL_FLOAT, GL_FALSE
      , sizeof(bsp_vertex), (void*)(long)offsetof(bsp_vertex, lightmap));

  for (size_t i = 0; i < _faces.size(); i++) {
    if (!_visible_faces[i])
      continue;
    if (_faces[i].type == (int)face::polygon
        || _faces[i].type == (int)face::mesh
        || _faces[i].type == (int)face::patch) {
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, texture_ids[_faces[i].texture]);
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, _lightmap_texture_ids[_faces[i].lm_index]);
      glDrawElements(GL_TRIANGLES, _faces[i].n_meshverts, GL_UNSIGNED_INT
          , (void*)(long)(_faces[i].meshvert * sizeof(GLuint)));
    }
  }

  sp.dont_use_this_prog();
}

void bsp::trace_sphere(trace_result *tr, const glm::vec3 &start
    , const glm::vec3 &end, float radius) {
  trace_description t { trace_type::sphere, radius };
  _trace(tr, t, start, end);
}

void bsp::_trace(trace_result *tr, const trace_description &td
    , const glm::vec3 &start, const glm::vec3 &end) {
  tr->fraction = 1.0f;
  tr->start_solid = tr->all_solid = false;

  _check_node(tr, td, 0, 0.0f, 1.0f, start, end);

  if (tr->fraction == 1.0f)
    tr->end = end;
  else
    tr->end = start + tr->fraction * (end - start);
}

// Q3: CM_TraceThroughTree
void bsp::_check_node(trace_result *tr, const trace_description &td
    , int node_index, float start_fraction, float end_fraction
    , const glm::vec3 &start, const glm::vec3 &end) {
  if (tr->fraction <= start_fraction)
    return;

  if (node_index < 0) {
    bsp_leaf *leaf = &_leaves[-(node_index + 1)];
    for (int i = 0; i < leaf->n_leafbrushes; ++i) {
      bsp_brush *b = &_brushes[_leafbrushes[leaf->leafbrush + i].brush];
      if (b->n_brushsides > 0)
        _check_brush(tr, td, b, start, end);
    }
    return;
  }

  bsp_node *node = &_nodes[node_index];
  bsp_plane *plane = &_planes[node->plane];
  float start_dist = glm::dot(start, plane->normal) - plane->dist
    , end_dist = glm::dot(end, plane->normal) - plane->dist, offset = 0;

  switch (td.type) {
    case trace_type::ray:
      offset = 0;
      break;
    case trace_type::sphere:
      offset = td.radius;
      break;
    default:
      offset = 0;
      break;
  }

  if (start_dist >= offset + 1 && end_dist >= offset + 1) {
    _check_node(tr, td, node->front, start_fraction, end_fraction, start, end);
    return;
  } else if (start_dist < -offset - 1 && end_dist < -offset - 1) {
    _check_node(tr, td, node->back, start_fraction, end_fraction, start, end);
    return;
  }

  int side;
  float fraction1, fraction2;
  if (start_dist < end_dist) {
    float inverse_dist = 1.0f / (start_dist - end_dist); // TODO
    side = 1; // back
    fraction1 = (start_dist - offset + surface_clip_eps) * inverse_dist;
    fraction2 = (start_dist + offset + surface_clip_eps) * inverse_dist;
  } else if (start_dist > end_dist) {
    float inverse_dist = 1.0f / (start_dist - end_dist);
    side = 0; // front
    fraction1 = (start_dist + offset + surface_clip_eps) * inverse_dist;
    fraction2 = (start_dist - offset - surface_clip_eps) * inverse_dist;
  } else {
    side = 0;
    fraction1 = 1.0f;
    fraction2 = 0.0f;
  }

  // TODO
  // fraction1 = clamp(fraction1, 0.f, 1.f);
  // fraction2 = clamp(fraction2, 0.f, 1.f);
  if (fraction1 < 0.0f)
    fraction1 = 0.0f;
  else if (fraction1 > 1.0f)
    fraction1 = 1.0f;
  if (fraction2 < 0.0f)
    fraction2 = 0.0f;
  else if (fraction2 > 1.0f)
    fraction2 = 1.0f;

  float middle_fraction = start_fraction + (end_fraction - start_fraction)
    * fraction1;

  glm::vec3 middle = start + fraction1 * (end - start);

  if (side == 0)
    _check_node(tr, td, node->front, start_fraction, middle_fraction, start
        , middle);
  else
    _check_node(tr, td, node->back, start_fraction, middle_fraction, start
        , middle);

  middle_fraction = start_fraction + (end_fraction - start_fraction) * fraction2;
  middle = start + fraction2 * (end - start);

  if (side == 0)
    _check_node(tr, td, node->back, middle_fraction, end_fraction, middle, end);
  else
    _check_node(tr, td, node->front, middle_fraction, end_fraction, middle, end);
}

// Q3: CM_TraceThroughBrush
void bsp::_check_brush(trace_result *tr, const trace_description &td
    , bsp_brush *b, const glm::vec3 &input_start, const glm::vec3 &input_end) {
  float start_fraction = -1.0f, end_fraction = 1.0f;
  bool get_out = false, starts_out = 0;
  for (int i = 0; i < b->n_brushsides; ++i) {
    bsp_brushside *brushside = &_brushsides[b->brushside + i];
    bsp_plane *plane = &_planes[brushside->plane];
    glm::vec3 offset;
    float start_dist, end_dist;

    switch (td.type) {
      case trace_type::ray:
        start_dist = glm::dot(input_start, plane->normal) - plane->dist;
        end_dist = glm::dot(input_end, plane->normal) - plane->dist;
        break;
      case trace_type::sphere:
        start_dist = glm::dot(input_start, plane->normal) - plane->dist - td.radius;
        end_dist = glm::dot(input_end, plane->normal) - plane->dist - td.radius;
        break;
      default:
        start_dist = glm::dot(input_start, plane->normal) - plane->dist;
        end_dist = glm::dot(input_end, plane->normal) - plane->dist;
        break;
    }

    if (end_dist > 0)
      get_out = 1;

    if (start_dist > 0)
      starts_out = 1;

    if (start_dist > 0 && (end_dist >= surface_clip_eps || end_dist >= start_dist))
      return;

    if (start_dist <= 0 && end_dist <= 0)
      continue;

    if (start_dist > end_dist) {
      float fraction = (start_dist - surface_clip_eps) / (start_dist - end_dist);
      if (fraction < 0)
        fraction = 0;
      if (fraction > start_fraction) {
        start_fraction = fraction;
        tr->clip_plane_normal = plane->normal;
      }
    } else {
      float fraction = (start_dist + surface_clip_eps) / (start_dist - end_dist);
      if (fraction > 1)
        fraction = 1;
      if (fraction < end_fraction)
        end_fraction = fraction;
    }
  }

  if (!starts_out) {
    tr->start_solid = true;
    if (!get_out) {
      tr->all_solid = true;
      tr->fraction = 0;
    }
    return;
  }

  if (start_fraction < end_fraction)
    if (start_fraction > -1 && start_fraction < tr->fraction) {
      if (start_fraction < 0)
        start_fraction = 0;
      tr->fraction = start_fraction;
    }
}

