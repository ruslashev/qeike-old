#include "bsp.hh"
#include "utils.hh"
#include "shaders.hh"
#include <fstream>
#include <glm/gtc/type_ptr.hpp>

namespace qke {

const float feps = 1e-4f, surface_clip_eps = 0.125f
    , grid_wrap_point_epsilon = 0.1f, subdivide_distance = 4
    , point_epsilon = 0.1f, plane_eps = 0.000001f, plane_tri_eps = 0.1f
    , normal_eps = 0.0001f, dist_eps = 0.02f, chop_eps = 0.1f;

const int max_patch_verts = 2048, bogus_range = 65535
    , max_points_on_winding = 96, max_map_bounds = 65535;

bsp::bsp_plane::bsp_plane()
  : normal({ 0.f, 0.f, 0.f })
  , dist(0.0f) {
}

bsp::bsp_vertex bsp::bsp_vertex::operator+(const bsp_vertex &v) const {
  bsp_vertex res;
  res.position = position + v.position;
  res.lightmap = lightmap + v.lightmap;
  return res;
}

bsp::bsp_vertex bsp::bsp_vertex::operator*(float factor) const {
  bsp_vertex res;
  res.position = position * factor;
  res.lightmap = lightmap * factor;
  return res;
}

enum class lump {
  entities = 0,
  shaders,
  planes,
  nodes,
  leaves,
  leaffaces, // Q3: leafsurfaces
  leafbrushes,
  models,
  brushes,
  brushsides,
  vertices,
  meshverts,
  fogs,
  faces,
  lightmaps,
  lightgrid,
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

struct bsp_raw_node {
  int plane;
  int front;
  int back;
  glm::ivec3 mins;
  glm::ivec3 maxs;
};

struct bsp_raw_leaf {
  int cluster;
  int area;
  glm::ivec3 mins;
  glm::ivec3 maxs;
  int leafface;
  int n_leaffaces;
  int leafbrush;
  int n_leafbrushes;
};

struct bsp_raw_vertex {
  glm::vec3 position;
  glm::vec2 decal;
  glm::vec2 lightmap;
  glm::vec3 normal;
  unsigned char color[4];
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

void bsp::_load_file(const char *filename, int tesselation_level) {
  std::ifstream ifs(filename, std::ios::binary);
  if (!ifs)
    die("failed to open map \"%s\"", filename);

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

  load_lump(ifs, &header, lump::shaders, _shaders);

  load_lump(ifs, &header, lump::planes, _planes);
  for (bsp_plane &p : _planes) {
    std::swap(p.normal.y, p.normal.z);
    p.normal.z = -p.normal.z;
  }

  std::vector<bsp_raw_node> raw_nodes;
  load_lump(ifs, &header, lump::nodes, raw_nodes);
  for (const bsp_raw_node &n : raw_nodes) {
    bsp_node conv_node;
    conv_node.plane = n.plane;
    conv_node.front = n.front;
    conv_node.back = n.back;
    conv_node.mins = glm::vec3(n.mins.x, n.mins.z, -n.mins.y);
    conv_node.maxs = glm::vec3(n.maxs.x, n.maxs.z, -n.maxs.y);
    _nodes.push_back(conv_node);
  }

  std::vector<bsp_raw_leaf> raw_leaves;
  load_lump(ifs, &header, lump::leaves, raw_leaves);
  for (const bsp_raw_leaf &l : raw_leaves) {
    bsp_leaf conv_leaf;
    conv_leaf.cluster = l.cluster;
    conv_leaf.area = l.area;
    conv_leaf.mins = glm::vec3(l.mins.x, l.mins.z, -l.mins.y);
    conv_leaf.maxs = glm::vec3(l.maxs.x, l.maxs.z, -l.maxs.y);
    conv_leaf.leafface = l.leafface;
    conv_leaf.n_leaffaces = l.n_leaffaces;
    conv_leaf.leafbrush = l.leafbrush;
    conv_leaf.n_leafbrushes = l.n_leafbrushes;
    _leaves.push_back(conv_leaf);
  }

  load_lump(ifs, &header, lump::leaffaces, _leaffaces);

  load_lump(ifs, &header, lump::leafbrushes, _leafbrushes);

  load_lump(ifs, &header, lump::brushes, _brushes);

  load_lump(ifs, &header, lump::brushsides, _brushsides);

  std::vector<bsp_raw_vertex> raw_vertices;
  load_lump(ifs, &header, lump::vertices, raw_vertices);
  for (const bsp_raw_vertex &v : raw_vertices) {
    bsp_vertex conv_vertex;
    conv_vertex.position = v.position;
    std::swap(conv_vertex.position.y, conv_vertex.position.z);
    conv_vertex.position.z = -conv_vertex.position.z;
    conv_vertex.lightmap = v.lightmap;
    _vertices.push_back(conv_vertex);
  }

  load_lump(ifs, &header, lump::meshverts, _meshverts);

  load_lump(ifs, &header, lump::faces, _faces);
  _visible_faces.resize(_faces.size(), 0);

  _create_patches(tesselation_level);

  for (size_t i = 0; i < _meshverts.size(); i += 3) {
    if (i + 2 >= _meshverts.size()) {
      warning("meshvert winding swap routine ended too soon");
      break;
    }
    std::swap(_meshverts[i + 1], _meshverts[i + 2]);
  }

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

void bsp::_create_patches(int tesselation_level) {
  int patch_count = 0
    , patch_size = (tesselation_level + 1) * (tesselation_level + 1)
    , patch_index_size = tesselation_level * tesselation_level * 6;
  for (const bsp_face &f : _faces)
    if (f.type == (int)face::patch)
      patch_count += ((f.size[0] - 1) / 2) * ((f.size[0] - 1) / 2);

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

void bsp::_set_visible_faces(const glm::vec3 &camera_pos) {
  int leaf_index = _find_leaf(camera_pos);
  std::fill(_visible_faces.begin(), _visible_faces.end(), 0);
  for (bsp_leaf &l : _leaves)
    if (_cluster_visible(_leaves[leaf_index].cluster, l.cluster))
      for (int j = 0; j < l.n_leaffaces; j++)
        _visible_faces[_leaffaces[l.leafface + j].face] = 1;
}

bsp::bsp(const char *filename, int tesselation_level)
  : _sp(shaders::map_vert, shaders::map_frag) {
  _load_file(filename, tesselation_level);

  _vao.bind();

  _sp.use_this_prog();
  _vertex_pos_attr = _sp.bind_attrib("vertex_pos");
  _lightmap_coord_attr = _sp.bind_attrib("lightmap_coord");
  _mvp_mat_unif = _sp.bind_uniform("mvp");
  glUniform1i(_sp.bind_uniform("lightmap_sampler"), 1);

  _vbo.bind();
  _vbo.upload(sizeof(_vertices[0]) * _vertices.size(), &_vertices[0]);

  _ebo.bind();
  _ebo.upload(sizeof(_meshverts[0]) * _meshverts.size(), &_meshverts[0]);

  glEnableVertexAttribArray(_vertex_pos_attr);
  glEnableVertexAttribArray(_lightmap_coord_attr);

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
  glVertexAttribPointer(_vertex_pos_attr, 3, GL_FLOAT, GL_FALSE
      , sizeof(bsp_vertex), (void*)offsetof(bsp_vertex, position));
  glVertexAttribPointer(_lightmap_coord_attr, 2, GL_FLOAT, GL_FALSE
      , sizeof(bsp_vertex), (void*)offsetof(bsp_vertex, lightmap));

  _sp.dont_use_this_prog();
  _vao.unbind();
  _vbo.unbind();
  _ebo.unbind();
}

bsp::~bsp() {
  glDeleteTextures(_lightmaps.size(), _lightmap_texture_ids.data());
}

void bsp::draw(const glm::vec3 &position, const glm::mat4 &mvp) {
  _set_visible_faces(position);

  _vao.bind();
  _sp.use_this_prog();

  glUniformMatrix4fv(_mvp_mat_unif, 1, GL_FALSE, glm::value_ptr(mvp));

  glActiveTexture(GL_TEXTURE1);
  for (size_t i = 0; i < _faces.size(); i++) {
    if (!_visible_faces[i])
      continue;
    if (_faces[i].type == (int)face::polygon
        || _faces[i].type == (int)face::mesh
        || _faces[i].type == (int)face::patch) {
      glBindTexture(GL_TEXTURE_2D, _lightmap_texture_ids[_faces[i].lm_index]);
      glDrawElements(GL_TRIANGLES, _faces[i].n_meshverts, GL_UNSIGNED_INT
          , (void*)(_faces[i].meshvert * sizeof(GLuint)));
    }
  }

  _vao.unbind();
  _sp.dont_use_this_prog();
}

void bsp::trace_sphere(trace_result *tr, const glm::vec3 &start
      , const glm::vec3 &end, float radius) {
  trace_description t { trace_type::sphere, radius };
  _trace(tr, t, start, end);
}

void bsp::_trace(trace_result *tr, const trace_description &td
    , const glm::vec3 &start, const glm::vec3 &end) {
  tr->clip_plane_normal = glm::vec3(0, 0, 0);
  tr->fraction = 1.0f;
  tr->end = end;
  tr->start_solid = tr->all_solid = false;

  _trace_node(tr, td, 0, 0.0f, 1.0f, start, end);

  if (tr->fraction == 1.0f)
    tr->end = end;
  else
    tr->end = start + tr->fraction * (end - start);
}

// Q3: CM_TraceThroughTree
void bsp::_trace_node(trace_result *tr, const trace_description &td
    , int node_index, float start_fraction, float end_fraction
    , const glm::vec3 &start, const glm::vec3 &end) {
  if (tr->fraction <= start_fraction)
    return;

  if (node_index < 0) {
    _trace_leaf(tr, td, &_leaves[-node_index - 1], start, end);
    return;
  }

  bsp_node *node = &_nodes[node_index];
  bsp_plane *plane = &_planes[node->plane];
  // TODO: plane->type not respected
  float start_dist = glm::dot(start, plane->normal) - plane->dist
    , end_dist = glm::dot(end, plane->normal) - plane->dist, offset = 0;

  offset = td.radius;

  if (start_dist >= offset + 1 && end_dist >= offset + 1) {
    _trace_node(tr, td, node->front, start_fraction, end_fraction, start, end);
    return;
  } else if (start_dist < -offset - 1 && end_dist < -offset - 1) {
    _trace_node(tr, td, node->back, start_fraction, end_fraction, start, end);
    return;
  }

  int side;
  float fraction1, fraction2;
  if (start_dist < end_dist) {
    const float inverse_dist = 1.0f / (start_dist - end_dist);
    side = 1; // back
    fraction1 = (start_dist - offset + surface_clip_eps) * inverse_dist;
    fraction2 = (start_dist + offset + surface_clip_eps) * inverse_dist;
  } else if (start_dist > end_dist) {
    const float inverse_dist = 1.0f / (start_dist - end_dist);
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
  if (fraction1 < 0.f)
    fraction1 = 0.f;
  else if (fraction1 > 1.f)
    fraction1 = 1.f;
  if (fraction2 < 0.f)
    fraction2 = 0.f;
  else if (fraction2 > 1.f)
    fraction2 = 1.f;

  float middle_fraction = start_fraction + (end_fraction - start_fraction)
    * fraction1;

  glm::vec3 middle = start + fraction1 * (end - start);

  if (side == 0)
    _trace_node(tr, td, node->front, start_fraction, middle_fraction, start
        , middle);
  else
    _trace_node(tr, td, node->back, start_fraction, middle_fraction, start
        , middle);

  middle_fraction = start_fraction + (end_fraction - start_fraction) * fraction2;
  middle = start + fraction2 * (end - start);

  if (side == 0)
    _trace_node(tr, td, node->back, middle_fraction, end_fraction, middle, end);
  else
    _trace_node(tr, td, node->front, middle_fraction, end_fraction, middle, end);
}

// Q3: CM_TraceThroughBrush
// TODO: put input_start and end to trace_description
void bsp::_trace_brush(trace_result *tr, const trace_description &td
    , const bsp_brush *b, const glm::vec3 &input_start
    , const glm::vec3 &input_end) {
  bsp_plane *clip_plane = nullptr;
  float start_fraction = -1.0f, end_fraction = 1.0f;
  if (!b->n_brushsides)
    return;
  bool get_out = false, starts_out = false;
  for (int i = 0; i < b->n_brushsides; ++i) {
    bsp_brushside *brushside = &_brushsides[b->brushside + i];
    bsp_plane *plane = &_planes[brushside->plane];
    float start_dist, end_dist;

    start_dist = glm::dot(input_start, plane->normal) - plane->dist - td.radius;
    end_dist = glm::dot(input_end, plane->normal) - plane->dist - td.radius;

    if (end_dist > 0)
      get_out = true;

    if (start_dist > 0)
      starts_out = true;

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
        clip_plane = plane;
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
      tr->clip_plane_normal = clip_plane->normal;
    }
}

void bsp::_trace_leaf(trace_result *tr, const trace_description &td
    , const bsp_leaf *l, const glm::vec3 &start
    , const glm::vec3 &end) {
  // trace line against all brushes in the leaf
  for (int k = 0; k < l->n_leafbrushes; ++k) {
    bsp_brush *b = &_brushes[_leafbrushes[l->leafbrush + k]];
    // if (b->checkcount == cm.checkcount)
    //   continue; // already checked this brush in another leaf
    // b->checkcount = cm.checkcount;
    if (!(_shaders[b->shader].contents & 1))
      continue;
    _trace_brush(tr, td, b, start, end);
    if (!tr->fraction)
      return;
  }
}

};

