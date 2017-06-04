#include "bsp.hh"
#include "utils.hh"
#include "shaders.hh"
#include <fstream>
#include <glm/gtc/type_ptr.hpp>

const float feps = 1e-4f;

bsp_plane::bsp_plane()
  : normal({ 0.f, 0.f, 0.f })
  , dist(0.0f) {
}

bsp_biquadratic_patch::bsp_biquadratic_patch()
  : _tesselation_level(0)
  , _triangles_per_row(nullptr)
  , _row_index_pointers(nullptr) {
}

bsp_biquadratic_patch::~bsp_biquadratic_patch() {
  delete [] _triangles_per_row;
  delete [] _row_index_pointers;
}

void bsp_biquadratic_patch::tesselate(int level) {
  _tesselation_level = level;
  vertices.resize((_tesselation_level + 1) * (_tesselation_level + 1));

  for (int i = 0; i <= _tesselation_level; ++i) {
    float a = (float)i / (float)_tesselation_level, b = 1.f - a;
    vertices[i] = control_points[0] * (b * b)
      + control_points[3] * (2 * b * a)
      + control_points[6] * (a * a);
  }

  for (int i = 1; i <= _tesselation_level; ++i) {
    float a = (float)i / (float)_tesselation_level, b = 1.f - a;
    bsp_vertex temp[3];
    for (int j = 0, k = 0; j < 3; ++j, k = 3 * j)
      temp[j] = control_points[k + 0] * b * b
        + control_points[k + 1] * 2 * b * a
        + control_points[k + 2] * a * a;
    for (int j = 0; j <= _tesselation_level; ++j) {
      float ta = (float)j / (float)_tesselation_level, tb = 1.f - ta;
      vertices[i * (_tesselation_level + 1) + j] = temp[0] * tb * tb
        + temp[1] * 2 * tb * ta
        + temp[2] * ta * ta;
    }
  }

  _indices.resize(_tesselation_level * (_tesselation_level + 1) * 2);
  for (int row = 0; row < _tesselation_level; ++row)
    for (int col = 0; col <= _tesselation_level; ++col) {
      _indices[(row * (_tesselation_level + 1) + col) * 2 + 1]
        =  row      * (_tesselation_level + 1) + col;
      _indices[(row * (_tesselation_level + 1) + col) * 2]
        = (row + 1) * (_tesselation_level + 1) + col;
    }

  _triangles_per_row = new int [_tesselation_level];
  _row_index_pointers = new unsigned int* [_tesselation_level];
  for (int row = 0; row < _tesselation_level; ++row) {
    _triangles_per_row[row] = 2 * (_tesselation_level + 1);
    _row_index_pointers[row] = &_indices[row * 2 * (_tesselation_level + 1)];
  }
}

void bsp_biquadratic_patch::render() const {
  for (int row = 0; row < _tesselation_level; ++row)
    glDrawElements(GL_TRIANGLE_STRIP, 2 * (_tesselation_level + 1), GL_UNSIGNED_INT,
        &_indices[row * 2 * (_tesselation_level + 1)]);
}

#if 0
bsp_plane::bsp_plane(glm::vec3 n_normal, float n_intercept)
  : normal(n_normal)
  , intercept(n_intercept) {
}

bsp_plane::bsp_plane(const bsp_plane &rhs) {
  normal = rhs.normal;
  intercept = rhs.intercept;
}

void bsp_plane::set_normal(const glm::vec3 &rhs) {
  normal = rhs;
}

void bsp_plane::set_intercept(float n_intercept) {
  intercept = n_intercept;
}

void bsp_plane::calculate_intercept(const glm::vec3 &point_on_plane) {
  intercept = -glm::dot(normal, point_on_plane);
}

void bsp_plane::set_from_points(const glm::vec3 &p0, const glm::vec3 &p1
    , const glm::vec3 &p2) {
  normal = glm::normalize(glm::cross(p1 - p0, p2 - p0));
  calculate_intercept(p0);
}

void bsp_plane::normalize() {
  float normal_len = normal.length();
  normal /= normal_len;
  intercept /= normal_len;
}

bool bsp_plane::intersect3(const bsp_plane &p2, const bsp_plane &p3
    , glm::vec3 &result) {
  float denominator = glm::dot(normal, glm::cross(p2.normal, p3.normal));
  if (abs(denominator) <= feps)
    return false;
  glm::vec3 temp1 = glm::cross(p2.normal, p3.normal) * intercept
    , temp2 = glm::cross(p3.normal, normal) * p2.intercept
    , temp3 = glm::cross(normal, p2.normal) * p3.intercept;
  result = (temp1 + temp2 + temp3) / (-denominator);
  return true;
}

float bsp_plane::get_distance(const glm::vec3 &point) const {
  return point.x * normal.x + point.y * normal.y + point.z * normal.z
    + intercept;
}

int bsp_plane::classify_point(const glm::vec3 &point) const {
  float distance = get_distance(point);
  if (distance > feps)
    return 1;
  if (distance < -feps)
    return -1;
  return 0;
}

bsp_plane bsp_plane::lerp(const bsp_plane &p2, float factor) {
  bsp_plane result;
  result.normal = glm::normalize(normal * (1.0f - factor) + p2.normal * factor);
  result.intercept = intercept * (1.0f - factor) + p2.intercept * factor;
  return result;
}

bool bsp_plane::operator==(const bsp_plane &rhs) const {
  return (glm::all(glm::lessThan(glm::abs(normal - rhs.normal), glm::vec3(feps)))
      && abs(intercept - rhs.intercept) < feps);
}

bool bsp_plane::operator!=(const bsp_plane &rhs) const {
  return !((*this) == rhs);
}

bsp_plane bsp_plane::operator-() const {
  return bsp_plane(-normal, -intercept);
}

bsp_plane bsp_plane::operator+() const {
  return (*this);
}
#endif

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

  load_lump(ifs, &header, lump::vertices, _vertices);
  for (bsp_vertex &v : _vertices) {
    v.position /= world_scale;
    std::swap(v.position.y, v.position.z);
    v.position.z = -v.position.z;
  }

  load_lump(ifs, &header, lump::meshverts, _meshverts);

  load_lump(ifs, &header, lump::faces, _faces);
  _visible_faces.resize(_faces.size(), 0);
  for (const bsp_face &f : _faces)
    if (f.type == (int)face::patch)
      _create_patch(f, tesselation_level);

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

void bsp::_create_patch(const bsp_face &f, int tesselation_level) {
  bsp_patch *new_patch = new bsp_patch;
  new_patch->texture_idx = f.texture;
  new_patch->lightmap_idx = f.lm_index;
  new_patch->width = f.size[0];
  new_patch->height = f.size[1];
  int num_patches_width  = (new_patch->width - 1) >> 1
    , num_patches_height = (new_patch->height - 1) >> 1;

  new_patch->quadratic_patches.resize(num_patches_width * num_patches_height);

  for (int y = 0; y < num_patches_height; ++y)
    for (int x = 0; x < num_patches_width; ++x) {
      for (int row = 0; row < 3; ++row)
        for (int col = 0; col < 3; ++col) {
          int patch_idx = y * num_patches_width + x, cp_idx = row * 3 + col
            , vert_idx = f.vertex + (y * 2 * new_patch->width + x * 2)
            + row * new_patch->width + col;
          new_patch->quadratic_patches[patch_idx].control_points[cp_idx]
            = _vertices[vert_idx];
        }
      new_patch->quadratic_patches[y * num_patches_width + x].tesselate(tesselation_level);
    }

  _patches.push_back(new_patch);
}

bsp::bsp(const char *filename, float world_scale, int tesselation_level)
  : sp(shaders::map_vert, shaders::map_frag) {
  _load_file(filename, world_scale, tesselation_level);

  sp.use_this_prog();
  _vertex_pos_attr = sp.bind_attrib("vertex_pos");
  _texture_coord_attr = sp.bind_attrib("texture_coord");
  _lightmap_coord_attr = sp.bind_attrib("lightmap_coord");
  _mvp_mat_unif = sp.bind_uniform("mvp");
  glUniform1i(glGetUniformLocation(sp.id, "texture_sampler"), 0);
  glUniform1i(glGetUniformLocation(sp.id, "lightmap_sampler"), 1);

  vbo.bind();
  vbo.upload(sizeof(_vertices[0]) * _vertices.size(), &_vertices[0]);

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
  /*
     glEnable(GL_CULL_FACE);
     glCullFace(GL_FRONT);
  */

  _set_visible_faces(position);

  sp.use_this_prog();

  glUniformMatrix4fv(_mvp_mat_unif, 1, GL_FALSE, glm::value_ptr(mvp));

  glEnableVertexAttribArray(_vertex_pos_attr);
  glEnableVertexAttribArray(_texture_coord_attr);
  glEnableVertexAttribArray(_lightmap_coord_attr);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  for (size_t i = 0; i < _faces.size(); i++) {
    if (!_visible_faces[i])
      continue;
    if (_faces[i].type == (int)face::polygon
        || _faces[i].type == (int)face::mesh) {
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, texture_ids[_faces[i].texture]);
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, _lightmap_texture_ids[_faces[i].lm_index]);
      glVertexAttribPointer(_vertex_pos_attr, 3, GL_FLOAT, GL_FALSE
          , sizeof(bsp_vertex), &_vertices[_faces[i].vertex].position);
      glVertexAttribPointer(_texture_coord_attr, 2, GL_FLOAT, GL_FALSE
          , sizeof(bsp_vertex), &_vertices[_faces[i].vertex].decal);
      glVertexAttribPointer(_lightmap_coord_attr, 2, GL_FLOAT, GL_FALSE
          , sizeof(bsp_vertex), &_vertices[_faces[i].vertex].lightmap);
      glDrawElements(GL_TRIANGLES, _faces[i].n_meshverts, GL_UNSIGNED_INT
          , &_meshverts[_faces[i].meshvert].offset);
    } else if (_faces[i].type == (int)face::patch) {
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, texture_ids[_patches[i]->texture_idx]);
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D
          , _lightmap_texture_ids[_patches[i]->lightmap_idx]);
      for (const bsp_biquadratic_patch &p : _patches[i]->quadratic_patches) {
        glVertexAttribPointer(_vertex_pos_attr, 3, GL_FLOAT, GL_FALSE
            , sizeof(bsp_vertex)
            , &p.vertices[0].position);
        glVertexAttribPointer(_texture_coord_attr, 2, GL_FLOAT, GL_FALSE
            , sizeof(bsp_vertex)
            , &p.vertices[0].decal);
        glVertexAttribPointer(_lightmap_coord_attr, 2, GL_FLOAT, GL_FALSE
            , sizeof(bsp_vertex)
            , &p.vertices[0].lightmap);
        p.render();
      }
    }
  }

  glDisableVertexAttribArray(_vertex_pos_attr);
  glDisableVertexAttribArray(_texture_coord_attr);
  glDisableVertexAttribArray(_lightmap_coord_attr);

  sp.dont_use_this_prog();
}

