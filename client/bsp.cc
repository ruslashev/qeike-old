#include "bsp.hh"
#include "utils.hh"
#include "shaders.hh"
#include <fstream>
#include <glm/gtc/type_ptr.hpp>

const float feps = 1e-4f;
const float scale = 16.f;

bsp_plane::bsp_plane()
  : normal({0.f, 0.f, 0.f})
  , dist(0.0f) {
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

bsp_vertex bsp_vertex::operator+(const bsp_vertex &v) const {
  bsp_vertex res;
  res.position_x = position_x + v.position_x;
  res.position_y = position_y + v.position_y;
  res.position_z = position_z + v.position_z;
  res.decal_s = decal_s + v.decal_s;
  res.decal_t = decal_t + v.decal_t;
  res.lightmap_s = lightmap_s + v.lightmap_s;
  res.lightmap_t = lightmap_t + v.lightmap_t;
  res.normal_x = normal_x + v.normal_x;
  res.normal_y = normal_y + v.normal_y;
  res.normal_z = normal_z + v.normal_z;
  return res;
}

bsp_vertex bsp_vertex::operator*(float factor) const {
  bsp_vertex res;
  res.position_x = position_x * factor;
  res.position_y = position_y * factor;
  res.position_z = position_z * factor;
  res.decal_s = decal_s * factor;
  res.decal_t = decal_t * factor;
  res.lightmap_s = lightmap_s * factor;
  res.lightmap_t = lightmap_t * factor;
  res.normal_x = normal_x * factor;
  res.normal_y = normal_y * factor;
  res.normal_z = normal_z * factor;
  return res;
}
#endif

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

void bsp::_load_file(const char *filename) {
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
    p.dist /= scale;
  }

  load_lump(ifs, &header, lump::nodes, _nodes);
  for (bsp_node &n : _nodes) {
    n.mins /= scale;
    std::swap(n.mins.y, n.mins.z);
    n.mins.z = -n.mins.z;
    n.maxs /= scale;
    std::swap(n.maxs.y, n.maxs.z);
    n.maxs.z = -n.maxs.z;
  }

  load_lump(ifs, &header, lump::leaves, _leaves);
  for (bsp_leaf &l : _leaves) {
    l.mins /= scale;
    std::swap(l.mins.y, l.mins.z);
    l.mins.z = -l.mins.z;
    l.maxs /= scale;
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
    v.position /= scale;
    std::swap(v.position.y, v.position.z);
    v.position.z = -v.position.z;
  }

  load_lump(ifs, &header, lump::meshverts, _meshverts);

  load_lump(ifs, &header, lump::faces, _faces);
  _visible_faces.resize(_faces.size(), 0);

  load_lump(ifs, &header, lump::lightmaps, _lightmaps);
  _lightmap_texture_ids.resize(_lightmaps.size());

  ifs.seekg(header.direntries[(int)lump::visdata].offset);
  ifs.read((char*)&_visdata.n_vecs, sizeof(int));
  ifs.read((char*)&_visdata.sz_vecs, sizeof(int));
  int size = _visdata.n_vecs * _visdata.sz_vecs;
  _visdata.vecs.reserve(size);
  ifs.read((char*)&_visdata.vecs[0], size * sizeof(unsigned char));

  ifs.close();
}

bsp::bsp(const char *filename)
  : sp(shaders::map_vert, shaders::map_frag) {
  _load_file(filename);

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

int bsp::find_leaf(glm::vec3 position) {
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

int bsp::cluster_visible(int vis_cluster, int test_cluster) {
  if (vis_cluster < 0 || _visdata.vecs.size() == 0)
    return 1;
  int i = vis_cluster * _visdata.sz_vecs + (test_cluster >> 3);
  unsigned char vis_set = _visdata.vecs[i];
  return vis_set & (1 << (test_cluster & 7));
}

void bsp::set_visible_faces(glm::vec3 camera_pos) {
  int leaf_index = find_leaf(camera_pos);
  std::fill(_visible_faces.begin(), _visible_faces.end(), 0);
  for (bsp_leaf &l : _leaves)
    if (cluster_visible(_leaves[leaf_index].cluster, l.cluster))
      for (int j = 0; j < l.n_leaffaces; j++)
        _visible_faces[_leaffaces[l.leafface + j].face] = 1;
}

void bsp::render(const glm::mat4 &mvp) {
  sp.use_this_prog();

  glUniformMatrix4fv(_mvp_mat_unif, 1, GL_FALSE, glm::value_ptr(mvp));

  glEnableVertexAttribArray(_vertex_pos_attr);
  glEnableVertexAttribArray(_texture_coord_attr);
  glEnableVertexAttribArray(_lightmap_coord_attr);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  for (size_t i = 0; i < _faces.size(); i++) {
    if (!_visible_faces[i] || (_faces[i].type != (int)face::polygon
          && _faces[i].type != (int)face::mesh))
      continue;
    glVertexAttribPointer(_vertex_pos_attr, 3, GL_FLOAT, GL_FALSE
        , sizeof(bsp_vertex), &_vertices[_faces[i].vertex].position);
    glVertexAttribPointer(_texture_coord_attr, 2, GL_FLOAT, GL_FALSE
        , sizeof(bsp_vertex), &_vertices[_faces[i].vertex].decal);
    glVertexAttribPointer(_lightmap_coord_attr, 2, GL_FLOAT, GL_FALSE
        , sizeof(bsp_vertex), &_vertices[_faces[i].vertex].lightmap);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texture_ids[_faces[i].texture]);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, _lightmap_texture_ids[_faces[i].lm_index]);
    glDrawElements(GL_TRIANGLES, _faces[i].n_meshverts, GL_UNSIGNED_INT
        , &_meshverts[_faces[i].meshvert].offset);
  }

  glDisableVertexAttribArray(_vertex_pos_attr);
  glDisableVertexAttribArray(_texture_coord_attr);
  glDisableVertexAttribArray(_lightmap_coord_attr);

  sp.dont_use_this_prog();
}

