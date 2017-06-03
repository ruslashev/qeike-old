#include "bsp.hh"
#include "utils.hh"
#include <cstring> // memset
#include <fstream>

const float feps = 1e-4f;
const float scale = 16.f;

#if 0
bsp_biquadratic_patch::bsp_biquadratic_patch()
  : vertices(nullptr)
  , indices(nullptr) {
}

bsp_biquadratic_patch::~bsp_biquadratic_patch() {
  if (vertices)
    delete [] vertices;
  if (indices)
    delete [] indices;
}

void bsp_biquadratic_patch::tesselate(int n_subdivisions) {
  subdivisions = n_subdivisions;
  float px, py;
  bsp_vertex temp[3];
  vertices = new bsp_vertex [(subdivisions + 1) * (subdivisions + 1)];
  for (int v = 0; v <= subdivisions; ++v) {
    px = (float)v / subdivisions;
    vertices[v] = control_points[0] * ((1.0f - px) * (1.0f - px))
      + control_points[3] * ((1.0f - px) * px * 2)
      + control_points[6] * (px * px);
  }
  for (int u = 1; u <= subdivisions; ++u) {
    py = (float)u / subdivisions;
    temp[0] = control_points[0] * (1.0f - py) * (1.0f - py)
      + control_points[1] * (1.0f - py) * py * 2
      + control_points[2] * py * py;
    temp[1] = control_points[3] * (1.0f - py) * (1.0f - py)
      + control_points[4] * (1.0f - py) * py * 2
      + control_points[5] * py * py;
    temp[2] = control_points[6] * (1.0f - py) * (1.0f - py)
      + control_points[7] * (1.0f - py) * py * 2
      + control_points[8] * py * py;
    for (int v = 0; v <= subdivisions; ++v) {
      px = (float)v / subdivisions;
      vertices[u * (subdivisions + 1) + v] = temp[0] * (1.0f - px) * (1.0f - px)
        + temp[1] * (1.0f - px) * px * 2 + temp[2] * px * px;
    }
  }
  indices = new unsigned int [subdivisions * (subdivisions + 1) * 2];
  for (int row = 0; row < subdivisions; ++row)
    for (int point = 0; point <= subdivisions; ++point) {
      // reverse winding
      indices[(row * (subdivisions + 1) + point) * 2 + 1]
        = row * (subdivisions + 1) + point;
      indices[(row * (subdivisions + 1) + point) * 2]
        = (row + 1) * (subdivisions + 1) + point;
    }

  triangles_per_row = new int [subdivisions];
  row_index_pointers = new unsigned int* [subdivisions];

  for (int row = 0; row < subdivisions; ++row) {
    triangles_per_row[row] = 2 * (subdivisions + 1);
    row_index_pointers[row] = &indices[row * 2 * (subdivisions + 1)];
  }
}
#endif

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

bsp::bsp(const char *filename) {
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

  load_lump(ifs, &header, lump::vertices, vertices);
  for (bsp_vertex &v : vertices) {
    v.position /= scale;
    std::swap(v.position.y, v.position.z);
    v.position.z = -v.position.z;
  }

  load_lump(ifs, &header, lump::meshverts, meshverts);

  load_lump(ifs, &header, lump::faces, faces);
  visible_faces.resize(faces.size(), 0);

  load_lump(ifs, &header, lump::lightmaps, _lightmaps);
  lightmap_texture_ids.resize(_lightmaps.size());

  ifs.seekg(header.direntries[(int)lump::visdata].offset);
  ifs.read((char*)&_visdata.n_vecs, sizeof(int));
  ifs.read((char*)&_visdata.sz_vecs, sizeof(int));
  int size = _visdata.n_vecs * _visdata.sz_vecs;
  _visdata.vecs.reserve(size);
  ifs.read((char*)&_visdata.vecs[0], size * sizeof(unsigned char));

#if 0
  _patches = new bsp_patch [_num_patches];
  for (int curr_patch = 0, i = 0; i < num_faces; ++i) {
    if (read_faces[i].type != patch)
      continue;
    _patches[curr_patch].texture_idx = read_faces[i].texture;
    _patches[curr_patch].lightmap_idx = read_faces[i].lm_index;
    _patches[curr_patch].width = read_faces[i].size[0];
    _patches[curr_patch].height = read_faces[i].size[1];
    face_directory[i].face_type = patch;
    face_directory[i].face_number = curr_patch;

    int num_patches_width = (_patches[curr_patch].width - 1) / 2
      , num_patches_height = (_patches[curr_patch].height - 1) / 2;

    _patches[curr_patch].num_quadratic_patches = num_patches_width
      * num_patches_height;
    _patches[curr_patch].quadratic_patches =
      new bsp_biquadratic_patch [_patches[curr_patch].num_quadratic_patches];
    for (int y = 0; y < num_patches_height; ++y)
      for (int x = 0; x < num_patches_width; ++x) {
        for (int row = 0; row < 3; ++row)
          for (int col = 0; col < 3; ++col)
            _patches[curr_patch].quadratic_patches[y
              * num_patches_width+x].control_points[row * 3 + col]
              = vertices[read_faces[i].vertex
              + (y * 2 * _patches[curr_patch].width + x * 2)
              + row * _patches[curr_patch].width + col];
        _patches[curr_patch].quadratic_patches[y * num_patches_width + x]
          .tesselate(10);
      }
    ++curr_patch;
  }
  delete [] read_faces;
#endif

  ifs.close();

  vbo.bind();
  vbo.upload(sizeof(vertices[0]) * vertices.size(), &vertices[0]);

  glActiveTexture(GL_TEXTURE0);
  for (size_t i = 0; i < _textures.size(); ++i)
    texture_ids[i] = 0;

  glActiveTexture(GL_TEXTURE1);
  glGenTextures(_lightmaps.size(), &lightmap_texture_ids[0]);
  for (size_t i = 0; i < _lightmaps.size(); ++i) {
    glBindTexture(GL_TEXTURE_2D, lightmap_texture_ids[i]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 128, 128, 0, GL_RGB
        , GL_UNSIGNED_BYTE, _lightmaps[i].map);
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }
}

bsp::~bsp() {
  glDeleteTextures(_textures.size(), texture_ids.data());
  glDeleteTextures(_lightmaps.size(), lightmap_texture_ids.data());
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
  std::fill(visible_faces.begin(), visible_faces.end(), 0);
  for (bsp_leaf &l : _leaves)
    if (cluster_visible(_leaves[leaf_index].cluster, l.cluster))
      for (int j = 0; j < l.n_leaffaces; j++)
        visible_faces[_leaffaces[l.leafface + j].face] = 1;
}

