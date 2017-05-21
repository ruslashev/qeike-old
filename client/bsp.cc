#if 0

#include <GL/glew.h>
#include <map>
#include <algorithm>

enum face_type { POLYGON = 1, PATCH, MESH, BILLBOARD };

#include "utils.hh"
#include <fstream>

struct bsp_entities {
  char *ents;
};

struct bsp_texture {
  char name[64];
  int flags;
  int contents;
};

struct bsp_leafbrush {
  int brush;
};

struct bsp_model {
  float mins[3];
  float maxs[3];
  int face;
  int num_faces;
  int brush;
  int num_brushes;
};

struct bsp_brush {
  int brushside;
  int num_brushsides;
  int texture;
};

struct bsp_brushside {
  int plane;
  int texture;
};


struct bsp_effect {
  char name[64];
  int brush;
  int unknown;
};

struct bsp_lightmap {
  unsigned char map[128][128][3];
};

struct bsp_lightvol {
  unsigned char ambient[3];
  unsigned char directional[3];
  unsigned char dir[2];
};

class bezier {
public:
  bezier();
  ~bezier();

  bsp_vertex calculate_quadratic_bezier(float t
      , const bsp_vertex* control_vertexes);
  void tessellate(int subdivisions);

  bsp_vertex* m_vertexes;
  bsp_vertex m_control_vertexes[9];
  unsigned int* m_indexes;
  unsigned int* m_row_indexes[10];
  unsigned int m_tri_per_row[10];

  unsigned int m_vertex_offset;
  unsigned int m_index_offset;
};

bezier::bezier() : m_vertexes(nullptr), m_indexes(nullptr) {
}

bezier::~bezier() {
  if (m_vertexes)
    delete[] m_vertexes;
  if (m_indexes)
    delete[] m_indexes;
}


// TODO: for adding color interpolation maybe check for t < 0.5 and use
// 1+2 or 2+3 vertex accordingly
bsp_vertex bezier::calculate_quadratic_bezier(float t
    , const bsp_vertex* control_vertexes) {
  return (control_vertexes[0]) * (1 - t) * (1 - t) + (control_vertexes[1]) * 2
    * t * (1 - t) + (control_vertexes[2]) * t * t;
}

void bezier::tessellate(int subdivisions) {
  bsp_vertex temp[3];
  int subdivisions1 = subdivisions + 1;

  if (m_vertexes != nullptr)
    delete [] m_vertexes;
  m_vertexes = new bsp_vertex[subdivisions1 * subdivisions1];

  for (int i = 0; i <= subdivisions; ++i) {
    float l = (float)i/subdivisions;

    for (int j = 0; j < 3; ++j) {
      int k = 3 * j;
      temp[j] = calculate_quadratic_bezier(l, &(m_control_vertexes[k]));
    }

    int col = 0;
    for(int j = 0; j <= subdivisions; ++j) {
      float a = (float)j / subdivisions;
      m_vertexes[i * subdivisions1 + j] = calculate_quadratic_bezier(a, temp);
      m_vertexes[i * subdivisions1 + j].color[0] = 0xff;
      m_vertexes[i * subdivisions1 + j].color[1] = 0xff;
      m_vertexes[i * subdivisions1 + j].color[2] = 0xff;
      m_vertexes[i * subdivisions1 + j].color[3] = 0xff;
    }
  }

  if (m_indexes != NULL)
    delete [] m_indexes;
  m_indexes = new unsigned int[subdivisions * subdivisions1 * 2];

  // maybe use degenerated triangle strips to merge
  for (int row = 0; row < subdivisions; ++row) {
    for(int col = 0; col <= subdivisions; ++col) {
      int g = (row * (subdivisions1) + col) * 2 + 1;
      int h = (row * (subdivisions1) + col) * 2;
      m_indexes[g] = row       * subdivisions1 + col;
      m_indexes[h] = (row + 1) * subdivisions1 + col;
    }
  }

  for (int row = 0; row < subdivisions; ++row) {
    m_tri_per_row[row] = 2 * subdivisions1;
    m_row_indexes[row] = &(m_indexes[row * 2 * subdivisions1]);
  }
}

_bsp::_bsp(const char *filename) {
  std::ifstream ifs(filename, std::ios::binary);

  bsp_header header;
  ifs.read((char*)&header, sizeof(bsp_header));
  std::string magic = std::string(header.magic, 4);
  assertf(magic == "IBSP", "_bsp \"%s\": invalid magic \"%s\"", filename
      , magic.c_str());
  assertf(header.version == 0x2f, "_bsp \"%s\": invalid version %d = 0x%x"
      , filename, header.version, header.version);

  int num_textures = header.direntries[LUMP_TEXTURES].length
    / sizeof(bsp_texture);
  int num_planes = header.direntries[LUMP_PLANES].length / sizeof(bsp_plane);
  int num_nodes = header.direntries[LUMP_NODES].length / sizeof(bsp_node);
  _num_leafs = header.direntries[LUMP_LEAFS].length / sizeof(bsp_leaf);
  int num_leaffaces =
      header.direntries[LUMP_LEAFFACES].length / sizeof(bsp_leafface);
  int num_leafbrushes =
      header.direntries[LUMP_LEAFBRUSHES].length / sizeof(bsp_leafbrush);
  int num_models = header.direntries[LUMP_MODELS].length / sizeof(bsp_model);
  int num_brushes = header.direntries[LUMP_BRUSHES].length / sizeof(bsp_brush);
  int num_brushsides =
      header.direntries[LUMP_BRUSHSIDES].length / sizeof(bsp_brushside);
  int num_vertexes = header.direntries[LUMP_VERTICES].length
    / sizeof(bsp_vertex);
  int num_meshverts =
      header.direntries[LUMP_MESHVERTS].length / sizeof(bsp_meshvert);
  int num_effects = header.direntries[LUMP_EFFECTS].length / sizeof(bsp_effect);
  int num_faces = header.direntries[LUMP_FACES].length / sizeof(bsp_face);
  int num_lightmaps =
      header.direntries[LUMP_LIGHTMAPS].length / sizeof(bsp_lightmap);
  int num_lightvols =
      header.direntries[LUMP_LIGHTVOLS].length / sizeof(bsp_lightvol);

  bsp_entities *entities = new bsp_entities;
  entities->ents = new char[header.direntries[LUMP_ENTITIES].length];
  bsp_texture *textures = new bsp_texture[num_textures];
  _nodes = new bsp_node[num_nodes];
  _planes = new bsp_plane[num_planes];
  _leafs = new bsp_leaf[_num_leafs];
  _leaffaces = new bsp_leafface[num_leaffaces];
  bsp_leafbrush *leafbrushes = new bsp_leafbrush[num_leafbrushes];
  bsp_model *models = new bsp_model[num_models];
  bsp_brush *brushes = new bsp_brush[num_brushes];
  bsp_brushside *brushsides = new bsp_brushside[num_brushsides];
  bsp_vertex *vertexes = new bsp_vertex[num_vertexes];
  bsp_meshvert *meshverts = new bsp_meshvert[num_meshverts];
  bsp_effect *effects = new bsp_effect[num_effects];
  _faces = new bsp_face[num_faces];
  bsp_lightmap *lightmaps = new bsp_lightmap[num_lightmaps];
  bsp_lightvol *lightvols = new bsp_lightvol[num_lightvols];
  _visdata = new bsp_visdata;

  ifs.seekg(header.direntries[LUMP_ENTITIES].offset);
  ifs.read((char*)entities->ents, header.direntries[LUMP_ENTITIES].length);

  ifs.seekg(header.direntries[LUMP_TEXTURES].offset);
  ifs.read((char*)textures, header.direntries[LUMP_TEXTURES].length);

  ifs.seekg(header.direntries[LUMP_PLANES].offset);
  ifs.read((char*)_planes, header.direntries[LUMP_PLANES].length);

  ifs.seekg(header.direntries[LUMP_NODES].offset);
  ifs.read((char*)_nodes, header.direntries[LUMP_NODES].length);

  ifs.seekg(header.direntries[LUMP_LEAFS].offset);
  ifs.read((char*)_leafs, header.direntries[LUMP_LEAFS].length);

  ifs.seekg(header.direntries[LUMP_LEAFFACES].offset);
  ifs.read((char*)_leaffaces, header.direntries[LUMP_LEAFFACES].length);

  ifs.seekg(header.direntries[LUMP_LEAFBRUSHES].offset);
  ifs.read((char*)leafbrushes, header.direntries[LUMP_LEAFBRUSHES].length);

  ifs.seekg(header.direntries[LUMP_MODELS].offset);
  ifs.read((char*)models, header.direntries[LUMP_MODELS].length);

  ifs.seekg(header.direntries[LUMP_BRUSHES].offset);
  ifs.read((char*)brushes, header.direntries[LUMP_BRUSHES].length);

  ifs.seekg(header.direntries[LUMP_BRUSHSIDES].offset);
  ifs.read((char*)brushsides, header.direntries[LUMP_BRUSHSIDES].length);

  ifs.seekg(header.direntries[LUMP_VERTICES].offset);
  ifs.read((char*)vertexes, header.direntries[LUMP_VERTICES].length);

  ifs.seekg(header.direntries[LUMP_MESHVERTS].offset);
  ifs.read((char*)meshverts, header.direntries[LUMP_MESHVERTS].length);

  ifs.seekg(header.direntries[LUMP_EFFECTS].offset);
  ifs.read((char*)effects, header.direntries[LUMP_EFFECTS].length);

  ifs.seekg(header.direntries[LUMP_FACES].offset);
  ifs.read((char*)_faces, header.direntries[LUMP_FACES].length);

  ifs.seekg(header.direntries[LUMP_LIGHTMAPS].offset);
  ifs.read((char*)lightmaps, header.direntries[LUMP_LIGHTMAPS].length);

  ifs.seekg(header.direntries[LUMP_LIGHTVOLS].offset);
  ifs.read((char*)lightvols, header.direntries[LUMP_LIGHTVOLS].length);

  ifs.seekg(header.direntries[LUMP_VISDATA].offset);
  ifs.read((char*)_visdata, 2 * sizeof(int));
  _visdata->vecs = new unsigned char[_visdata->n_vecs * _visdata->sz_vecs];
  ifs.read((char*)_visdata->vecs,
           sizeof(unsigned char) * _visdata->n_vecs * _visdata->sz_vecs);

  ifs.close();

  std::map<bsp_face*, std::vector<bezier*>> patches;

  for (int i = 0; i < num_faces; ++i)
    if (_faces[i].type == PATCH) {
      bsp_face *face = &(_faces[i]);

      int width = _faces[i].size[0];
      int height = _faces[i].size[1];
      int widthCount = (width - 1) / 2;
      int heightCount = (height - 1) / 2;

      patches[face].resize(widthCount * heightCount);
      for (int j = 0; j < widthCount * heightCount; ++j)
        patches[face][j] = new bezier();

      for (int y = 0; y < heightCount; y++)
        for (int x = 0; x < widthCount; x++) {
          for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++) {
              patches[face][y * widthCount + x]
                  ->m_control_vertexes[row * 3 + col] =
                  vertexes[_faces[i].vertex + (y * 2 * width + x * 2) +
                            row * width + col];
            }
          patches[face][y * widthCount + x]->tessellate(10);
        }
    }

  int num_bezier_vertexes = 0, num_bezier_indexes = 0;

  for (auto it = patches.begin(); it != patches.end(); ++it) {
    num_bezier_vertexes += 11 * 11 * it->second.size();
    num_bezier_indexes += 10 * 11 * 2 * it->second.size();
  }

  glGenBuffers(1, &_vbo_id);
  glBindBuffer(GL_ARRAY_BUFFER, _vbo_id);
  glBufferData(GL_ARRAY_BUFFER, header.direntries[LUMP_VERTICES].length +
      num_bezier_vertexes * sizeof(bsp_vertex), NULL, GL_STATIC_DRAW);
  glBufferSubData(GL_ARRAY_BUFFER, 0, header.direntries[LUMP_VERTICES].length,
      vertexes);

  glGenBuffers(1, &_ebo_id);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _ebo_id);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, header.direntries[LUMP_MESHVERTS].length
      + num_bezier_indexes * sizeof(unsigned int), NULL, GL_STATIC_DRAW);
  glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0,
      header.direntries[LUMP_MESHVERTS].length, meshverts);

  int offset_verts = 0, offset_idx = 0;
  for (auto it = patches.begin(); it != patches.end(); ++it)
    for (unsigned int j = 0; j < it->second.size(); ++j) {
      glBufferSubData(GL_ARRAY_BUFFER, header.direntries[LUMP_VERTICES].length
          + offset_verts, 11 * 11 * sizeof(bsp_vertex)
          , it->second[j]->m_vertexes);

      glBufferSubData(GL_ELEMENT_ARRAY_BUFFER,
          header.direntries[LUMP_MESHVERTS].length + offset_idx, 10 * 11 * 2
          * sizeof(unsigned int), it->second[j]->m_indexes);

      it->second[j]->m_vertex_offset = header.direntries[LUMP_VERTICES].length
        + offset_verts;
      it->second[j]->m_index_offset = header.direntries[LUMP_MESHVERTS].length
        + offset_idx;

      offset_verts += sizeof(bsp_vertex) * 11 * 11;
      offset_idx += sizeof(unsigned int) * 10 * 11 * 2;
    }

  delete [] entities->ents;
  delete entities;
  delete [] textures;
  delete [] leafbrushes;
  delete [] models;
  delete [] brushes;
  delete [] brushsides;
  delete [] vertexes;
  delete [] meshverts;
  delete [] effects;
  delete [] lightmaps;
  delete [] lightvols;
}

_bsp::~_bsp() {
  delete [] _leafs;
  delete [] _leaffaces;
  delete [] _nodes;
  delete [] _planes;
  delete [] _faces;
  delete _visdata;
}

int _bsp::find_leaf(const glm::vec4 &camera_position) {
  int index = 0;

  while (index >= 0) {
    const bsp_node &node = _nodes[index];
    const bsp_plane &plane = _planes[node.plane];

    // TODO: multiplicate the plane with our transformation matrices (inverse
    // and transpose needed for plane transformation)
    // TODO: maybe do this during load time to enhance framerate!

    glm::vec4 pos = camera_position;

    /*  if (plane.type < 3) // type < 3 -> axial plane
    {
    const float distance = pos[plane->type] - plane.distance;
    }
    else
    */
    const float distance =
        glm::dot(plane.normal, glm::vec3(pos)) - plane.dist;

    if (distance >= 0)
      index = node.front;
    else
      index = node.back;
  }

  return -index - 1;
}

bool _bsp::cluster_is_visible(int cluster, int test_cluster) {
  if ((_visdata->vecs == NULL) || (cluster < 0))
    return true;

  int i = (cluster * _visdata->sz_vecs) + (test_cluster >> 3);
  unsigned char visSet = _visdata->vecs[i];

  if (!(visSet & (1 << (test_cluster & 7))))
    return false;
  return true;
}

std::vector<bsp_face*>
_bsp::compute_visible_faces(const glm::vec4 &camera_position) {
  std::bitset<10000> already_visible;
  std::vector<bsp_face *> visible_faces;

  int leafindex = find_leaf(camera_position);
  int cluster = _leafs[leafindex].cluster;

  for (int i = _num_leafs - 1; i >= 0; --i) {
    if (!cluster_is_visible(cluster, _leafs[i].cluster))
      continue;

    for (int j = _leafs[i].leafface + _leafs[i].n_leaffaces - 1;
         j >= _leafs[i].leafface; --j) {
      int face = _leaffaces[j].face;
      if (already_visible.test(face))
        continue;
      already_visible.set(face);

      visible_faces.push_back(&(_faces[face]));
    }
  }

  std::sort(visible_faces.begin(), visible_faces.end()
      , [](const bsp_face *left, const bsp_face *right) {
        if (left->texture == right->texture) {
          if (left->lm_index < right->lm_index)
            return true;
          return false;
        } else if (left->texture < right->texture)
          return true;
        return false;
      });

  return visible_faces;
}

struct bsp_node {
  int plane;
  int front;
  int back;
  int mins[3];
  int maxs[3];
};

struct bsp_plane {
  glm::vec3 normal;
  float dist;
};

struct bsp_visdata {
  int n_vecs;
  int sz_vecs;
  unsigned char *vecs; // [n_vecs * sz_vecs]
};

#endif

//////////////////////////////

#include "bsp.hh"
#include "utils.hh"
#include <cstring> // memset
#include <fstream>

bitset::bitset() : num_bytes(0), bits(nullptr) {
}

bitset::~bitset() {
  if (bits)
    delete [] bits;
}

void bitset::init(int number_of_bits) {
  if (bits)
    delete [] bits;
  bits = nullptr;
  num_bytes = (number_of_bits >> 3) + 1;
  bits = new unsigned char [num_bytes];
  clear_all();
}

void bitset::clear_all() {
  memset(bits, 0, num_bytes);
}

void bitset::set_all() {
  memset(bits, 0xFF, num_bytes);
}

void bitset::clear(int bit) {
  bits[bit >> 3] &= ~(1 << (bit & 7));
}

void bitset::set(int bit) {
  bits[bit >> 3] |= 1 << (bit & 7);
}

unsigned char bitset::is_set(int bit) {
  return bits[bit >> 3] & 1 << (bit & 7);
}

bsp_plane::bsp_plane()
  : normal(glm::vec3(0))
  , intercept(0.0f) {
}

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

enum lump {
  LUMP_ENTITIES = 0,
  LUMP_TEXTURES,
  LUMP_PLANES,
  LUMP_NODES,
  LUMP_LEAVES,
  LUMP_LEAFFACES,
  LUMP_LEAFBRUSHES,
  LUMP_MODELS,
  LUMP_BRUSHES,
  LUMP_BRUSHSIDES,
  LUMP_VERTICES,
  LUMP_MESHVERTS,
  LUMP_EFFECTS,
  LUMP_FACES,
  LUMP_LIGHTMAPS,
  LUMP_LIGHTVOLS,
  LUMP_VISDATA
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

struct bsp_read_vertex {
  float position_x;
  float position_y;
  float position_z;
  float decal_s;
  float decal_t;
  float lightmap_s;
  float lightmap_t;
  float normal_x;
  float normal_y;
  float normal_z;
  unsigned char color[4];
}; // TODO merge bsp_read_vertex and bsp_vertex

struct bsp_meshvert {
  int offset;
};

struct bsp_read_face {
  int texture;
  int effect;
  int type;
  int vertex;
  int n_vertices;
  int meshvert;
  int n_meshverts;
  int lm_index;
  int lm_start[2];
  int lm_size[2];
  float lm_origin_x;
  float lm_origin_y;
  float lm_origin_z;
  float lm_s_x;
  float lm_s_y;
  float lm_s_z;
  float lm_t_x;
  float lm_t_y;
  float lm_t_z;
  float normal_x;
  float normal_y;
  float normal_z;
  int size[2];
};

enum bsp_face_type {
  polygon = 1,
  patch,
  mesh,
  billboard
};

struct bsp_face_direntry {
  bsp_face_type face_type;
  int face_number;
};

struct bsp_read_leaf {
  int cluster;
  int area;
  int mins[3];
  int maxs[3];
  int leafface;
  int n_leaffaces;
  int leafbrush;
  int n_leafbrushes;
};

bsp_visdata::~bsp_visdata() {
  delete [] bitset;
}

bsp::bsp(const char *filename)
  : _num_polygon_faces(0)
  , _num_mesh_faces(0)
  , _num_patches(0)
  , _num_leaves(0)
  , _num_planes(0)
  , _num_nodes(0)
  , _faces(nullptr)
  , _mesh_faces(nullptr)
  , _patches(nullptr)
  , _leaves(nullptr)
  , _leaffaces(nullptr)
  , _planes(nullptr)
  , _nodes(nullptr)
  , _visdata(nullptr) {
  std::ifstream ifs(filename, std::ios::binary);

  bsp_header header;
  ifs.read((char*)&header, sizeof(bsp_header));
  std::string magic = std::string(header.magic, 4);
  // TODO: file stream does not close on failures
  assertf(magic == "IBSP", "\"%s\": invalid magic \"%s\"", filename
      , magic.c_str());
  assertf(header.version == 0x2f, "\"%s\": invalid version %d = 0x%x"
      " (expected 0x2f)", filename, header.version, header.version);

  int num_vertices = header.direntries[LUMP_VERTICES].length
    / sizeof(bsp_vertex);
  bsp_read_vertex *read_vertices = new bsp_read_vertex [num_vertices];
  ifs.seekg(header.direntries[LUMP_VERTICES].offset);
  ifs.read((char*)read_vertices, header.direntries[LUMP_VERTICES].length);
  bsp_vertex *vertices = new bsp_vertex [num_vertices];
  // convert bsp_read_vertex to bsp_vertex
  for (int i = 0; i < num_vertices; ++i) {
    // swap y and z and negate z
    vertices[i].position_x =  read_vertices[i].position_x;
    vertices[i].position_y =  read_vertices[i].position_z;
    vertices[i].position_z = -read_vertices[i].position_y;
    // scale down
    vertices[i].position_x /= 64.f;
    vertices[i].position_y /= 64.f;
    vertices[i].position_z /= 64.f;
    // invert texture coordinates
    vertices[i].decal_s =  read_vertices[i].decal_s;
    vertices[i].decal_t = -read_vertices[i].decal_t;
    // Transfer lightmap coordinates
    vertices[i].lightmap_s = read_vertices[i].lightmap_s;
    vertices[i].lightmap_t = read_vertices[i].lightmap_t;
  }
  delete [] read_vertices;

  int num_meshverts = header.direntries[LUMP_MESHVERTS].length
    / sizeof(bsp_meshvert);
  bsp_meshvert *meshverts = new bsp_meshvert [num_meshverts];
  ifs.seekg(header.direntries[LUMP_MESHVERTS].offset);
  ifs.read((char*)meshverts, header.direntries[LUMP_MESHVERTS].length);

  int num_faces = header.direntries[LUMP_FACES].length / sizeof(bsp_face);
  bsp_read_face *read_faces = new bsp_read_face [num_faces];
  ifs.seekg(header.direntries[LUMP_FACES].offset);
  ifs.read((char*)read_faces, header.direntries[LUMP_FACES].length);

  bsp_face_direntry *face_directory = new bsp_face_direntry [num_faces];
  memset(face_directory, 0, num_faces * sizeof(bsp_face_direntry));

  _faces_to_draw.init(num_faces);

  for (int i = 0; i < num_faces; ++i)
    switch (read_faces[i].type) {
      case polygon: ++_num_polygon_faces; break;
      case patch:   ++_num_patches; break;
      case mesh:    ++_num_mesh_faces; break;
      default: break;
    }

  _faces = new bsp_face [_num_polygon_faces];
  // convert bsp_read_face to bsp_face
  for (int curr_face = 0, i = 0; i < num_faces; ++i) {
    if (read_faces[i].type != polygon)
      continue;
    _faces[curr_face].first_vertex_idx = read_faces[i].vertex;
    _faces[curr_face].num_vertices = read_faces[i].n_vertices;
    _faces[curr_face].texture_idx = read_faces[i].texture;
    _faces[curr_face].lightmap_idx = read_faces[i].lm_index;
    face_directory[i].face_type = polygon;
    face_directory[i].face_number = curr_face;
    ++curr_face; // TODO move to for loop incr
  }

  _mesh_faces = new bsp_mesh_face [_num_mesh_faces];
  for (int curr_face = 0, i = 0; i < num_faces; ++i) {
    if (read_faces[i].type != mesh)
      continue;
    _mesh_faces[curr_face].first_vertex_idx = read_faces[i].vertex;
    _mesh_faces[curr_face].num_vertices = read_faces[i].n_vertices;
    _mesh_faces[curr_face].texture_idx = read_faces[i].texture;
    _mesh_faces[curr_face].lightmap_idx = read_faces[i].lm_index;
    _mesh_faces[curr_face].first_mesh_idx = read_faces[i].meshvert;
    _mesh_faces[curr_face].num_mesh_indices = read_faces[i].n_meshverts;
    face_directory[i].face_type = mesh;
    face_directory[i].face_number = curr_face;
    ++curr_face;
  }

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

  _num_leaves = header.direntries[LUMP_LEAVES].length / sizeof(bsp_read_leaf);
  bsp_read_leaf *read_leaves = new bsp_read_leaf [_num_leaves];
  ifs.seekg(header.direntries[LUMP_LEAVES].offset);
  ifs.read((char*)read_leaves, header.direntries[LUMP_LEAVES].length);
  _leaves = new bsp_leaf [_num_leaves];
  for (int i = 0; i < _num_leaves; ++i) {
    _leaves[i].cluster = read_leaves[i].cluster;
    _leaves[i].first_leaf_face = read_leaves[i].leafface;
    _leaves[i].num_faces = read_leaves[i].n_leaffaces;
    _leaves[i].bounding_box_vertices_x[0] =  (float)read_leaves[i].mins[0];
    _leaves[i].bounding_box_vertices_y[0] =  (float)read_leaves[i].mins[2];
    _leaves[i].bounding_box_vertices_z[0] = -(float)read_leaves[i].mins[1];
    _leaves[i].bounding_box_vertices_x[1] =  (float)read_leaves[i].mins[0];
    _leaves[i].bounding_box_vertices_y[1] =  (float)read_leaves[i].mins[2];
    _leaves[i].bounding_box_vertices_z[1] = -(float)read_leaves[i].maxs[1];
    _leaves[i].bounding_box_vertices_x[2] =  (float)read_leaves[i].mins[0];
    _leaves[i].bounding_box_vertices_y[2] =  (float)read_leaves[i].maxs[2];
    _leaves[i].bounding_box_vertices_z[2] = -(float)read_leaves[i].mins[1];
    _leaves[i].bounding_box_vertices_x[3] =  (float)read_leaves[i].mins[0];
    _leaves[i].bounding_box_vertices_y[3] =  (float)read_leaves[i].maxs[2];
    _leaves[i].bounding_box_vertices_z[3] = -(float)read_leaves[i].maxs[1];
    _leaves[i].bounding_box_vertices_x[4] =  (float)read_leaves[i].maxs[0];
    _leaves[i].bounding_box_vertices_y[4] =  (float)read_leaves[i].mins[2];
    _leaves[i].bounding_box_vertices_z[4] = -(float)read_leaves[i].mins[1];
    _leaves[i].bounding_box_vertices_x[5] =  (float)read_leaves[i].maxs[0];
    _leaves[i].bounding_box_vertices_y[5] =  (float)read_leaves[i].mins[2];
    _leaves[i].bounding_box_vertices_z[5] = -(float)read_leaves[i].maxs[1];
    _leaves[i].bounding_box_vertices_x[6] =  (float)read_leaves[i].maxs[0];
    _leaves[i].bounding_box_vertices_y[6] =  (float)read_leaves[i].maxs[2];
    _leaves[i].bounding_box_vertices_z[6] = -(float)read_leaves[i].mins[1];
    _leaves[i].bounding_box_vertices_x[7] =  (float)read_leaves[i].maxs[0];
    _leaves[i].bounding_box_vertices_y[7] =  (float)read_leaves[i].maxs[2];
    _leaves[i].bounding_box_vertices_z[7] = -(float)read_leaves[i].maxs[1];
    for (int j = 0; j < 8; ++j) {
      _leaves[i].bounding_box_vertices_x[j] /= 64.f;
      _leaves[i].bounding_box_vertices_y[j] /= 64.f;
      _leaves[i].bounding_box_vertices_z[j] /= 64.f;
    }
  }
  delete [] read_leaves;

  int num_leaffaces = header.direntries[LUMP_LEAFFACES].length
    / sizeof(bsp_leafface);
  _leaffaces = new bsp_leafface [num_leaffaces];
  ifs.seekg(header.direntries[LUMP_LEAFFACES].offset);
  ifs.read((char*)_leaffaces, header.direntries[LUMP_LEAFFACES].length);

  _num_planes = header.direntries[LUMP_PLANES].length / sizeof(bsp_plane);
  _planes = new bsp_plane [_num_planes];
  ifs.seekg(header.direntries[LUMP_PLANES].offset);
  ifs.read((char*)_planes, header.direntries[LUMP_PLANES].length);
  // reverse the intercept on the planes and convert to OpenGL coordinates
  for (int i = 0; i < _num_planes; ++i) {
    //swap y and z and negate z
    float temp = _planes[i].normal.y;
    _planes[i].normal.y = _planes[i].normal.z;
    _planes[i].normal.z = -temp;
    _planes[i].intercept = -_planes[i].intercept;
    _planes[i].intercept /= 64.f; //scale down
  }

  _num_nodes = header.direntries[LUMP_NODES].length / sizeof(bsp_node);
  _nodes = new bsp_node [_num_nodes];
  ifs.seekg(header.direntries[LUMP_NODES].offset);
  ifs.read((char*)_nodes, header.direntries[LUMP_NODES].length);

  _visdata = new bsp_visdata;
  ifs.seekg(header.direntries[LUMP_VISDATA].offset);
  ifs.read((char*)_visdata, 2 * sizeof(int));
  _visdata->bitset = new unsigned char[_visdata->n_clusters
    * _visdata->bytes_per_cluster];
  ifs.read((char*)_visdata->bitset, sizeof(unsigned char)
      * _visdata->n_clusters * _visdata->bytes_per_cluster);

  puts("ok");
  ifs.close();
}

bsp::~bsp() {
  if (_faces)
    delete [] _faces;
  if (_mesh_faces)
    delete [] _mesh_faces;
  if (_patches)
    delete [] _patches;
  if (_leaves)
    delete [] _leaves;
  if (_leaffaces)
    delete [] _leaffaces;
  if (_planes)
    delete [] _planes;
  if (_nodes)
    delete [] _nodes;
  if (_visdata)
    delete _visdata;
}

