#pragma once

// #include "ogl.hh"
// #include <vector>
#include <glm/glm.hpp>
#include <vector>

#if 0
class bitset {
  int num_bytes;
  unsigned char *bits;
public:
  bitset();
  ~bitset();
  void init(int number_of_bits);
  void clear_all();
  void set_all();
  void clear(int bit);
  void set(int bit);
  unsigned char is_set(int bit);
};

struct bsp_biquadratic_patch {
  bsp_vertex control_points[9];
  int subdivisions;
  bsp_vertex *vertices;
  unsigned int *indices;
  int *triangles_per_row;
  unsigned int **row_index_pointers;

  bsp_biquadratic_patch();
  ~bsp_biquadratic_patch();
  void tesselate(int subdivisions);
  void draw();
};
#endif

struct bsp_plane {
  glm::vec3 normal;
  float dist;

  bsp_plane();
#if 0
  bsp_plane(glm::vec3 newNormal, float n_intercept);
  bsp_plane(const bsp_plane &rhs);
  void set_normal(const glm::vec3 &rhs);
  void set_intercept(float n_intercept);
  void set_from_points(const glm::vec3 &p0, const glm::vec3 &p1
      , const glm::vec3 &p2);
  void calculate_intercept(const glm::vec3 &point_on_plane);
  void normalize(void);
  bool intersect3(const bsp_plane &p2, const bsp_plane &p3, glm::vec3 &result);
  float get_distance(const glm::vec3 &point) const;
  int classify_point(const glm::vec3 &point) const;
  bsp_plane lerp(const bsp_plane &p2, float factor);
  bool operator==(const bsp_plane &rhs) const;
  bool operator!=(const bsp_plane &rhs) const;
  bsp_plane operator-() const;
  bsp_plane operator+() const;
#endif
};

struct bsp_node {
  int plane;
  int front;
  int back;
  int mins[3];
  int maxs[3];
};

struct bsp_leaf {
  int cluster;
  int area;
  int mins[3];
  int maxs[3];
  int leafface;
  int n_leaffaces;
  int leafbrush;
  int n_leafbrushes;
};

struct bsp_leafface {
  int face;
};

struct bsp_vertex {
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

  bsp_vertex operator+(const bsp_vertex &v) const;
  bsp_vertex operator*(float factor) const;
};

struct bsp_face {
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

struct bsp_lightmap {
  char map[128][128][3];
};

struct bsp_visdata {
  int n_vecs;
  int sz_vecs;
  std::vector<unsigned char> vecs;
};

class bsp {
  std::vector<bsp_plane> _planes;
  std::vector<bsp_node> _nodes;
  std::vector<bsp_leaf> _leaves;
  std::vector<bsp_leafface> _leaffaces;
  std::vector<bsp_vertex> _vertices;
  std::vector<bsp_face> _faces;
  std::vector<int> _visible_faces;
  std::vector<bsp_lightmap> _lightmaps;
  std::vector<unsigned int> _lightmap_texture_ids;
  bsp_visdata _visdata;
public:
  bsp(const char *filename);
  int find_leaf(glm::vec3 position);
  int cluster_visible(int vis_cluster, int test_cluster);
  void set_visible_faces(glm::vec3 camera_pos);
};

