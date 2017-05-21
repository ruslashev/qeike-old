#pragma once

// #include "ogl.hh"
// #include <vector>
#include <glm/glm.hpp>

const float feps = 1e-4f;

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

struct bsp_plane {
  bsp_plane();
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
  glm::vec3 normal;
  float intercept;
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

  bsp_vertex operator+(const bsp_vertex &v) const {
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

  bsp_vertex operator*(float factor) const {
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
};

struct bsp_face {
  int first_vertex_idx;
  int num_vertices;
  int texture_idx;
  int lightmap_idx;
};

struct bsp_mesh_face {
  int first_vertex_idx;
  int num_vertices;
  int texture_idx;
  int lightmap_idx;
  int first_mesh_idx;
  int num_mesh_indices;
};

struct bsp_biquadratic_patch {
  bsp_vertex control_points[9];
  int subdivisions;
  bsp_vertex *vertices;
  unsigned int *indices;
  int *triangles_per_row;
  unsigned int **row_index_pointers;

  bsp_biquadratic_patch() : vertices(nullptr), indices(nullptr) {
  }
  ~bsp_biquadratic_patch() {
    if (vertices)
      delete [] vertices;
    if (indices)
      delete [] indices;
  }
  void tesselate(int subdivisions);
  void draw();
};

struct bsp_patch {
  int texture_idx;
  int lightmap_idx;
  int width, height;
  int num_quadratic_patches;
  bsp_biquadratic_patch *quadratic_patches;
};

struct bsp_leaf {
  int cluster;
  float bounding_box_vertices_x[8];
  float bounding_box_vertices_y[8];
  float bounding_box_vertices_z[8];
  int first_leaf_face;
  int num_faces;
};

struct bsp_leafface {
  int face;
};

struct bsp_node {
  int plane;
  int front;
  int back;
  int mins[3];
  int maxs[3];
};

struct bsp_visdata {
  int n_clusters;
  int bytes_per_cluster;
  unsigned char *bitset;
  ~bsp_visdata();
};

class bsp {
  bitset _faces_to_draw;
  int _num_polygon_faces, _num_mesh_faces, _num_patches, _num_leaves
    , _num_planes, _num_nodes;
  bsp_face *_faces;
  bsp_mesh_face *_mesh_faces;
  bsp_patch *_patches;
  bsp_leaf *_leaves;
  bsp_leafface *_leaffaces;
  bsp_plane *_planes;
  bsp_node *_nodes;
  bsp_visdata *_visdata;
public:
  bsp(const char *filename);
  ~bsp();
};

