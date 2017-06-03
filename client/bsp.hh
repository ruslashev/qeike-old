#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <map>
#include "ogl.hh"

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
  glm::ivec3 mins;
  glm::ivec3 maxs;
};

struct bsp_texture {
  char name[64];
  int flags;
  int contents;
};

struct bsp_leaf {
  int cluster;
  int area;
  glm::ivec3 mins;
  glm::ivec3 maxs;
  int leafface;
  int n_leaffaces;
  int leafbrush;
  int n_leafbrushes;
};

struct bsp_leafface {
  int face;
};

struct bsp_vertex {
  glm::vec3 position;
  glm::vec2 decal;
  glm::vec2 lightmap;
  glm::vec3 normal;
  unsigned char color[4];
#if 0
  bsp_vertex operator+(const bsp_vertex &v) const;
  bsp_vertex operator*(float factor) const;
#endif
};

struct bsp_meshvert {
  int offset;
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
  glm::ivec2 lm_start;
  glm::ivec2 lm_size;
  glm::vec3 lm_origin;
  glm::vec3 lm_s;
  glm::vec3 lm_t;
  glm::vec3 normal;
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
  std::vector<bsp_texture> _textures;
  std::vector<bsp_plane> _planes;
  std::vector<bsp_node> _nodes;
  std::vector<bsp_leaf> _leaves;
  std::vector<bsp_leafface> _leaffaces;
  std::vector<bsp_vertex> _vertices;
  std::vector<bsp_meshvert> _meshverts;
  std::vector<bsp_face> _faces;
  std::vector<int> _visible_faces;
  std::vector<bsp_lightmap> _lightmaps;
  std::vector<GLuint> _lightmap_texture_ids;
  bsp_visdata _visdata;
  GLint _vertex_pos_attr, _texture_coord_attr, _lightmap_coord_attr
    , _mvp_mat_unif;

  void _load_file(const char *filename, float world_scale);
  int _find_leaf(glm::vec3 position);
  int _cluster_visible(int vis_cluster, int test_cluster);
  void _set_visible_faces(glm::vec3 camera_pos);
public:
  std::vector<GLuint> texture_ids;
  shader_program sp;
  array_buffer vbo;

  bsp(const char *filename, float world_scale);
  ~bsp();
  void render(glm::vec3 position, const glm::mat4 &mvp);
};

