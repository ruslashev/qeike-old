#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <map>
#include "ogl.hh"

struct bsp_texture {
  char name[64];
  int flags;
  int contents;
};

struct bsp_plane {
  glm::vec3 normal;
  float dist;
  bsp_plane();
};

struct bsp_node {
  int plane;
  int front;
  int back;
  glm::ivec3 mins;
  glm::ivec3 maxs;
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

struct bsp_leafbrush {
  int brush;
};

struct bsp_brush {
  int brushside;
  int n_brushsides;
  int texture;
};

struct bsp_brushside {
  int plane;
  int texture;
};

struct bsp_vertex {
  glm::vec3 position;
  glm::vec2 decal;
  glm::vec2 lightmap;
  glm::vec3 normal;
  unsigned char color[4];
  bsp_vertex operator+(const bsp_vertex &v) const;
  bsp_vertex operator*(float factor) const;
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

class bsp_biquadratic_patch {
  int _tesselation_level;
  std::vector<GLuint> _indices;
  // store as pointer arrays for easier access by GL functions
  int *_triangles_per_row;
  unsigned int **_row_index_pointers;
public:
  bsp_vertex control_points[9];
  std::vector<bsp_vertex> vertices;

  bsp_biquadratic_patch();
  ~bsp_biquadratic_patch();
  void tesselate(int level);
  void render() const;
};

struct bsp_patch {
  int texture_idx;
  int lightmap_idx;
  int width;
  int height;
  std::vector<bsp_biquadratic_patch> quadratic_patches;
};

enum class trace_type {
  ray,
  sphere,
};

struct trace_description {
  trace_type type;
  float radius;
};

struct trace_result {
  glm::vec3 clip_plane_normal;
  float fraction;
  glm::vec3 end;
  bool start_solid;
  bool all_solid;
};

class bsp {
  std::vector<bsp_texture> _textures;
  std::vector<bsp_plane> _planes;
  std::vector<bsp_node> _nodes;
  std::vector<bsp_leaf> _leaves;
  std::vector<bsp_leafface> _leaffaces;
  std::vector<bsp_leafbrush> _leafbrushes;
  std::vector<bsp_brush> _brushes;
  std::vector<bsp_brushside> _brushsides;
  std::vector<bsp_vertex> _vertices;
  std::vector<bsp_meshvert> _meshverts;
  std::vector<bsp_face> _faces;
  std::vector<int> _visible_faces;
  std::vector<bsp_lightmap> _lightmaps;
  std::vector<GLuint> _lightmap_texture_ids;
  bsp_visdata _visdata;
  std::vector<bsp_patch*> _patches;
  GLint _vertex_pos_attr, /* _texture_coord_attr, */ _lightmap_coord_attr
    , _mvp_mat_unif;
  vertex_array_object _vao;
  element_array_buffer _ebo;
  array_buffer _vbo;

  void _tesselate(int level, int controlOffset, int controlWidth, int vOffset, int iOffset);
  void _load_file(const char *filename, float world_scale
      , int tesselation_level);
  int _find_leaf(glm::vec3 position);
  int _cluster_visible(int vis_cluster, int test_cluster);
  void _set_visible_faces(glm::vec3 camera_pos);
  void _create_patch(const bsp_face &f, int tesselation_level);
  void _trace(trace_result *tr, const trace_description &td
      , const glm::vec3 &start, const glm::vec3 &end);
  void _check_node(trace_result *tr, const trace_description &td, int node_index
      , float start_fraction, float end_fraction, const glm::vec3 &start
      , const glm::vec3 &end);
  void _check_brush(trace_result *tr, const trace_description &td, bsp_brush *b
      , const glm::vec3 &input_start, const glm::vec3 &input_end);
public:
  std::vector<GLuint> texture_ids;
  shader_program sp;

  bsp(const char *filename, float world_scale, int tesselation_level);
  ~bsp();
  void render(glm::vec3 position, const glm::mat4 &mvp);
  void trace_sphere(trace_result *tr, const glm::vec3 &start
    , const glm::vec3 &end, float radius);
};

