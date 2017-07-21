#pragma once

#include <glm/glm.hpp>
#include <vector>
#include "frustum.hh"
#include "ogl.hh"

struct trace_result {
  glm::vec3 clip_plane_normal;
  float fraction; // 1 - didn't hit anything
  glm::vec3 end;
  bool start_solid; // true - initial point was in a solid area
  bool all_solid; // true - plane is not valid TODO: rename to not_valid
};

class bsp {
  struct bsp_shader {
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
    glm::vec3 mins;
    glm::vec3 maxs;
  };

  struct bsp_leaf {
    int cluster;
    int area;
    glm::vec3 mins;
    glm::vec3 maxs;
    int leafface;
    int n_leaffaces;
    int leafbrush;
    int n_leafbrushes;
  };

  struct bsp_leafface {
    int face;
  };

  typedef int bsp_leafbrush;

  struct bsp_brush {
    int brushside;
    int n_brushsides;
    int shader;
  };

  struct bsp_brushside {
    int plane;
    int shader;
  };

  struct bsp_vertex {
    glm::vec3 position;
    glm::vec2 lightmap;
    bsp_vertex operator+(const bsp_vertex &v) const;
    bsp_vertex operator*(float factor) const;
  };

  struct bsp_meshvert {
    int offset;
  };

  struct bsp_face {
    int shader;
    int fog;
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

  enum class trace_type {
    ray,
    sphere,
  };

  struct trace_description {
    trace_type type;
    float radius;
  };

  std::vector<bsp_shader> _shaders;
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

  GLint _vertex_pos_attr, _lightmap_coord_attr, _mvp_mat_unif;
  shader_program _sp;
  array_buffer _vbo;
  element_array_buffer _ebo;
  vertex_array_object _vao;

  void _load_file(const char *filename, int tesselation_level);
  void _create_patches(int tesselation_level);
  void _tesselate(int tesselation_level, int control_offset, int control_width
      , int voff, int eoff);

  int _find_leaf(glm::vec3 position);
  int _cluster_visible(int vis_cluster, int test_cluster);
  void _set_visible_faces(const glm::vec3 &camera_pos, const frustum &f);

  void _trace(trace_result *tr, const trace_description &td
      , const glm::vec3 &start, const glm::vec3 &end);
  void _trace_node(trace_result *tr, const trace_description &td, int node_index
      , float start_fraction, float end_fraction, const glm::vec3 &start
      , const glm::vec3 &end);
  void _trace_brush(trace_result *tr, const trace_description &td
      , const bsp_brush *b, const glm::vec3 &input_start
      , const glm::vec3 &input_end);
  void _trace_leaf(trace_result *tr, const trace_description &td
      , const bsp_leaf *l, const glm::vec3 &start, const glm::vec3 &end);
  bool _check_facet_plane(glm::vec3 plane, glm::vec3 start, glm::vec3 end
      , float *enter_frac, float *leave_frac, bool *hit);
public:
  bsp(const char *filename, float world_scale, int tesselation_level);
  ~bsp();
  void draw(const glm::vec3 &position, const glm::mat4 &mvp, const frustum &f);
  void trace_sphere(trace_result *tr, const glm::vec3 &start
      , const glm::vec3 &end, float radius);
};

