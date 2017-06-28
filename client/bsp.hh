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

  struct bsp_patch_plane {
    glm::vec4 plane;
    // signx + (signy << 1) + (signz << 2), used as lookup during collision
    int sign_bits;
  };

  struct bsp_facet {
    int surface_plane;
    int n_borders; // 3 or four + 6 axial bevels + 4 or 3 * 4 edge bevels
    int border_planes[4 + 6 + 16];
    int border_inward[4 + 6 + 16];
    bool border_no_adjust[4 + 6 + 16];
  };

  struct bsp_patchcollide {
    bool valid;
    glm::vec3 bounds[2];
    int n_planes; // surface planes plus edge planes
    bsp_patch_plane *planes;
    int n_facets;
    bsp_facet *facets;
    // std::vector<bsp_patch_plane> patch_planes;
    // std::vector<bsp_facet> facets;
    void clear_bounds();
    void add_point_to_bounds(const glm::vec3 &v);
  };

#define MAX_GRID_SIZE 129
  struct bsp_grid {
    int width;
    int height;
    bool wrap_width;
    bool wrap_height;
    glm::vec3 points[MAX_GRID_SIZE][MAX_GRID_SIZE]; // [width][height]

    void set_wrap_width();
    void subdivide_columns();
    void remove_degenerate_columns();
    void transpose();
  };

  struct bsp_winding {
    std::vector<glm::vec3> p;
    int numpoints; // can be less than p.size()

    bsp_winding(int size);
    void bounds(glm::vec3 & mins, glm::vec3 & maxs);
  };

  enum edge_name {
    EN_TOP,
    EN_RIGHT,
    EN_BOTTOM,
    EN_LEFT
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
  std::vector<bsp_patchcollide> _patchcollides;
  bsp_visdata _visdata;

#define MAX_PATCH_PLANES 2048
#define MAX_FACETS 1024
  int numPlanes;
  bsp_patch_plane planes[MAX_PATCH_PLANES];
  int numFacets;
  bsp_facet facets[MAX_PATCH_PLANES]; //maybe MAX_FACETS ??

  GLint _vertex_pos_attr, _lightmap_coord_attr, _mvp_mat_unif;
  shader_program _sp;
  array_buffer _vbo;
  element_array_buffer _ebo;
  vertex_array_object _vao;

  void _load_file(const char *filename, float world_scale
      , int tesselation_level);
  void _create_patches(int tesselation_level);
  void _tesselate(int tesselation_level, int control_offset, int control_width
      , int voff, int eoff);

  void _create_patch_collides();
  void _create_patch_collide_from_grid(bsp_grid *grid, bsp::bsp_patchcollide *pf);
  int find_plane( glm::vec3 p1, glm::vec3 p2, glm::vec3 p3 );
  int _edge_plane_for_num(bsp_grid *grid
      , int gridPlanes[MAX_GRID_SIZE][MAX_GRID_SIZE][2], int i, int j, int k );
  void set_border_inward( bsp_facet *facet, bsp_grid *grid, int gridPlanes[MAX_GRID_SIZE][MAX_GRID_SIZE][2], int i, int j, int which );
  int point_on_plane_side( glm::vec3 p, int planeNum );
  bool validate_facet( bsp_facet *facet );
  bsp_winding *base_winding_for_plane (glm::vec3 normal, float dist);
  void chop_winding_in_place (bsp_winding **inout, glm::vec3 normal, float dist, float epsilon);
  void add_facet_bevels( bsp_facet *facet );
  int plane_equal(bsp_patch_plane *p, glm::vec4 plane, int *flipped);
  int find_plane2(glm::vec4 plane, int *flipped);

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
  void _trace_patch(trace_result *tr, const trace_description &td, bsp_patchcollide *p, glm::vec3 input_start, glm::vec3 input_end);
  bool _check_facet_plane(glm::vec3 plane, glm::vec3 start, glm::vec3 end, float *enterFrac, float *leaveFrac, bool *hit);
public:
  bsp(const char *filename, float world_scale, int tesselation_level);
  ~bsp();
  void draw(const glm::vec3 &position, const glm::mat4 &mvp, const frustum &f);
  void trace_sphere(trace_result *tr, const glm::vec3 &start
      , const glm::vec3 &end, float radius);
};

