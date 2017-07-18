#pragma once

#include <string>
#include <vector>
#include <glm/glm.hpp>
#include "ogl.hh"
#include "frustum.hh"

struct d3_vertex {
  glm::vec3 pos;
  glm::vec2 texcoord;
  glm::vec3 normal;
};

class d3_surface {
  vertex_array_object _vao;
  int _num_elements;
public:
  d3_surface(const std::vector<d3_vertex> &vertices
      , const std::vector<unsigned int> &elements);
  ~d3_surface();
  void draw() const;
};

class d3_portal;
class d3map;

class d3_model {
  std::vector<d3_surface*> _surfaces;
public:
  std::vector<d3_portal*> portals;
  std::string name;
  int index; // TODO private members?

  d3_model(const std::string &n_name, int n_index);
  void read_from_file(std::ifstream &file);
  void draw(const glm::ivec2 &min, const glm::ivec2 &max, const frustum &f
      , const glm::mat4 &projection, d3map *map);
};

class d3_portal {
  enum frustum_visibility {
    FRUSTUM_VISIBILITY_OUTSIDE,
    FRUSTUM_VISIBILITY_INSIDE,
    FRUSTUM_VISIBILITY_INTERSECTS
  };

  int _model_pos, _model_neg; // TODO needed as member?
  std::vector<glm::vec3> _points;
  qkmath::plane _plane;
  glm::ivec2 _transformed_min, _transformed_max;

  void _create_plane_from_points();
  frustum_visibility _check_visibility(const frustum &f) const;
  void _transform_points(const glm::mat4 &projection);
public:
  void read_from_file(std::ifstream &file, d3map *map);
  void draw_from_model(int idx, const glm::ivec2 &min
      , const glm::ivec2 &max, const frustum &f
      , const glm::mat4 &projection, d3map *map);
};

struct d3_node {
  qkmath::plane plane;
  int positive_child, negative_child;
};

class d3map {
  std::vector<d3_model> _models;
  std::vector<d3_portal> _portals;
  std::vector<d3_node> _nodes;

  GLint _vertex_pos_attr, _mvp_mat_unif;
  shader_program _sp;

  void _load_proc(const std::string &filename);
  int _get_model_idx_by_pos(const glm::vec3 &position);
public:
  d3map(const std::string &filename);
  void add_portal_to_model(d3_portal *p, int idx);
  int get_model_idx_by_name(const std::string &name); // TODO std::map
  d3_model* get_model_by_idx(int idx);
  void draw(const glm::mat4 &projection, const glm::mat4 &view
      , const frustum &f);
};

