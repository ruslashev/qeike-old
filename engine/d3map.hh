#pragma once

#include <string>
#include <vector>
#include <glm/glm.hpp>
#include "ogl.hh"
#include "frustum.hh"

namespace qke {

struct d3_vertex {
  glm::vec3 pos;
  glm::vec2 texcoord;
  glm::vec3 normal;
};

struct d3_plane { // TODO use own
  glm::vec3 normal;
  float dist;

  bool in_front(const glm::vec3 &point);
};

class d3_surface {
  vertex_array_object _vao;
  int _num_elements;
public:
  d3_surface(const std::vector<d3_vertex> &vertices
      , const std::vector<unsigned int> &elements);
  void draw() const;
};

class d3_portal;
class d3map;

class d3_model {
  std::vector<d3_surface*> _surfaces;
  unsigned long long int _rendered_frame_idx;
public:
  std::vector<d3_portal*> portals;
  std::string name;
  int index;

  d3_model(const std::string &n_name, int n_index);
  void read_from_file(std::ifstream &file);
  void draw(const glm::mat4 &mvp, glm::ivec2 min, glm::ivec2 max, d3map *map);
};

class d3map;

class d3_portal {
  int _model_pos, _model_neg;
  std::vector<glm::vec3> _points;
  std::vector<glm::ivec2> _transformed_points;
  glm::ivec2 _transformed_min, _transformed_max;
  bool _visible;
  unsigned long long int _rendered_frame_idx;
public:
  d3_portal();
  void read_from_file(std::ifstream &file, d3map *map);
  void draw_from_model(int idx, const glm::mat4 &mvp, glm::ivec2 min
      , glm::ivec2 max, d3map *map);
  void transform_points(const glm::mat4 &mvp);
  bool project(const glm::vec4 &vec, int &x, int &y, const glm::mat4 &mvp);
};

struct d3_node {
  d3_plane plane;
  int positive_child, negative_child;
};

class d3map {
  std::vector<d3_model> _models;
  std::vector<d3_portal*> _portals; // TODO rework this
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
  void draw(const glm::vec3 &position, const glm::mat4 &mvp, const frustum &f);
};

};

