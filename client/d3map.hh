#pragma once

#include <string>
#include <vector>
#include <glm/glm.hpp>
#include "ogl.hh"
#include "frustum.hh"

class d3_lexer {
public:
  d3_lexer(std::string filename);
};

// TODO space waste?
struct d3_vertex {
  glm::vec3 pos;
  glm::vec2 texcoord;
  glm::vec3 normal;
};

struct d3_plane {
  glm::vec3 normal;
  float dist;

  void normalize();
  bool in_front(const glm::vec3 &point);
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

class d3_model {
  std::vector<d3_surface*> _surfaces;
public:
  std::vector<d3_portal*> portals;
  std::string name;
  int index;

  d3_model(const std::string &n_name, int n_index);
  void read_from_file(std::ifstream &file);
  void draw(const glm::vec3 &position) const;
};

class d3map;

class d3_portal {
  int _model_pos, _model_neg;
  std::vector<glm::vec3> _points;
  bool _visible;
public:
  void read_from_file(std::ifstream &file, d3map *map);
  void render_from_model(const glm::vec3 &position, int idx);
};

struct d3_node {
  d3_plane plane;
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
  void bind();
  void add_portal_to_model(d3_portal *p, int idx);
  int get_model_idx_by_name(const std::string &name); // TODO std::map
  d3map(const std::string &filename);
  void draw(const glm::vec3 &position, const glm::mat4 &mvp, const frustum &f);
};

