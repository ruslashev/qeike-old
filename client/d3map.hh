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

class d3map;

class d3_portal {
  int _area_pos, _area_neg;
  std::vector<glm::vec3> _points;
  bool _visible;
public:
  void read_from_file(std::ifstream &file, d3map *map);
  void render_from_area(const glm::vec3 &position, int idx);
};

class d3_area {
public:
  std::vector<d3_portal*> portals;
  std::string name;
  int index;

  d3_area(const std::string &n_name, int n_index);
  void read_from_file(std::ifstream &file, std::vector<d3_vertex> &vertices
      , std::vector<unsigned int> &elements);
  void draw(const glm::vec3 &position) const;
};

struct d3_node {
  d3_plane plane;
  int pos_child, neg_child;
};

class d3map {
  std::vector<d3_area> _areas;
  std::vector<d3_portal> _portals;
  std::vector<d3_node> _nodes;

  GLint _vertex_pos_attr, _mvp_mat_unif;
  shader_program _sp;
  vertex_array_object _vao;
  array_buffer _vbo;
  element_array_buffer _ebo;

  void _load_proc(const std::string &filename);
  int _get_area_idx_by_pos(const glm::vec3 &position);
public:
  void add_portal_to_area(d3_portal *p, int idx);
  int get_area_idx_by_name(const std::string &name); // TODO std::map
  d3map(const std::string &filename);
  void draw(const glm::vec3 &position, const glm::mat4 &mvp, const frustum &f);
};

