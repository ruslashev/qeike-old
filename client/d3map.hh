#pragma once

#include <string>
#include <vector>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/glm.hpp>

class d3_lexer {
public:
  d3_lexer(std::string filename);
};

struct d3_vertex {
  glm::vec3 pos;
  glm::vec2 texcoord;
  glm::vec3 normal;
};

struct d3_plane {
  glm::vec3 normal;
  float dist;

  void normalize();
};

class d3map;

class d3_portal {
  int _area_pos, _area_neg;
  std::vector<glm::vec3> _points;
public:
  void read_from_file(std::ifstream &file, d3map *map);
};

class d3_area {
public:
  std::vector<d3_portal*> portals;
  std::string name;
  int index;

  d3_area(const std::string &n_name, int n_index);
  void read_from_file(std::ifstream &file);
};

struct d3_node {
  d3_plane plane;
  int pos_child, neg_child;
};

class d3map {
  std::vector<d3_area> _areas;
  std::vector<d3_portal> _portals;
  std::vector<d3_node> _nodes;

  void _load_proc(const std::string &filename);
public:
  void add_portal_to_area(d3_portal *p, int idx);
  int get_area_index_by_name(const std::string &name); // TODO std::map
  d3map(const std::string &filename);
};

