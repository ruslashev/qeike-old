#include "d3map.hh"
#include "utils.hh"
#include <fstream>
#include <sstream>

d3_lexer::d3_lexer(std::string filename) {
  // TODO
}

void d3_plane::normalize() {
  // TODO make faster
  float inv_length = 1.0f / glm::length(normal);
  normal *= inv_length;
  dist *= inv_length;
}

std::string proc_get_next_value(std::ifstream &file) {
  std::string s;
  while (file >> s) {
    if (s == "/*") {
      // skip comments
      while (s != "*/")
        file >> s;
    } else if (s == "{" || s == "}") {
      // skip braces
    } else if (s == "(" || s == ")") {
      // skip parentheses
    } else
      return s;
  }
  die("unexpected eof");
}

std::string proc_get_next_string(std::ifstream &file) {
  std::string s = proc_get_next_value(file);
  return s.substr(1, s.size() - 2);
}

#define proc_get_next_float(x) atof(proc_get_next_value(x).c_str())
#define proc_get_next_int(x) atoi(proc_get_next_value(x).c_str())

void d3_portal::read_from_file(std::ifstream &file, d3map *map) {
  const int num_points = proc_get_next_int(file)
    , pos_area = proc_get_next_int(file), neg_area = proc_get_next_int(file);
  for (int i = 0; i < num_points; ++i) {
    glm::vec3 tmp;
    tmp.x = proc_get_next_float(file);
    tmp.y = proc_get_next_float(file);
    tmp.z = proc_get_next_float(file);
    _points.push_back(tmp); // TODO std::move
  }

  std::stringstream name_pos, name_neg;
  name_pos << "_area" << pos_area;
  _area_pos = map->get_area_index_by_name(name_pos.str());
  name_neg << "_area" << neg_area;
  _area_neg = map->get_area_index_by_name(name_neg.str());

  if (_area_pos >= 0)
    map->add_portal_to_area(this, _area_pos);
  if (_area_neg >= 0)
    map->add_portal_to_area(this, _area_neg);
}

d3_area::d3_area(const std::string &n_name, int n_index)
  : name(n_name)
  , index(n_index) {
}

void d3_area::read_from_file(std::ifstream &file) {
  const int num_surfaces = proc_get_next_int(file);
  for (int i = 0; i < num_surfaces; ++i) {
    std::vector<d3_vertex> vertices;
    std::vector<unsigned int> indices;

    std::string texture_filename = proc_get_next_string(file) + ".tga";

    const int num_verts = proc_get_next_int(file)
      , num_ind = proc_get_next_int(file);

    for (int j = 0; j < num_verts; ++j) {
      d3_vertex v;
      v.pos.x = proc_get_next_float(file);
      v.pos.y = proc_get_next_float(file);
      v.pos.z = proc_get_next_float(file);
      v.texcoord.x = proc_get_next_float(file);
      v.texcoord.y = proc_get_next_float(file);
      v.normal.x = proc_get_next_float(file);
      v.normal.y = proc_get_next_float(file);
      v.normal.z = proc_get_next_float(file);
      vertices.push_back(v); // TODO std::move
    }
    for (int j = 0; j < num_ind; ++j)
      indices.push_back(proc_get_next_int(file));
  }
}

void d3map::_load_proc(const std::string &filename) {
  std::ifstream ifs(filename.c_str(), std::ios::binary);
  if (!ifs)
    die("failed to open \"%s\"", filename.c_str());

  std::string t;
  while (ifs >> t) {
    if (t == "model") {
      std::string name = proc_get_next_string(ifs);
      d3_area a(name, _areas.size());
      a.read_from_file(ifs);
      _areas.push_back(a); // TODO std::move
    } else if (t == "interAreaPortals") {
      const int num_areas = proc_get_next_int(ifs),
            num_portals = proc_get_next_int(ifs);
      for (int i = 0; i < num_portals; ++i) {
        d3_portal portal;
        portal.read_from_file(ifs, this);
        _portals.push_back(portal); // TODO std::move
      }
    } else if (t == "nodes") {
      const int num_nodes = proc_get_next_int(ifs);
      for (int i = 0; i < num_nodes; ++i) {
        d3_node node;
        node.plane.normal.x = proc_get_next_float(ifs);
        node.plane.normal.y = proc_get_next_float(ifs);
        node.plane.normal.z = proc_get_next_float(ifs);
        node.plane.dist = proc_get_next_float(ifs);
        node.pos_child = proc_get_next_int(ifs);
        node.neg_child = proc_get_next_int(ifs);
        if (node.pos_child < 0) {
          std::stringstream name;
          name << "_area" << (-1 - node.pos_child);
          node.pos_child = -1 - get_area_index_by_name(name.str());
        }
        if (node.neg_child < 0) {
          std::stringstream name;
          name << "_area" << (-1 - node.neg_child);
          node.neg_child = -1 - get_area_index_by_name(name.str());
        }
        _nodes.push_back(node); // TODO std::move
      }
    } else {
      // token skipped
    }
  }
  ifs.close();
}

void d3map::add_portal_to_area(d3_portal *p, int idx) {
  _areas[idx].portals.push_back(p);
}

int d3map::get_area_index_by_name(const std::string &name) {
  for (size_t i = 0; i < _areas.size(); ++i)
    if (_areas[i].name == name)
      return i;
  return -1;
}

d3map::d3map(const std::string &filename) {
  std::string proc_filename = filename + ".proc";
  this->_load_proc(proc_filename);
}

