#pragma once

#include "ogl.hh"
#include <string>
#include <vector>

typedef std::vector<std::vector<std::vector<int>>> map_t;

class map {
  int _width, _height, _depth;
  map_t _data;
  vertex_array_object _vao;
  array_buffer _vbo;
  element_array_buffer _ebo;
  shader_program _sp;
  GLint _vertex_pos_attr, _mvp_mat_unif;
  int _num_elements;

  void _generate_mesh();
public:
  map(std::string filename);
  void draw(const glm::mat4 &mvp);
};

