#include "map.hh"
#include "utils.hh"
#include "shaders.hh"
#include <fstream>
#include <cstring>
#include <glm/gtc/type_ptr.hpp>

struct header {
  char magic[9];
  uint8_t version;
  uint32_t tiles;
  uint16_t width, height, depth;
};

struct vertex {
  glm::vec3 position;
  glm::vec3 normal;
};

void map::_generate_mesh() {
  const std::vector<glm::vec3> base_verts = {
    glm::vec3(0, 0, 1), glm::vec3(1, 0, 1), glm::vec3(1, 1, 1), glm::vec3(0, 1, 1),
    glm::vec3(0, 0, 0), glm::vec3(1, 0, 0), glm::vec3(1, 1, 0), glm::vec3(0, 1, 0)
  };
  const std::vector<GLushort> base_elements = {
    0, 1, 2, 2, 3, 0, 1, 5, 6, 6, 2, 1, 7, 6, 5, 5, 4, 7,
    4, 0, 3, 3, 7, 4, 4, 5, 1, 1, 0, 4, 3, 2, 6, 6, 7, 3,
  };
  int element_offset = 0;
  std::vector<vertex> vertices;
  std::vector<GLushort> elements;

  for (int z = 0; z < _depth; ++z)
    for (int y = 0; y < _height; ++y)
      for (int x = 0; x < _width; ++x)
        if (_data[z][y][x]) {
          for (const glm::vec3 &v : base_verts) {
            vertex transf_vert;
            transf_vert.normal = glm::vec3(0, 0, 0);
            transf_vert.position = v + glm::vec3(x, y, z);
            vertices.push_back(transf_vert);
          }
          for (GLushort e : base_elements)
            elements.push_back(e + element_offset);
          element_offset += 8;
        }

  _vao.bind();

  _sp.use_this_prog();
  _vertex_pos_attr = _sp.bind_attrib("vertex_pos");
  _mvp_mat_unif = _sp.bind_uniform("mvp");

  _vbo.bind();
  _vbo.upload(vertices.size() * sizeof(vertices[0]), &vertices[0]);
  _ebo.bind();
  _ebo.upload(elements);
  _num_elements = elements.size();

  glEnableVertexAttribArray(_vertex_pos_attr);
  glVertexAttribPointer(_vertex_pos_attr, 3, GL_FLOAT, GL_FALSE
      , sizeof(vertex), (void*)offsetof(vertex, position));

  _sp.dont_use_this_prog();
  _vao.unbind();
  _vbo.unbind();
  _ebo.unbind();
}

map::map(std::string filename)
  : _width(0)
  , _height(0)
  , _depth(0)
  , _sp(shaders::simple_vert, shaders::simple_frag) {
  std::ifstream ifs(filename, std::ios::binary);
  if (!ifs)
    die("failed to open map \"%s\"", filename.c_str());

  header h;
  ifs.read((char*)&h, sizeof(header));
  if (strcmp(h.magic, "QEIKEMAP") != 0 || h.version != 1) {
    ifs.close();
    die("invalid magic at \"%s\": %s %d", filename.c_str(), h.magic, h.version);
  }

  _width = h.width;
  _height = h.height;
  _depth = h.depth;

  _data.resize(_depth, std::vector<std::vector<int>>(_height
        , std::vector<int>(_width, 0)));

  for (unsigned int i = 0; i < h.tiles; ++i) {
    uint8_t x, y, z;
    ifs.read((char*)&x, sizeof(x));
    ifs.read((char*)&y, sizeof(y));
    ifs.read((char*)&z, sizeof(z));
    _data[z][y][x] = 1;
  }

  ifs.close();

  _generate_mesh();
}

void map::draw(const glm::mat4 &mvp) {
  _vao.bind();
  _sp.use_this_prog();

  glUniformMatrix4fv(_mvp_mat_unif, 1, GL_FALSE, glm::value_ptr(mvp));

  glDrawElements(GL_TRIANGLES, _num_elements, GL_UNSIGNED_SHORT, 0);

  _vao.unbind();
  _sp.dont_use_this_prog();
}

