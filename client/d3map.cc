#include "d3map.hh"
#include "utils.hh"
#include "shaders.hh"
#include <fstream>
#include <glm/gtc/type_ptr.hpp>

bool d3_plane::in_front(const glm::vec3 &point) {
  return glm::dot(point, normal) + dist > 0;
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
    , positive_model = proc_get_next_int(file)
    , negative_model = proc_get_next_int(file);
  for (int i = 0; i < num_points; ++i) {
    glm::vec3 tmp;
    tmp.x = proc_get_next_float(file);
    tmp.y = proc_get_next_float(file);
    tmp.z = proc_get_next_float(file);
    _points.push_back(std::move(tmp));
  }

  std::string name_positive = "_area" + std::to_string(positive_model)
    , name_negative = "_area" + std::to_string(negative_model);
  _model_pos = map->get_model_idx_by_name(name_positive);
  _model_neg = map->get_model_idx_by_name(name_negative);
  if (_model_pos >= 0)
    map->add_portal_to_model(this, _model_pos);
  if (_model_neg >= 0)
    map->add_portal_to_model(this, _model_neg);
}

void d3_portal::render_from_model(const glm::vec3 &position, int idx) {
#if 0
  if (!(_visible = check_visibility(camera)))
    return; // portal is outside frustrum
  transform_points();

  if(!_visible)
    return;
  else if (_visible < 0) {
    // intersection of portal and front plane of frustum
    // set min and max to renderport
    m_transformed_min = glm::ivec2(0, 0);
    m_transformed_max = glm::ivec2(base->get_window_width(), base->get_window_height());
  }

  // clip min and max
  if (min.x < m_transformed_min.x) min.x = m_transformed_min.x;
  if (max.x > m_transformed_max.x) max.x = m_transformed_max.x;

  if (min.y < m_transformed_min.y) min.y = m_transformed_min.y;
  if (max.y > m_transformed_max.y) max.y = m_transformed_max.y;

  // render model if visible
  if ((max.x > min.x) && (max.y > min.y)) {
    if (index == m_model_pos)
      m_scene->get_model(m_model_neg)->render(camera, min, max);
    else
      m_scene->get_model(m_model_pos)->render(camera, min, max);
  }
#endif
}

d3_surface::d3_surface(const std::vector<d3_vertex> &vertices
    , const std::vector<unsigned int> &elements) {
  _vao.bind();

  array_buffer vbo;
  vbo.bind();
  glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertices[0])
      , vertices.data(), GL_STATIC_DRAW);

  element_array_buffer ebo;
  ebo.bind();
  ebo.upload(elements);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE
      , sizeof(d3_vertex), 0);
  glEnableVertexAttribArray(0);

  _vao.unbind();

  _num_elements = elements.size();
}

void d3_surface::draw() const {
  _vao.bind();
  glDrawElements(GL_TRIANGLES, _num_elements, GL_UNSIGNED_INT, 0);
}

d3_model::d3_model(const std::string &n_name, int n_index)
  : name(n_name)
  , index(n_index) {
}

void d3_model::read_from_file(std::ifstream &file) {
  const int num_surfaces = proc_get_next_int(file);
  for (int i = 0; i < num_surfaces; ++i) {
    std::vector<d3_vertex> vertices;
    std::vector<unsigned int> elements;

    std::string texture_filename = proc_get_next_string(file) + ".tga";

    const int num_verts = proc_get_next_int(file)
      , num_ind = proc_get_next_int(file);
    for (int j = 0; j < num_verts; ++j) {
      d3_vertex v;
      v.pos.z = proc_get_next_float(file);
      v.pos.x = proc_get_next_float(file);
      v.pos.y = proc_get_next_float(file);
      v.texcoord.x = proc_get_next_float(file);
      v.texcoord.y = proc_get_next_float(file);
      v.normal.x = proc_get_next_float(file);
      v.normal.y = proc_get_next_float(file);
      v.normal.z = proc_get_next_float(file);
      vertices.push_back(std::move(v));
    }
    for (int j = 0; j < num_ind; ++j)
      elements.push_back(proc_get_next_int(file));

    d3_surface *n_surf = new d3_surface(vertices, elements);
    _surfaces.push_back(n_surf);
  }
}

void d3_model::draw(const glm::vec3 &position) const {
  for (const d3_surface *s : _surfaces)
    s->draw();
  // for (d3_portal *p : portals)
  //   p->render_from_model(position, index);
}

void d3map::_load_proc(const std::string &filename) {
  std::ifstream ifs(filename.c_str(), std::ios::binary);
  if (!ifs)
    die("failed to open \"%s\"", filename.c_str());

  std::string t;
  while (ifs >> t) {
    if (t == "model") {
      std::string name = proc_get_next_string(ifs);
      d3_model m(name, _models.size());
      m.read_from_file(ifs);
      _models.push_back(std::move(m));
    } else if (t == "interAreaPortals") {
      const int num_models = proc_get_next_int(ifs),
            num_portals = proc_get_next_int(ifs);
      for (int i = 0; i < num_portals; ++i) {
        d3_portal portal;
        portal.read_from_file(ifs, this);
        _portals.push_back(std::move(portal));
      }
    } else if (t == "nodes") {
      const int num_nodes = proc_get_next_int(ifs);
      for (int i = 0; i < num_nodes; ++i) {
        d3_node node;
        node.plane.normal.z = proc_get_next_float(ifs);
        node.plane.normal.x = proc_get_next_float(ifs);
        node.plane.normal.y = proc_get_next_float(ifs);
        node.plane.dist = proc_get_next_float(ifs);
        node.positive_child = proc_get_next_int(ifs);
        node.negative_child = proc_get_next_int(ifs);
        if (node.positive_child < 0) {
          std::string name = "_area" + std::to_string(-1 - node.positive_child);
          node.positive_child = -1 - get_model_idx_by_name(name);
        }
        if (node.negative_child < 0) {
          std::string name = "_area" + std::to_string(-1 - node.negative_child);
          node.negative_child = -1 - get_model_idx_by_name(name);
        }
        _nodes.push_back(std::move(node));
      }
    } else {
      // token skipped
    }
  }
  ifs.close();
}

int d3map::_get_model_idx_by_pos(const glm::vec3 &position) {
  if (!_nodes.size())
    return 0;
  d3_node *node = &_nodes[0];
  while (true)
    if (node->plane.in_front(position)) { // in front
      if (node->positive_child > 0)
        node = &_nodes[node->positive_child];
      else
        return node->positive_child;
    } else { // behind
      if (node->negative_child > 0)
        node = &_nodes[node->negative_child];
      else
        return node->negative_child;
    }
}

void d3map::add_portal_to_model(d3_portal *p, int idx) {
  _models[idx].portals.push_back(p);
}

int d3map::get_model_idx_by_name(const std::string &name) {
  for (size_t i = 0; i < _models.size(); ++i)
    if (_models[i].name == name)
      return i;
  return -1;
}

d3map::d3map(const std::string &filename)
  : _sp(shaders::simple_vert, shaders::simple_frag) {

  _sp.use_this_prog();
  _vertex_pos_attr = _sp.bind_attrib("vertex_pos");
  _mvp_mat_unif = _sp.bind_uniform("mvp");

  this->_load_proc(filename + ".proc");

  _sp.dont_use_this_prog();
}

void d3map::draw(const glm::vec3 &position, const glm::mat4 &mvp
    , const frustum &f) {
  // _set_visible_faces(position, f);

  _sp.use_this_prog();

  glUniformMatrix4fv(_mvp_mat_unif, 1, GL_FALSE, glm::value_ptr(mvp));

  int start_model = _get_model_idx_by_pos(position);
  if (start_model >= 0)
    // position is in the void
    for (const d3_model &m : _models)
      m.draw(position);
  else
    _models[-1 - start_model].draw(position);

  _sp.dont_use_this_prog();
}

