#include "d3map.hh"
#include "utils.hh"
#include "shaders.hh"
#include "screen.hh"
#include <fstream>
#include <glm/gtc/type_ptr.hpp>

namespace qke {

bool d3_plane::in_front(const glm::vec3 &point) {
  return glm::dot(point, normal) + dist > 0;
}

std::string proc_get_next_value(std::ifstream &file) {
  std::string s;
  while (file >> s) {
    if (s == "/*") { // skip comments
      while (s != "*/")
        file >> s;
    } else if (s == "{" || s == "}") { // skip braces
    } else if (s == "(" || s == ")") { // skip parentheses
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

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(d3_vertex), 0);
  glEnableVertexAttribArray(0);

  _vao.unbind();

  _num_elements = elements.size();
}

void d3_surface::draw() const {
  _vao.bind();
  glDrawElements(GL_TRIANGLES, _num_elements, GL_UNSIGNED_INT, 0);
}

d3_model::d3_model(const std::string &n_name, int n_index)
  : _rendered_frame_idx(-1)
  , name(n_name)
  , index(n_index) {
}

void d3_model::read_from_file(std::ifstream &file) {
  const int num_surfaces = proc_get_next_int(file);
  for (int i = 0; i < num_surfaces; ++i) {
    std::vector<d3_vertex> vertices;
    std::vector<unsigned int> elements;

    std::string texture_filename = proc_get_next_string(file) + ".tga";

    const int num_verts = proc_get_next_int(file)
      , num_elem = proc_get_next_int(file);
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
    for (int j = 0; j < num_elem; ++j)
      elements.push_back(proc_get_next_int(file));

    // TODO std::move(vertices)?
    d3_surface *n_surf = new d3_surface(vertices, elements);
    _surfaces.push_back(n_surf);
  }
}

void d3_model::draw(const glm::mat4 &mvp, glm::ivec2 min, glm::ivec2 max
    , d3map *map) {
  if (_rendered_frame_idx == g_screen->get_frame_idx())
    return;
  _rendered_frame_idx = g_screen->get_frame_idx();
  // TODO try without
  // glScissor(min.x, min.y, max.x - min.x + 1, max.y - min.y + 1);
  for (const d3_surface *s : _surfaces)
    s->draw();
  for (size_t i = 0; i < portals.size(); ++i)
    portals[i]->draw_from_model(i, mvp, min, max, map);
}

d3_portal::d3_portal()
  : _rendered_frame_idx(-1) {
}

void d3_portal::read_from_file(std::ifstream &file, d3map *map) {
  const int num_points = proc_get_next_int(file)
    , positive_model = proc_get_next_int(file)
    , negative_model = proc_get_next_int(file);
  for (int i = 0; i < num_points; ++i) {
    glm::vec3 tmp_point;
    tmp_point.x = proc_get_next_float(file);
    tmp_point.y = proc_get_next_float(file);
    tmp_point.z = proc_get_next_float(file);
    _points.push_back(std::move(tmp_point));
    _transformed_points.push_back(glm::ivec2(0, 0));
  }

  std::string name_positive = "_area" + std::to_string(positive_model)
    , name_negative = "_area" + std::to_string(negative_model);
  _model_pos = map->get_model_idx_by_name(name_positive);
  _model_neg = map->get_model_idx_by_name(name_negative);
  if (_model_pos >= 0)
    map->add_portal_to_model(this, _model_pos);
  if (_model_neg >= 0)
    map->add_portal_to_model(this, _model_neg);

  _visible = true; // TODO add frustum culling and remove this line
}

void d3_portal::draw_from_model(int idx, const glm::mat4 &mvp, glm::ivec2 min
    , glm::ivec2 max, d3map *map) {
  // this check is crucial, otherwise infinite loops will occur
  // or check in dr_model
  if (_rendered_frame_idx != g_screen->get_frame_idx()) {
    _rendered_frame_idx = g_screen->get_frame_idx();
    // TODO
    // if (!(_visible = check_visibility(camera)))
    //   return; // portal is outside frustrum
    transform_points(mvp);
  }

  // TODO
  // if(!_visible)
  //   return;
  // else if (_visible < 0) {
    // intersection of portal and front plane of frustum
    // set min and max to renderport
    // _transformed_min = glm::ivec2(0, 0);
    // _transformed_max = glm::ivec2(g_screen->get_window_width(), g_screen->get_window_height());
  // }

  // TODO: clamp
  if (min.x < _transformed_min.x)
    min.x = _transformed_min.x;
  if (max.x > _transformed_max.x)
    max.x = _transformed_max.x;
  if (min.y < _transformed_min.y)
    min.y = _transformed_min.y;
  if (max.y > _transformed_max.y)
    max.y = _transformed_max.y;

  if ((max.x > min.x) && (max.y > min.y)) {
    if (idx == _model_pos) // TODO ternary
      map->get_model_by_idx(_model_neg)->draw(mvp, min, max, map);
    else
      map->get_model_by_idx(_model_pos)->draw(mvp, min, max, map);
  }
}

bool d3_portal::project(const glm::vec4 &vec, int &x, int &y, const glm::mat4 &mvp) {
  int viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  glm::vec4 cvec = mvp * vec;
  if (cvec.w == 0) {
    x = y = 0;
    return false;
  }
  cvec.x /= cvec.w;
  cvec.y /= cvec.w;
  cvec.x *= 0.5f;
  cvec.y *= 0.5f;
  x = (int)((cvec.x + 0.5f) * (float)viewport[2]);
  y = (int)((cvec.y + 0.5f) * (float)viewport[3]);
  return true;
}

void d3_portal::transform_points(const glm::mat4 &mvp) {
  _transformed_min = glm::ivec2(99999, 99999);
  _transformed_max = glm::ivec2(-99999, -99999);
  for (size_t i = 0; i < _points.size(); ++i)
    if (project(glm::vec4(_points[i], 1.f), _transformed_points[i].x
          , _transformed_points[i].y, mvp)) {
      if (_transformed_points[i].x > _transformed_max.x)
        _transformed_max.x = _transformed_points[i].x;
      if (_transformed_points[i].x < _transformed_min.x)
        _transformed_min.x = _transformed_points[i].x;
      if (_transformed_points[i].y > _transformed_max.y)
        _transformed_max.y = _transformed_points[i].y;
      if (_transformed_points[i].y < _transformed_min.y)
        _transformed_min.y = _transformed_points[i].y;
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
      d3_model m(name, _models.size());
      m.read_from_file(ifs);
      _models.push_back(std::move(m));
    } else if (t == "interAreaPortals") {
      const int num_models = proc_get_next_int(ifs)
        , num_portals = proc_get_next_int(ifs);
      for (int i = 0; i < num_portals; ++i) {
        d3_portal *portal = new d3_portal;
        portal->read_from_file(ifs, this);
        _portals.push_back(portal);
      }
    } else if (t == "nodes") {
      const int num_nodes = proc_get_next_int(ifs);
      _nodes.resize(num_nodes);
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
        // _nodes.push_back(std::move(node));
        _nodes[i] = node;
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
    if (node->plane.in_front(position)) {
      if (node->positive_child > 0)
        node = &_nodes[node->positive_child];
      else
        return node->positive_child;
    } else {
      if (node->negative_child > 0)
        node = &_nodes[node->negative_child];
      else
        return node->negative_child;
    }
}

d3map::d3map(const std::string &filename)
  : _sp(shaders::simple_vert, shaders::simple_frag) {

  _sp.use_this_prog();
  _vertex_pos_attr = _sp.bind_attrib("vertex_pos");
  _mvp_mat_unif = _sp.bind_uniform("mvp");

  this->_load_proc(filename + ".proc");

  _sp.dont_use_this_prog();
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

d3_model* d3map::get_model_by_idx(int idx) {
  return &_models[idx];
}

void d3map::draw(const glm::vec3 &position, const glm::mat4 &mvp
    , const frustum &f) {
  _sp.use_this_prog();

  glUniformMatrix4fv(_mvp_mat_unif, 1, GL_FALSE, glm::value_ptr(mvp));

  glm::ivec2 min = glm::ivec2(0, 0)
    , max = glm::ivec2(g_screen->get_window_width()
        , g_screen->get_window_height());
#if 0
  int start_model = _get_model_idx_by_pos(position);
  if (start_model < 0)
    _models[-1 - start_model].draw(mvp, min, max, this);
  else // position is in the void
    for (d3_model &m : _models)
      m.draw(mvp, min, max, this);
#else
  // NOTE: I admit that I failed to implement proper portal engine rendering
  // of Doom 3 levels (in given deadlines). Someday I probably will revisit
  // this then-forgotten file and this:
  // https://www.iddevnet.com/doom3/visportals.php
  for (d3_model &m : _models)
    m.draw(mvp, min, max, this);
#endif

  // TODO at the end of draw?
  // glScissor(0, 0, g_screen->get_window_width(), g_screen->get_window_height());

  _sp.dont_use_this_prog();
}

};

