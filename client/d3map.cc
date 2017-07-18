#include "d3map.hh"
#include "utils.hh"
#include "math.hh"
#include "shaders.hh"
#include <fstream>
#include <glm/gtc/type_ptr.hpp>

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
    _points.push_back(tmp);
  }
  _create_plane_from_points();

  std::string name_positive = "_area" + std::to_string(positive_model)
    , name_negative = "_area" + std::to_string(negative_model);
  _model_pos = map->get_model_idx_by_name(name_positive);
  _model_neg = map->get_model_idx_by_name(name_negative);
  if (_model_pos >= 0)
    map->add_portal_to_model(this, _model_pos);
  if (_model_neg >= 0)
    map->add_portal_to_model(this, _model_neg);
}

void d3_portal::_create_plane_from_points() {
  glm::vec3 center = glm::vec3(0);
  for (const glm::vec3 &p : _points)
    center += p;
  center *= (1.0f / _points.size());
  glm::vec3 v1 = _points[0] - center, v2 = _points[1] - center;
  _plane.normal = glm::cross(v2, v1);
  _plane.fit_through_point(_points[0]);
  _plane.fast_normalize();
}

d3_portal::frustum_visibility d3_portal::_check_visibility(const frustum &f)
  const {
  // TODO this is shit
  for (int frustum_idx = 0; frustum_idx < 6; ++frustum_idx) {
    bool l_all_back = true;
    for (const glm::vec3 &p : _points)
      if (f.point_in_front_of_plane(frustum_idx, p)) {
        l_all_back = false;
        break;
      }
    if (l_all_back)
      return FRUSTUM_VISIBILITY_OUTSIDE;
  }
  return FRUSTUM_VISIBILITY_INSIDE;
}

void d3_portal::_transform_points(const glm::mat4 &projection) {
  _transformed_min = glm::ivec2( 99999, 99999);
  _transformed_max = glm::ivec2(-99999,-99999);

  for (const glm::vec3 &p : _points) {
    glm::vec2 projected;
    // TODO
    if (qkmath::project_point(p, projection, 700, 525, projected)) {
      glm::ivec2 projected_int = glm::ivec2(projected.x, projected.y);
      if (projected_int.x > _transformed_max.x)
        _transformed_max.x = projected_int.x;
      if (projected_int.x < _transformed_min.x)
        _transformed_min.x = projected_int.x;

      if (projected_int.y > _transformed_max.y)
        _transformed_max.y = projected_int.y;
      if (projected_int.y < _transformed_min.y)
        _transformed_min.y = projected_int.y;
    }
  }
}

void d3_portal::draw_from_model(int idx, const glm::ivec2 &min
      , const glm::ivec2 &max, const frustum &f
      , const glm::mat4 &projection, d3map *map) {
  frustum_visibility v = _check_visibility(f);
  if (v == FRUSTUM_VISIBILITY_OUTSIDE)
    return;
  _transform_points(projection);

  if (v == FRUSTUM_VISIBILITY_INTERSECTS) {
    _transformed_min = glm::ivec2(0, 0);
    _transformed_max = glm::ivec2(700, 525);
  }

  glm::ivec2 n_min = min, n_max = max;
  if (n_min.x < _transformed_min.x)
    n_min.x = _transformed_min.x;
  if (n_max.x > _transformed_max.x)
    n_max.x = _transformed_max.x;
  if (n_min.y < _transformed_min.y)
    n_min.y = _transformed_min.y;
  if (n_max.y > _transformed_max.y)
    n_max.y = _transformed_max.y;

  if ((n_max.x > n_min.x) && (n_max.y > n_min.y)) {
    if (idx == _model_pos)
      map->get_model_by_idx(_model_neg)->draw(n_min, n_max, f, projection, map);
    else
      map->get_model_by_idx(_model_pos)->draw(n_min, n_max, f, projection, map);
  }
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

void d3_model::draw(const glm::ivec2 &min, const glm::ivec2 &max
    , const frustum &f, const glm::mat4 &projection, d3map *map) {
  // TODO
  // glScissor(min.x, min.y, max.x - min.x + 1, max.y - min.y + 1);
  for (const d3_surface *s : _surfaces)
    s->draw();
  for (d3_portal *p : portals)
    p->draw_from_model(index, min, max, f, projection, map);
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
        node.plane.distance = proc_get_next_float(ifs);
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
    if (node->plane.point_in_front_of_plane(position, 0.1f)) {
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

void d3map::draw(const glm::mat4 &projection, const glm::mat4 &view
    , const frustum &f) {
  _sp.use_this_prog();

  glUniformMatrix4fv(_mvp_mat_unif, 1, GL_FALSE
      , glm::value_ptr(projection * view));

  int start_model = _get_model_idx_by_pos(f.position);
  if (start_model >= 0)
    // position is in the void
    for (d3_model &m : _models)
      m.draw(glm::ivec2(0, 0), glm::ivec2(700, 525), f, projection, this);
  else
    _models[-1 - start_model].draw(glm::ivec2(0, 0), glm::ivec2(700, 525), f
        , projection, this);

  _sp.dont_use_this_prog();
}

