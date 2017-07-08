#include "d3map.hh"
#include "utils.hh"
#include "shaders.hh"
#include <fstream>
#include <sstream>
#include <glm/gtc/type_ptr.hpp>

void d3_plane::normalize() {
  // TODO fastinvsqrt
  float inv_length = 1.0f / glm::length(normal);
  normal *= inv_length;
  dist *= inv_length;
}

bool d3_plane::in_front(const glm::vec3 &point) {
  return (point.x * normal.x + point.y * normal.y + point.z * normal.z + dist
    > 0); // TODO glm::dot
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
    _points.push_back(std::move(tmp));
  }

  std::stringstream name_pos, name_neg;
  name_pos << "_area" << pos_area;
  _area_pos = map->get_area_idx_by_name(name_pos.str());
  name_neg << "_area" << neg_area;
  _area_neg = map->get_area_idx_by_name(name_neg.str());

  if (_area_pos >= 0)
    map->add_portal_to_area(this, _area_pos);
  if (_area_neg >= 0)
    map->add_portal_to_area(this, _area_neg);
}

void d3_portal::render_from_area(const glm::vec3 &position, int idx) {
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

  // render area if visible
  if ((max.x > min.x) && (max.y > min.y)) {
    if (index == m_area_pos)
      m_scene->get_area(m_area_neg)->render(camera, min, max);
    else
      m_scene->get_area(m_area_pos)->render(camera, min, max);
  }
#endif
}

d3_area::d3_area(const std::string &n_name, int n_index)
  : name(n_name)
  , index(n_index) {
}

void d3_area::read_from_file(std::ifstream &file, d3map *f) {
  const int num_surfaces = proc_get_next_int(file);
  for (int i = 0; i < num_surfaces; ++i) {
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
      f->vertices.push_back(std::move(v));
    }
    for (int j = 0; j < num_ind; ++j)
      f->elements.push_back(f->ebo_size + proc_get_next_int(file));

    f->ebo_size += num_verts;
  }
}

void d3_area::draw(const glm::vec3 &position) const {
  // for (size_t i = 0; i < _vbos.size(); ++i) {
  //   glDrawElements(GL_TRIANGLES, _ebo_sizes[i], GL_UNSIGNED_INT, 0);
  // }
  // for (d3_portal *p : portals)
  //   p->render_from_area(position, index);
}

void d3map::_load_proc(const std::string &filename) {
  std::ifstream ifs(filename.c_str(), std::ios::binary);
  if (!ifs)
    die("failed to open \"%s\"", filename.c_str());

  ebo_size = 0;
  std::string t;
  while (ifs >> t) {
    if (t == "model") {
      std::string name = proc_get_next_string(ifs);
      d3_area a(name, _areas.size());
      a.read_from_file(ifs, this);
      _areas.push_back(std::move(a));
    } else if (t == "interAreaPortals") {
      const int num_areas = proc_get_next_int(ifs),
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
        node.plane.normal.x = proc_get_next_float(ifs);
        node.plane.normal.y = proc_get_next_float(ifs);
        node.plane.normal.z = proc_get_next_float(ifs);
        node.plane.dist = proc_get_next_float(ifs);
        node.pos_child = proc_get_next_int(ifs);
        node.neg_child = proc_get_next_int(ifs);
        if (node.pos_child < 0) {
          std::stringstream name;
          name << "_area" << (-1 - node.pos_child);
          node.pos_child = -1 - get_area_idx_by_name(name.str());
        }
        if (node.neg_child < 0) {
          std::stringstream name;
          name << "_area" << (-1 - node.neg_child);
          node.neg_child = -1 - get_area_idx_by_name(name.str());
        }
        _nodes.push_back(std::move(node));
      }
    } else {
      // token skipped
    }
  }
  ifs.close();
}

int d3map::_get_area_idx_by_pos(const glm::vec3 &position) {
  if (!_nodes.size())
    return 0;
  d3_node *node = &_nodes[0];
  while (true)
    if (node->plane.in_front(position)) { // in front
      if (node->pos_child > 0)
        node = &_nodes[node->pos_child];
      else
        return node->pos_child;
    } else { // behind
      if (node->neg_child > 0)
        node = &_nodes[node->neg_child];
      else
        return node->neg_child;
    }
}

void d3map::add_portal_to_area(d3_portal *p, int idx) {
  _areas[idx].portals.push_back(p);
}

int d3map::get_area_idx_by_name(const std::string &name) {
  for (size_t i = 0; i < _areas.size(); ++i)
    if (_areas[i].name == name)
      return i;
  return -1;
}

d3map::d3map(const std::string &filename)
  : _sp(shaders::simple_vert, shaders::simple_frag) {

  _sp.use_this_prog();
  _vertex_pos_attr = _sp.bind_attrib("vertex_pos");
  _mvp_mat_unif = _sp.bind_uniform("mvp");

  this->_load_proc(filename + ".proc");

  vbo.bind();
  glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(d3_vertex)
      , vertices.data(), GL_STATIC_DRAW);
  ebo.bind();
  ebo.upload(elements);

  _sp.dont_use_this_prog();
}

void d3map::draw(const glm::vec3 &position, const glm::mat4 &mvp
    , const frustum &f) {
  // _set_visible_faces(position, f);

  _sp.use_this_prog();

  glUniformMatrix4fv(_mvp_mat_unif, 1, GL_FALSE, glm::value_ptr(mvp));

  glEnableVertexAttribArray(_vertex_pos_attr);
  glVertexAttribPointer(_vertex_pos_attr, 3, GL_FLOAT, GL_FALSE
      , sizeof(d3_vertex), 0);

  vbo.bind();
  ebo.bind();
  glDrawElements(GL_TRIANGLES, ebo_size, GL_UNSIGNED_INT, 0);

  // int start_area = _get_area_idx_by_pos(position);
  // if (start_area >= 0) {
    // position is in the void
    // for (const d3_area &a : _areas) {
    //   a.draw(position);
    // }
  // else
  //   _areas[-1 - start_area].draw(position);

  _sp.dont_use_this_prog();
}

