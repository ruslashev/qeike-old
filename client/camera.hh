#pragma once

#include <glm/gtc/matrix_transform.hpp>

class camera {
  float _pos_x, _pos_y, _pos_z, _pitch, _yaw, _roll;
public:
  camera();
  camera(float n_pos_x, float n_pos_y, float n_pos_z);
  void update(float xrel, float yrel);
  glm::mat4 compute_view_mat() const;
};

