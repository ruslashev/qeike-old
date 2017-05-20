#pragma once

#include <glm/gtc/matrix_transform.hpp>

class camera {
  float _pos_x, _pos_y, _pos_z, _pitch, _yaw, _roll;
public:
  camera();
  camera(float n_pos_x, float n_pos_y, float n_pos_z);
  void update_view_angles(float xrel, float yrel);
  void update_position(double dt, int move, int strafe);
  glm::mat4 compute_view_mat() const;
  float pos_x() const;
  float pos_y() const;
  float pos_z() const;
};

