#pragma once

#include <glm/gtc/matrix_transform.hpp>

class camera {
  float _pitch, _yaw, _roll;
public:
  glm::vec3 pos;
  camera();
  camera(float n_pos_x, float n_pos_y, float n_pos_z);
  camera(glm::vec3 n_pos);
  void update_view_angles(float xrel, float yrel);
  void update_position(double dt, int move, int strafe);
  glm::mat4 compute_view_mat() const;
};

