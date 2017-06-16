#pragma once

#include <glm/gtc/matrix_transform.hpp>

class entity {
public:
  float pitch, yaw, roll;
  glm::vec3 pos;

  entity();
  entity(float n_pos_x, float n_pos_y, float n_pos_z);
  entity(glm::vec3 n_pos);
};

