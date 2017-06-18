#pragma once

#include "entity.hh"
#include <glm/gtc/matrix_transform.hpp>

class camera {
  // entity *_entity_attached_to;
public:
  glm::vec3 pos;
  float pitch, yaw, roll;
  camera(float n_pos_x, float n_pos_y, float n_pos_z);
  void update_view_angles(float xrel, float yrel);
  void update_position(double dt, int move, int strafe);
  glm::mat4 compute_view_mat() const;
};

