#pragma once

#include "entity.hh"
#include <glm/gtc/matrix_transform.hpp>

class camera {
  entity *entity_attached_to;
public:
  glm::vec3 pos, vel;
  float pitch, yaw, roll;
  camera(entity *n_entity_attached_to);
  void update_view_angles(float xrel, float yrel);
  void update_position(double dt, int move, int strafe);
  glm::vec3 compute_view_dir() const;
  glm::mat4 compute_view_mat() const;
};

