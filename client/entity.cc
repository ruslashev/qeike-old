#include "entity.hh"

entity::entity()
  : entity(glm::vec3(0, 0, 0)) {
}

entity::entity(float n_pos_x, float n_pos_y, float n_pos_z)
  : entity(glm::vec3(n_pos_x, n_pos_y, n_pos_z)) {
}

entity::entity(glm::vec3 n_pos)
  : pitch(0)
  , yaw(0)
  , roll(0)
  , pos(n_pos) {
}

