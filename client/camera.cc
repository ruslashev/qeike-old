#include "camera.hh"
#include "utils.hh"
#include "config.hh"

camera::camera(entity *n_entity_attached_to)
  : entity_attached_to(n_entity_attached_to) {
}

void camera::update_view_angles(float xrel, float yrel) {
  entity_attached_to->yaw += xrel * config.sensitivity * config.m_yaw;
  entity_attached_to->pitch -= yrel * config.sensitivity * config.m_pitch;

  if (entity_attached_to->yaw >= 360.f)
    entity_attached_to->yaw -= 360.f;
  if (entity_attached_to->yaw < 0.f)
    entity_attached_to->yaw += 360.f;

  if (entity_attached_to->pitch >= config.pitch_max)
    entity_attached_to->pitch = config.pitch_max;
  if (entity_attached_to->pitch <= -config.pitch_max)
    entity_attached_to->pitch = -config.pitch_max;
}

void camera::update_position(double dt, int move, int strafe) {
  const float pitch_rad = to_radians(entity_attached_to->pitch)
    , yaw_rad = to_radians(entity_attached_to->yaw)
    , perp_yaw = yaw_rad + static_cast<float>(M_PI_2);
  float dist_move = 5.f * static_cast<float>(move * dt)
    , dist_strafe = 5.f * static_cast<float>(strafe * dt);
  if (move == strafe) {
    dist_move /= sqrtf(2);
    dist_strafe /= sqrtf(2);
  }
  entity_attached_to->pos.x += cosf(yaw_rad) * cosf(pitch_rad) * dist_move
    + cosf(perp_yaw) * dist_strafe;
  entity_attached_to->pos.y += sinf(pitch_rad) * dist_move;
  entity_attached_to->pos.z += sinf(yaw_rad) * cosf(pitch_rad) * dist_move
    + sinf(perp_yaw) * dist_strafe;
}

glm::mat4 camera::compute_view_mat() const {
  const float pitch_rad = to_radians(entity_attached_to->pitch)
    , yaw_rad = to_radians(entity_attached_to->yaw);
  const glm::vec3 view_dir = glm::vec3(cos(yaw_rad) * cos(pitch_rad)
      , sin(pitch_rad), sin(yaw_rad) * cos(pitch_rad))
    , look_at = entity_attached_to->pos + view_dir;
  return glm::lookAt(entity_attached_to->pos, look_at, glm::vec3(0, 1, 0));
}

