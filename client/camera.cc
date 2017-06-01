#include "camera.hh"
#include "utils.hh"
#include "config.hh"

camera::camera()
  : camera(glm::vec3(0, 0, 0)) {
}

camera::camera(float n_pos_x, float n_pos_y, float n_pos_z)
  : camera(glm::vec3(n_pos_x, n_pos_y, n_pos_z)) {
}

camera::camera(glm::vec3 n_pos)
  : _pos(n_pos)
  , _pitch(0)
  , _yaw(0)
  , _roll(0) {
}

void camera::update_view_angles(float xrel, float yrel) {
  _yaw += xrel * config.sensitivity * config.m_yaw;
  _pitch -= yrel * config.sensitivity * config.m_pitch;

  if (_yaw >= 360.f)
    _yaw -= 360.f;
  if (_yaw < 0.f)
    _yaw += 360.f;

  if (_pitch >= config.pitch_max)
    _pitch = config.pitch_max;
  if (_pitch <= -config.pitch_max)
    _pitch = -config.pitch_max;
}

void camera::update_position(double dt, int move, int strafe) {
  const float pitch_rad = to_radians(_pitch), yaw_rad = to_radians(_yaw)
    , perp_yaw = yaw_rad + static_cast<float>(M_PI_2);
  float dist_move = 5.f * static_cast<float>(move * dt)
    , dist_strafe = 5.f * static_cast<float>(strafe * dt);
  if (move == strafe) {
    dist_move /= sqrtf(2);
    dist_strafe /= sqrtf(2);
  }
  _pos.x += cosf(yaw_rad) * cosf(pitch_rad) * dist_move + cosf(perp_yaw) * dist_strafe;
  _pos.y += sinf(pitch_rad) * dist_move;
  _pos.z += sinf(yaw_rad) * cosf(pitch_rad) * dist_move + sinf(perp_yaw) * dist_strafe;
}

glm::mat4 camera::compute_view_mat() const {
  const float pitch_rad = to_radians(_pitch), yaw_rad = to_radians(_yaw);
  const glm::vec3 view_dir = glm::vec3(cos(yaw_rad) * cos(pitch_rad)
      , sin(pitch_rad)
      , sin(yaw_rad) * cos(pitch_rad))
    , look_at = _pos + view_dir;
  return glm::lookAt(_pos, look_at, glm::vec3(0, 1, 0));
}

glm::vec3 camera::pos() const {
  return _pos;
}

