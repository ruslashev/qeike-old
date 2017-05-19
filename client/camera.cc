#include "camera.hh"
#include "utils.hh"
#include "config.hh"

camera::camera()
  : camera(0, 0, 0) {
}

camera::camera(float n_pos_x, float n_pos_y, float n_pos_z)
  : _pos_x(n_pos_x)
  , _pos_y(n_pos_y)
  , _pos_z(n_pos_z)
  , _pitch(0)
  , _yaw(0)
  , _roll(0) {
}

void camera::update(float xrel, float yrel) {
  _yaw -= xrel * config.sensitivity * config.m_yaw;
  _pitch += yrel * config.sensitivity * config.m_pitch;

  if (_yaw >= 360.f)
    _yaw -= 360.f;
  if (_yaw < 0.f)
    _yaw += 360.f;

  if (_pitch >= config.pitch_max)
    _pitch = config.pitch_max;
  if (_pitch <= -config.pitch_max)
    _pitch = -config.pitch_max;
}

glm::mat4 camera::compute_view_mat() const {
  const float pitch_rad = to_radians(_pitch), yaw_rad = to_radians(_yaw);
  const glm::vec3 position = glm::vec3(_pos_x, _pos_y, _pos_z)
    , view_dir = glm::vec3(cos(yaw_rad) * cos(pitch_rad), sin(pitch_rad)
        , sin(yaw_rad) * cos(pitch_rad))
    , look_at = position + view_dir;
  return glm::lookAt(position, look_at, glm::vec3(0, 1, 0));
}

