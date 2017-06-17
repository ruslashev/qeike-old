#pragma once

namespace qkmath {

inline float fastinvsqrt(float x) {
  float xhalf = x * 0.5f;
  int i = *(int*)&x;
  i = 0x5f3759df - (i >> 1);
  x = *(float*)&i;
  x *= 1.5f - xhalf * x * x;
  return x;
}

inline void quat_get_angleaxis(const glm::quat &q, float *angle, glm::vec3 *axis) {
  const float sq_len = q.x * q.x + q.y * q.y + q.z * q.z;
  if (sq_len > 1e-5f * 1e-5f) {
    *angle = 2.0f * (float)acos(q.w);
    const float inv_len = fastinvsqrt(sq_len);
    *axis = glm::vec3(q.x * inv_len, q.y * inv_len, q.z * inv_len);
  } else {
    *angle = 0;
    *axis = glm::vec3(1, 0, 0);
  }
}

};

