#pragma once

#include <glm/gtc/quaternion.hpp>

namespace qkmath {

inline bool float_eq(float a, float b, float epsilon = 1e-5) {
  return (std::abs(a - b) < epsilon);
}

inline float fastinvsqrt(float x) {
  float xhalf = x * 0.5f;
  int i = *(int*)&x;
  i = 0x5f3759df - (i >> 1);
  x = *(float*)&i;
  x *= 1.5f - xhalf * x * x;
  return x;
}

inline void quat_get_angleaxis(const glm::quat &q, float *angle
    , glm::vec3 *axis) {
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

enum plane_side {
  PLANE_SIDE_IN_FRONT,
  PLANE_SIDE_BEHIND,
  PLANE_SIDE_INTERSECTS,
};

class plane {
public:
  glm::vec3 normal;
  float distance;

  void fast_normalize();
  void fit_through_point(const glm::vec3 &point);
  float distance_to_point(const glm::vec3 &point) const;
  plane_side plane_side_point(const glm::vec3 &point, const float epsilon) const;
  bool point_in_front_of_plane(const glm::vec3 &point, const float epsilon)
    const;
};

// TODO may be broken
inline bool project_point(const glm::vec3 &point, const glm::mat4 &projection
    , const int viewport_width, const int viewport_height
    , glm::vec2 &projected_point) {
  return project_point(glm::vec4(point, 1.f), projection, viewport_width
      , viewport_height, projected_point);
}

// TODO may be broken
inline bool project_point(const glm::vec4 &point, const glm::mat4 &projection
    , const int viewport_width, const int viewport_height
    , glm::vec2 &projected_point) {
  glm::vec4 projected_point4 = projection * point;

  if (float_eq(projected_point4.w, 0, 0.01)) {
    projected_point.x = 0;
    projected_point.y = 0;
    return false;
  }

  projected_point4.x /= projected_point4.w;
  projected_point4.y /= projected_point4.w;

  projected_point4.x *= 0.5f;
  projected_point4.y *= 0.5f;

  projected_point.x = (projected_point4.x + 0.5f) * (float)viewport_width;
  projected_point.y = (projected_point4.y + 0.5f) * (float)viewport_height;

  return true;
}

};

