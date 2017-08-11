#pragma once

#include <glm/gtc/quaternion.hpp>

namespace qke {

namespace math { // TODO remove TODO or not?

const float C_PI = 3.1415926535897932384626433832795f
  , C_INFINITY = 1e30f
  , C_PI_DIV_180  = C_PI / 180.f
  , C_180_DIV_PI  = 180.f / C_PI
  , C_FLT_EPSILON = 1.192092896e-07f;

inline float deg2rad(float deg) { return deg * C_PI_DIV_180; }
inline float rad2deg(float rad) { return rad * C_180_DIV_PI; }

// template<class T> inline T max(T x, T y) { return (x > y) ? x : y; }
// template<class T> inline T min(T x, T y) { return (x < y) ? x : y; }

inline float fabs(float x) {
  int i = *(int*)(&x);
  i &= 0x7FFFFFFF;
  return *(float*)(&i);
}

inline bool float_eq(float a, float b, float epsilon = 1e-5) {
  return (fabs(a - b) < epsilon);
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

class plane {
public:
  enum class side {
    in_front,
    behind,
    intersects,
  };

  glm::vec3 normal;
  float distance;

  void fast_normalize();
  void fit_through_point(const glm::vec3 &point);
  float distance_to_point(const glm::vec3 &point) const;
  side plane_side_point(const glm::vec3 &point, const float epsilon) const;
  bool point_in_front_of_plane(const glm::vec3 &point, const float epsilon)
    const;
};

};

};

