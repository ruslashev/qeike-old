#pragma once

#include <glm/gtc/quaternion.hpp>

namespace qkmath {

const float C_PI = 3.14159265358979323846f
  , C_INFINITY = 1e30f
  , C_DEG2RAD  = C_PI / 180.f
  , C_RAD2DEG  = 180.f / C_PI
  , C_FLT_EPSILON = 1.192092896e-07f;

#define DEG2RAD(a) ((a) * qkmath::C_DEG2RAD)
#define RAD2DEG(a) ((a) * qkmath::C_RAD2DEG)

inline float inv_sqrt(float x) {
  return 1.f / sqrtf(x);
}

inline float sqrt(float x) {
  return sqrtf(x);
}

template<class T> inline T max(T x, T y) { return (x > y) ? x : y; }
template<class T> inline T min(T x, T y) { return (x < y) ? x : y; }

inline float fabs(float x) {
  int tmp = *reinterpret_cast<int*>(&x);
  tmp &= 0x7FFFFFFF;
  return *reinterpret_cast<float*>(&tmp);
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

inline void sincos(float a, float &s, float &c) {
  s = sinf(a);
  c = cosf(a);
}

inline float acos(float a) {
  if (a <= -1.0f)
    return C_PI;
  if (a >= 1.0f)
    return 0.0f;
  return acosf(a);
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

};

