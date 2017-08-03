#include "math.hh"

using namespace qke::math;

void plane::fast_normalize() {
  const float normal_inverse_length = fastinvsqrt(normal.x * normal.x
      + normal.y * normal.y + normal.z * normal.z);
  normal.x *= normal_inverse_length;
  normal.y *= normal_inverse_length;
  normal.z *= normal_inverse_length;
  distance *= normal_inverse_length;
}

void plane::fit_through_point(const glm::vec3 &point) {
  distance = -glm::dot(normal, point);
}

float plane::distance_to_point(const glm::vec3 &point) const {
  return glm::dot(normal, point) + distance;
}

plane::side plane::plane_side_point(const glm::vec3 &point
    , const float epsilon) const {
  const float d = distance_to_point(point);
  if (d > epsilon)
    return side::in_front;
  if (d < -epsilon)
    return side::behind;
  return side::intersects;
}

bool plane::point_in_front_of_plane(const glm::vec3 &point
    , const float epsilon) const {
  return plane_side_point(point, epsilon) == side::in_front;
}

