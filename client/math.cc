#include "math.hh"

void qkmath::plane::fast_normalize() {
  const float normal_inverse_length = qkmath::fastinvsqrt(normal.x * normal.x
      + normal.y * normal.y + normal.z * normal.z);
  normal.x *= normal_inverse_length;
  normal.y *= normal_inverse_length;
  normal.z *= normal_inverse_length;
  distance *= normal_inverse_length;
}

void qkmath::plane::fit_through_point(const glm::vec3 &point) {
  distance = -glm::dot(normal, point);
}

float qkmath::plane::distance_to_point(const glm::vec3 &point) const {
  return glm::dot(normal, point) + distance;
}

qkmath::plane_side qkmath::plane::plane_side_point(const glm::vec3 &point
    , const float epsilon) const {
  const float d = distance_to_point(point);
  if (d > epsilon)
    return PLANE_SIDE_IN_FRONT;
  if (d < -epsilon)
    return PLANE_SIDE_BEHIND;
  return PLANE_SIDE_INTERSECTS;
}

bool qkmath::plane::point_in_front_of_plane(const glm::vec3 &point
    , const float epsilon) const {
  return plane_side_point(point, epsilon) == PLANE_SIDE_IN_FRONT;
}

