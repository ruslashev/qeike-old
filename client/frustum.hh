#pragma once

#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include "math.hh"

enum frustum_plane {
  FRUSTUM_PLANE_LEFT = 0,
  FRUSTUM_PLANE_RIGHT,
  FRUSTUM_PLANE_UP,
  FRUSTUM_PLANE_DOWN,
  FRUSTUM_PLANE_NEAR,
  FRUSTUM_PLANE_FAR
};

class frustum {
  qkmath::plane _planes[6];

  void _extract_plane(size_t plane_idx, const glm::mat4 &mvp, int row);
public:
  glm::vec3 position;

  void extract_planes(const glm::mat4 &mvp);
  bool box_in_frustum(const glm::vec3 &min, const glm::vec3 &max) const;
  bool point_in_front_of_plane(int plane_idx, const glm::vec3 &point) const;
};

