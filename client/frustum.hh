#pragma once

#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>

class frustum {
  glm::vec4 _planes[6];
  void _extract_plane(size_t plane_idx, const glm::mat4 &mvp, int row);
public:
  void extract_planes(const glm::mat4 &mvp);
  bool box_in_frustum(const glm::vec3 &min, const glm::vec3 &max) const;
};

