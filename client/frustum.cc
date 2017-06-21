#include "frustum.hh"
#include <glm/gtc/matrix_access.hpp>
#include "math.hh"

void frustum::extract_planes(const glm::mat4 &mvp) {
  _extract_plane(0, mvp,  1);
  _extract_plane(1, mvp, -1);
  _extract_plane(2, mvp,  2);
  _extract_plane(3, mvp, -2);
  _extract_plane(4, mvp,  3);
  _extract_plane(5, mvp, -3);
}

bool frustum::box_in_frustum(const glm::vec3 &min, const glm::vec3 &max) const {
  for (int i = 0; i < 6; ++i) {
    glm::vec3 plane_xyz = glm::vec3(_planes[i].x, _planes[i].y, _planes[i].z);
    if (glm::dot(plane_xyz, glm::vec3(min.x, min.y, min.z)) + _planes[i].w > 0)
      continue;
    if (glm::dot(plane_xyz, glm::vec3(min.x, min.y, max.z)) + _planes[i].w > 0)
      continue;
    if (glm::dot(plane_xyz, glm::vec3(min.x, max.y, min.z)) + _planes[i].w > 0)
      continue;
    if (glm::dot(plane_xyz, glm::vec3(min.x, max.y, max.z)) + _planes[i].w > 0)
      continue;
    if (glm::dot(plane_xyz, glm::vec3(max.x, min.y, min.z)) + _planes[i].w > 0)
      continue;
    if (glm::dot(plane_xyz, glm::vec3(max.x, min.y, max.z)) + _planes[i].w > 0)
      continue;
    if (glm::dot(plane_xyz, glm::vec3(max.x, max.y, min.z)) + _planes[i].w > 0)
      continue;
    if (glm::dot(plane_xyz, glm::vec3(max.x, max.y, max.z)) + _planes[i].w > 0)
      continue;
    return false;
  }
  return true;
}

void frustum::_extract_plane(size_t plane_idx, const glm::mat4 &mvp, int row) {
  int scale = (row < 0) ? -1 : 1;
  row = abs(row) - 1;

  _planes[plane_idx].x = mvp[0][3] + scale * mvp[0][row];
  _planes[plane_idx].y = mvp[1][3] + scale * mvp[1][row];
  _planes[plane_idx].z = mvp[2][3] + scale * mvp[2][row];
  _planes[plane_idx].w = mvp[3][3] + scale * mvp[3][row];

  float normal_inverse_length = qkmath::fastinvsqrt(
      _planes[plane_idx].x * _planes[plane_idx].x
      + _planes[plane_idx].y * _planes[plane_idx].y
      + _planes[plane_idx].z * _planes[plane_idx].z);

  _planes[plane_idx].x *= normal_inverse_length;
  _planes[plane_idx].y *= normal_inverse_length;
  _planes[plane_idx].z *= normal_inverse_length;
  _planes[plane_idx].w *= normal_inverse_length;
}

