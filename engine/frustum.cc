#include "frustum.hh"
#include <glm/gtc/matrix_access.hpp>
#include "math.hh"

using namespace qke;

void frustum::_extract_plane(size_t plane_idx, const glm::mat4 &mvp, int row) {
  const int scale = (row < 0) ? -1 : 1;
  row = abs(row) - 1;

  _planes[plane_idx].normal.x = mvp[0][3] + scale * mvp[0][row];
  _planes[plane_idx].normal.y = mvp[1][3] + scale * mvp[1][row];
  _planes[plane_idx].normal.z = mvp[2][3] + scale * mvp[2][row];
  _planes[plane_idx].distance = mvp[3][3] + scale * mvp[3][row];

  _planes[plane_idx].fast_normalize();
}

void frustum::extract_planes(const glm::mat4 &mvp) {
  _extract_plane(0, mvp,  1);
  _extract_plane(1, mvp, -1);
  _extract_plane(2, mvp,  2);
  _extract_plane(3, mvp, -2);
  _extract_plane(4, mvp,  3);
  _extract_plane(5, mvp, -3);
}

bool frustum::box_in_frustum(const glm::vec3 &min, const glm::vec3 &max) const {
  for (const math::plane &p : _planes) {
    if (p.point_in_front_of_plane(glm::vec3(min.x, min.y, min.z), 0.1f))
      continue;
    if (p.point_in_front_of_plane(glm::vec3(min.x, min.y, max.z), 0.1f))
      continue;
    if (p.point_in_front_of_plane(glm::vec3(min.x, max.y, min.z), 0.1f))
      continue;
    if (p.point_in_front_of_plane(glm::vec3(min.x, max.y, max.z), 0.1f))
      continue;
    if (p.point_in_front_of_plane(glm::vec3(max.x, min.y, min.z), 0.1f))
      continue;
    if (p.point_in_front_of_plane(glm::vec3(max.x, min.y, max.z), 0.1f))
      continue;
    if (p.point_in_front_of_plane(glm::vec3(max.x, max.y, min.z), 0.1f))
      continue;
    if (p.point_in_front_of_plane(glm::vec3(max.x, max.y, max.z), 0.1f))
      continue;
    return false;
  }
  return true;
}

bool frustum::point_in_front_of_plane(int plane_idx, const glm::vec3 &point)
  const {
  return _planes[plane_idx].point_in_front_of_plane(point, 0.1f);
}

