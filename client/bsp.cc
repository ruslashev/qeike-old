#include "bsp.hh"
#include "utils.hh"
#include "shaders.hh"
#include <fstream>
#include <glm/gtc/type_ptr.hpp>

const float feps = 1e-4f, surface_clip_eps = 0.125f
    , grid_wrap_point_epsilon = 0.1f, subdivide_distance = 4
    , point_epsilon = 0.1f, plane_eps = 0.000001f, plane_tri_eps = 0.1f
    , normal_eps = 0.0001f, dist_eps = 0.02f;

const int max_patch_verts = 2048, bogus_range = 65535
    , max_points_on_winding = 96, max_map_bounds = 65535;

bsp::bsp_plane::bsp_plane()
  : normal({ 0.f, 0.f, 0.f })
  , dist(0.0f) {
}

bsp::bsp_vertex bsp::bsp_vertex::operator+(const bsp_vertex &v) const {
  bsp_vertex res;
  res.position = position + v.position;
  res.lightmap = lightmap + v.lightmap;
  return res;
}

bsp::bsp_vertex bsp::bsp_vertex::operator*(float factor) const {
  bsp_vertex res;
  res.position = position * factor;
  res.lightmap = lightmap * factor;
  return res;
}

void bsp::bsp_patchcollide::clear_bounds() {
  bounds[0] = glm::vec3(99999);
  bounds[1] = glm::vec3(-99999);
}

void bsp::bsp_patchcollide::add_point_to_bounds(const glm::vec3 &v) {
  int		i;
  float	val;

  val = v.x;
  if (val < bounds[0].x)
    bounds[0].x = val;
  if (val > bounds[1].x)
    bounds[1].x = val;

  val = v.y;
  if (val < bounds[0].y)
    bounds[0].y = val;
  if (val > bounds[1].y)
    bounds[1].y = val;

  val = v.z;
  if (val < bounds[0].z)
    bounds[0].z = val;
  if (val > bounds[1].z)
    bounds[1].z = val;
}

void bsp::bsp_grid::set_wrap_width() {
  int i, j;
  for (i = 0 ; i < height ; i++ ) {
    for (j = 0 ; j < 3 ; j++ ) {
      float d = points[0][i][j] - points[width-1][i][j];
      if ( d < -grid_wrap_point_epsilon || d > grid_wrap_point_epsilon ) {
        break;
      }
    }
    if ( j != 3 ) {
      break;
    }
  }
  if ( i == height ) {
    wrap_width = true;
  } else {
    wrap_width = false;
  }
}

// needs_subdivision: returns true if the given quadratic curve is not flat
// enough for collision detection purposes
static bool needs_subdivision(glm::vec3 a, glm::vec3 b, glm::vec3 c) {
  glm::vec3 cmid;
  glm::vec3 lmid;
  glm::vec3 delta;
  float		dist;
  int			i;

  // linear midpoint
  lmid = 0.5f*(a + c);

  // exact curve midpoint
  cmid = 0.5f * ( 0.5f*(a + b) + 0.5f*(b + c) );

  // see if the curve is far enough away from the linear mid
  delta = cmid - lmid;
  dist = glm::length( delta );

  return dist >= subdivide_distance;
}

/// subdivide: a, b, and c are control points.
/// the subdivided sequence will be: a, out1, out2, out3, c
static void subdivide( const glm::vec3 & a, const glm::vec3 & b
    , const glm::vec3 & c, glm::vec3 & out1, glm::vec3 & out2
    , glm::vec3 & out3 ) {
  out1 = 0.5f * (a + b);
  out3 = 0.5f * (b + c);
  out2 = 0.5f * (out1 + out3);
}

static bool compare_points( glm::vec3 a, glm::vec3 b ) {
  float		d;
  d = a.x - b.x;
  if ( d < -point_epsilon || d > point_epsilon ) {
    return false;
  }
  d = a.y - b.y;
  if ( d < -point_epsilon || d > point_epsilon ) {
    return false;
  }
  d = a.z - b.z;
  if ( d < -point_epsilon || d > point_epsilon ) {
    return false;
  }
  return true;
}

void bsp::bsp_grid::subdivide_columns() {
  int		i, j, k;

  for ( i = 0 ; i < width - 2 ;  ) {
    // points[i][x] is an interpolating control point
    // points[i+1][x] is an aproximating control point
    // points[i+2][x] is an interpolating control point

    // first see if we can collapse the aproximating collumn away
    for ( j = 0 ; j < height ; j++ ) {
      if ( needs_subdivision( points[i][j], points[i+1][j], points[i+2][j] ) ) {
        break;
      }
    }
    if ( j == height ) {
      // all of the points were close enough to the linear midpoints
      // that we can collapse the entire column away
      for ( j = 0 ; j < height ; j++ ) {
        // remove the column
        for ( k = i + 2 ; k < width ; k++ ) {
          points[k-1][j] = points[k][j];
        }
      }

      width--;

      // go to the next curve segment
      i++;
      continue;
    }

    //
    // we need to subdivide the curve
    //
    for ( j = 0 ; j < height ; j++ ) {
      glm::vec3	prev, mid, next;

      // save the control points now
      prev = points[i][j];
      mid = points[i+1][j];
      next = points[i+2][j];

      // make room for two additional columns in the grid
      // columns i+1 will be replaced, column i+2 will become i+4
      // i+1, i+2, and i+3 will be generated
      for ( k = width - 1 ; k > i + 1 ; k-- ) {
        points[k+2][j] = points[k][j];
      }

      // generate the subdivided points
      subdivide( prev, mid, next, points[i+1][j], points[i+2][j], points[i+3][j] );
    }

    width += 2;

    // the new aproximating point at i+1 may need to be removed
    // or subdivided farther, so don't advance i
  }
}

void bsp::bsp_grid::remove_degenerate_columns() {
  int		i, j, k;

  for ( i = 0 ; i < width - 1 ; i++ ) {
    for ( j = 0 ; j < height ; j++ ) {
      if ( !compare_points( points[i][j], points[i+1][j] ) ) {
        break;
      }
    }

    if ( j != height ) {
      continue;	// not degenerate
    }

    for ( j = 0 ; j < height ; j++ ) {
      // remove the column
      for ( k = i + 2 ; k < width ; k++ ) {
        points[k-1][j] =  points[k][j];
      }
    }
    width--;

    // check against the next column
    i--;
  }
}

void bsp::bsp_grid::transpose() {
  int			i, j, l;
  glm::vec3 temp;
  bool	temp_wrap;

  if ( width > height ) {
    for ( i = 0 ; i < height ; i++ ) {
      for ( j = i + 1 ; j < width ; j++ ) {
        if ( j < height ) {
          // TODO: std::swap
          // swap the value
          temp =  points[i][j];
          points[i][j] =  points[j][i];
          points[j][i] =  temp;
        } else {
          // just copy
          points[i][j] =  points[j][i];
        }
      }
    }
  } else {
    for ( i = 0 ; i < width ; i++ ) {
      for ( j = i + 1 ; j < height ; j++ ) {
        if ( j < width ) {
          // swap the value
          temp =  points[j][i];
          points[j][i] =  points[i][j];
          points[i][j] =  temp;
        } else {
          // just copy
          points[j][i] =  points[i][j];
        }
      }
    }
  }

  l = width;
  width = height;
  height = l;

  temp_wrap = wrap_width;
  wrap_width = wrap_height;
  wrap_height = temp_wrap;
}

static bool plane_from_points( glm::vec4 &plane, glm::vec3 a, glm::vec3 b, glm::vec3 c ) {
  glm::vec3	d1, d2;
  d1 = b - a;
  d2 = c - a;
  glm::vec3 tplane = glm::cross(d2, d1);
  if (glm::length(tplane) < plane_eps)
    return false;
  tplane = glm::normalize(tplane);
  plane.x = tplane.x;
  plane.y = tplane.y;
  plane.z = tplane.z;
  plane.w = glm::dot(a, tplane);
  return true;
}

static int sign_bits_for_normal( glm::vec3 normal ) {
  int bits = 0;
  if (normal.x < 0)
    bits |= 1<<0;
  if (normal.y < 0)
    bits |= 1<<1;
  if (normal.z < 0)
    bits |= 1<<2;
  return bits;
}

int bsp::find_plane( glm::vec3 p1, glm::vec3 p2, glm::vec3 p3 ) {
  glm::vec4	plane;
  int		i;
  float	d;

  if ( !plane_from_points( plane, p1, p2, p3 ) ) {
    return -1;
  }

  // see if the points are close enough to an existing plane
  for ( i = 0 ; i < numPlanes ; i++ ) {
    if ( glm::dot( glm::vec3(plane), glm::vec3(planes[i].plane) ) < 0 ) {
      continue;	// allow backwards planes?
    }

    glm::vec3 comp_plane = glm::vec3(planes[i].plane);

    d = glm::dot( p1, comp_plane ) - planes[i].plane[3];
    if ( d < -plane_tri_eps || d > plane_tri_eps ) {
      continue;
    }

    d = glm::dot( p2, comp_plane ) - planes[i].plane[3];
    if ( d < -plane_tri_eps || d > plane_tri_eps ) {
      continue;
    }

    d = glm::dot( p3, comp_plane ) - planes[i].plane[3];
    if ( d < -plane_tri_eps || d > plane_tri_eps ) {
      continue;
    }

    // found it
    return i;
  }

  // add a new plane
  assertf(numPlanes != MAX_PATCH_PLANES, "things that shouldn't happen for 400");

  planes[numPlanes].plane = plane;
  planes[numPlanes].sign_bits = sign_bits_for_normal( plane );

  numPlanes++;

  return numPlanes-1;
}

static int grid_plane( int gridPlanes[MAX_GRID_SIZE][MAX_GRID_SIZE][2], int i, int j, int tri ) {
  int		p;
  p = gridPlanes[i][j][tri];
  if ( p != -1 ) {
    return p;
  }
  p = gridPlanes[i][j][!tri];
  if ( p != -1 ) {
    return p;
  }
  // should never happen
  warning_ln("grid plane unresolvable\n" );
  return -1;
}

int bsp::_edge_plane_for_num(bsp_grid *grid, int gridPlanes[MAX_GRID_SIZE][MAX_GRID_SIZE][2], int i, int j, int k ) {
  glm::vec3	p1, p2;
  glm::vec3		up;
  int			p;

  switch ( k ) {
    case 0:	// top border
      p1 = grid->points[i][j];
      p2 = grid->points[i+1][j];
      p = grid_plane( gridPlanes, i, j, 0 );
      up =  p1 +  glm::vec3(planes[ p ].plane) *  4.f;
      return find_plane( p1, p2, up );

    case 2:	// bottom border
      p1 = grid->points[i][j+1];
      p2 = grid->points[i+1][j+1];
      p = grid_plane( gridPlanes, i, j, 1 );
      up =  p1 +  glm::vec3(planes[ p ].plane) *  4.f;
      return find_plane( p2, p1, up );

    case 3: // left border
      p1 = grid->points[i][j];
      p2 = grid->points[i][j+1];
      p = grid_plane( gridPlanes, i, j, 1 );
      up =  p1 +  glm::vec3(planes[ p ].plane) *  4.f;
      return find_plane( p2, p1, up );

    case 1:	// right border
      p1 = grid->points[i+1][j];
      p2 = grid->points[i+1][j+1];
      p = grid_plane( gridPlanes, i, j, 0 );
      up =  p1 +  glm::vec3(planes[ p ].plane) *  4.f;
      return find_plane( p1, p2, up );

    case 4:	// diagonal out of triangle 0
      p1 = grid->points[i+1][j+1];
      p2 = grid->points[i][j];
      p = grid_plane( gridPlanes, i, j, 0 );
      up =  p1 +  glm::vec3(planes[ p ].plane) *  4.f;
      return find_plane( p1, p2, up );

    case 5:	// diagonal out of triangle 1
      p1 = grid->points[i][j];
      p2 = grid->points[i+1][j+1];
      p = grid_plane( gridPlanes, i, j, 1 );
      up =  p1 +  glm::vec3(planes[ p ].plane) *  4.f;
      return find_plane( p1, p2, up );
    default:
      warning_ln("bad k");
      return -1;
  }
}

#define	SIDE_FRONT		0
#define	SIDE_ON			2
#define	SIDE_BACK		1
#define	SIDE_CROSS		-2

int bsp::point_on_plane_side( glm::vec3 p, int planeNum ) {
  glm::vec4 plane;
  float	d;

  if ( planeNum == -1 ) {
    return SIDE_ON;
  }
  plane = planes[ planeNum ].plane;

  d = glm::dot( p, glm::vec3(plane) ) - plane[3];

  if ( d > plane_tri_eps ) {
    return SIDE_FRONT;
  }

  if ( d < -plane_tri_eps ) {
    return SIDE_BACK;
  }

  return SIDE_ON;
}

void bsp::set_border_inward( bsp_facet *facet, bsp_grid *grid, int gridPlanes[MAX_GRID_SIZE][MAX_GRID_SIZE][2], int i, int j, int which ) {
  int		k, l;
  glm::vec3 points[4];
  int		numPoints;

  switch ( which ) {
    case -1:
      points[0] = grid->points[i][j];
      points[1] = grid->points[i+1][j];
      points[2] = grid->points[i+1][j+1];
      points[3] = grid->points[i][j+1];
      numPoints = 4;
      break;
    case 0:
      points[0] = grid->points[i][j];
      points[1] = grid->points[i+1][j];
      points[2] = grid->points[i+1][j+1];
      numPoints = 3;
      break;
    case 1:
      points[0] = grid->points[i+1][j+1];
      points[1] = grid->points[i][j+1];
      points[2] = grid->points[i][j];
      numPoints = 3;
      break;
    default:
      warning_ln("bad parameter");
      numPoints = 0;
      break;
  }

  for ( k = 0 ; k < facet->n_borders ; k++ ) {
    int		front, back;

    front = 0;
    back = 0;

    for ( l = 0 ; l < numPoints ; l++ ) {
      int		side;

      side = point_on_plane_side( points[l], facet->border_planes[k] );
      if ( side == SIDE_FRONT ) {
        front++;
      } if ( side == SIDE_BACK ) {
        back++;
      }
    }

    if ( front && !back ) {
      facet->border_inward[k] = true;
    } else if ( back && !front ) {
      facet->border_inward[k] = false;
    } else if ( !front && !back ) {
      // flat side border
      facet->border_planes[k] = -1;
    } else {
      // bisecting side border
      warning_ln( "CM_SetBorderInward: mixed plane sides\n" );
      facet->border_inward[k] = false;
    }
  }
}

bsp::bsp_winding *bsp::base_winding_for_plane (glm::vec3 normal, float dist)
{
	int		i, x;
	float	max, v;
  glm::vec3	org, vright, vup;
	bsp_winding	*w = new bsp_winding;
	
// find the major axis

	max = -bogus_range;
	x = -1;
	for (i=0 ; i<3; i++)
	{
		v = fabs(normal[i]);
		if (v > max)
		{
			x = i;
			max = v;
		}
	}
	assertf (x!=-1, "no axis found");
		
  vup = glm::vec3(0, 0, 0);
	switch (x)
	{
	case 0:
	case 1:
		vup[2] = 1;
		break;		
	case 2:
		vup[0] = 1;
		break;		
	}

	v = glm::dot (vup, normal);
  vup += normal *  (-v);
  vup = glm::normalize(vup);
		
  org = normal * dist;
	
  vright = glm::cross(vup, normal);
	
  vup *= bogus_range;
  vright *= bogus_range;

// project a really big	axis aligned box onto the plane
	
  w->p[0] = org - vright;
  w->p[0] += vup;

  w->p[1] = org + vright;
  w->p[1] += vup;

  w->p[2] = org + vright;
  w->p[2] -= vup;

  w->p[3] = org - vright;
  w->p[3] -= vup;
	
	w->numpoints = 4;
	
	return w;	
}

void bsp::chop_winding_in_place (bsp_winding **inout, glm::vec3 normal, float dist, float epsilon)
{
  bsp_winding *in;
  float	dists[max_points_on_winding+4];
  int sides[max_points_on_winding+4];
  int counts[3];
  //MrElusive: DOH can't use statics when unsing multithreading!!!
  float dot;		// VC 4.2 optimizer bug if not static
  int i, j;
  glm::vec3 p1, p2;
  glm::vec3 mid;
  bsp_winding *f;
  int maxpts;

  in = *inout;
  counts[0] = counts[1] = counts[2] = 0;

  // determine sides for each point
  for (i=0 ; i<in->numpoints ; i++)
  {
    dot = glm::dot (in->p[i], normal);
    dot -= dist;
    dists[i] = dot;
    if (dot > epsilon)
      sides[i] = SIDE_FRONT;
    else if (dot < -epsilon)
      sides[i] = SIDE_BACK;
    else
    {
      sides[i] = SIDE_ON;
    }
    counts[sides[i]]++;
  }
  sides[i] = sides[0];
  dists[i] = dists[0];

  if (!counts[0])
  {
    // FreeWinding (in);
    *inout = NULL;
    return;
  }
  if (!counts[1])
    return;		// inout stays the same

  maxpts = in->numpoints+4;	// cant use counts[0]+2 because
  // of fp grouping errors

  f = new bsp_winding;
  f->numpoints = maxpts;

  for (i=0 ; i<in->numpoints ; i++)
  {
    p1 = in->p[i];

    if (sides[i] == SIDE_ON)
    {
      f->p[f->numpoints] = p1;
      f->numpoints++;
      continue;
    }

    if (sides[i] == SIDE_FRONT)
    {
      f->p[f->numpoints] = p1;
      f->numpoints++;
    }

    if (sides[i+1] == SIDE_ON || sides[i+1] == sides[i])
      continue;

    // generate a split point
    p2 = in->p[(i+1)%in->numpoints];

    dot = dists[i] / (dists[i]-dists[i+1]);
    for (j=0 ; j<3 ; j++)
    {	// avoid round off error when possible
      if (normal[j] == 1)
        mid[j] = dist;
      else if (normal[j] == -1)
        mid[j] = -dist;
      else
        mid[j] = p1[j] + dot*(p2[j]-p1[j]);
    }

    f->p[f->numpoints] = mid;
    f->numpoints++;
  }

  assertf(f->numpoints <= maxpts, "points exceeded estimate");
  assertf(f->numpoints <= max_points_on_winding, "max_points_on_winding");

  // FreeWinding (in);
  *inout = f;
}

void bsp::bsp_winding::bounds (glm::vec3 & mins, glm::vec3 & maxs)
{
	float	v;
	int		i,j;

	mins[0] = mins[1] = mins[2] = 99999;
	maxs[0] = maxs[1] = maxs[2] = -99999;

	for (i=0 ; i<numpoints ; i++)
	{
		for (j=0 ; j<3 ; j++)
		{
			v = p[i][j];
			if (v < mins[j])
				mins[j] = v;
			if (v > maxs[j])
				maxs[j] = v;
		}
	}
}

bool bsp::validate_facet( bsp_facet *facet ) {
  glm::vec4		plane;
  int			j;
  bsp_winding	*w;
  glm::vec3		bounds[2];

  if ( facet->surface_plane == -1 ) {
    return false;
  }

  plane = planes[ facet->surface_plane ].plane;
  w = base_winding_for_plane( plane,  plane[3] );
  for ( j = 0 ; j < facet->n_borders && w ; j++ ) {
    if ( facet->border_planes[j] == -1 ) {
      return false;
    }
    plane = planes[ facet->border_planes[j] ].plane;
    if ( !facet->border_inward[j] ) {
      plane = -plane;
      plane[3] = -plane[3];
    }
    chop_winding_in_place( &w, plane, plane[3], 0.1f );
  }

  if ( !w ) {
    return false;		// winding was completely chopped away
  }

  // see if the facet is unreasonably large
  w->bounds(bounds[0], bounds[1] );
  // FreeWinding( w );

  for ( j = 0 ; j < 3 ; j++ ) {
    if ( bounds[1][j] - bounds[0][j] > max_map_bounds ) {
      return false;		// we must be missing a plane
    }
    if ( bounds[0][j] >= max_map_bounds ) {
      return false;
    }
    if ( bounds[1][j] <= -max_map_bounds ) {
      return false;
    }
  }
  return true;		// winding is fine
}

int bsp::plane_equal(bsp_patch_plane *p, glm::vec4 plane, int *flipped) {
  glm::vec4 invplane;

  if (
      fabs(p->plane[0] - plane[0]) < normal_eps
      && fabs(p->plane[1] - plane[1]) < normal_eps
      && fabs(p->plane[2] - plane[2]) < normal_eps
      && fabs(p->plane[3] - plane[3]) < dist_eps)
  {
    *flipped = false;
    return true;
  }

  invplane = -plane;
  invplane[3] = -plane[3];

  if (
      fabs(p->plane[0] - invplane[0]) < normal_eps
      && fabs(p->plane[1] - invplane[1]) < normal_eps
      && fabs(p->plane[2] - invplane[2]) < normal_eps
      && fabs(p->plane[3] - invplane[3]) < dist_eps)
  {
    *flipped = true;
    return true;
  }

  return false;
}

static void snap_vector(glm::vec3 &normal) {
  int		i;

  for (i=0 ; i<3 ; i++)
  {
    if ( fabs(normal[i] - 1) < normal_eps )
    {
      normal = glm::vec3(0);
      normal[i] = 1;
      break;
    }
    if ( fabs(normal[i] - -1) < normal_eps )
    {
      normal = glm::vec3(0);
      normal[i] = -1;
      break;
    }
  }
}

int bsp::find_plane2(glm::vec4 plane, int *flipped) {
  int i;

  // see if the points are close enough to an existing plane
  for ( i = 0 ; i < numPlanes ; i++ ) {
    if (plane_equal(&planes[i], plane, flipped)) return i;
  }

  // add a new plane
  assertf ( numPlanes != MAX_PATCH_PLANES, "MAX_PATCH_PLANES");

  planes[numPlanes].plane = plane;
  planes[numPlanes].sign_bits = sign_bits_for_normal( plane );

  numPlanes++;

  *flipped = false;

  return numPlanes-1;
}

void bsp::add_facet_bevels( bsp_facet *facet ) {
  int i, j, k, l;
  int axis, dir, order, flipped;
  float d;
  glm::vec4 plane, newplane;
  bsp_winding *w, *w2;
  glm::vec3 mins, maxs, vec, vec2;

  plane = planes[ facet->surface_plane ].plane;

  w = base_winding_for_plane( plane,  plane[3] );
  for ( j = 0 ; j < facet->n_borders && w ; j++ ) {
    if (facet->border_planes[j] == facet->surface_plane) continue;
    plane = planes[ facet->border_planes[j] ].plane;

    if ( !facet->border_inward[j] ) {
      plane = -plane;
      plane[3] = -plane[3];
    }

    chop_winding_in_place( &w, plane, plane[3], 0.1f );
  }
  if ( !w ) {
    return;
  }

  w->bounds(mins, maxs);

  // add the axial planes
  order = 0;
  for ( axis = 0 ; axis < 3 ; axis++ )
  {
    for ( dir = -1 ; dir <= 1 ; dir += 2, order++ )
    {
      plane = glm::vec4(0, 0, 0, plane.w);
      plane[axis] = dir;
      if (dir == 1) {
        plane[3] = maxs[axis];
      }
      else {
        plane[3] = -mins[axis];
      }
      //if it's the surface plane
      if (plane_equal(&planes[facet->surface_plane], plane, &flipped)) {
        continue;
      }
      // see if the plane is allready present
      for ( i = 0 ; i < facet->n_borders ; i++ ) {
        if (plane_equal(&planes[facet->border_planes[i]], plane, &flipped))
          break;
      }

      if ( i == facet->n_borders ) {
        if (facet->n_borders > 4 + 6 + 16)
          warning_ln("too many bevels");
        facet->border_planes[facet->n_borders] = find_plane2(plane, &flipped);
        facet->border_no_adjust[facet->n_borders] = 0;
        facet->border_inward[facet->n_borders] = flipped;
        facet->n_borders++;
      }
    }
  }
  //
  // add the edge bevels
  //
  // test the non-axial plane edges
  for ( j = 0 ; j < w->numpoints ; j++ )
  {
    k = (j+1)%w->numpoints;
    vec = w->p[j] - w->p[k];
    //if it's a degenerate edge
    if (glm::length (vec) < 0.5)
      continue;
    vec = glm::normalize(vec);
    snap_vector(vec);
    for ( k = 0; k < 3 ; k++ )
      if ( vec[k] == -1 || vec[k] == 1 )
        break;	// axial
    if ( k < 3 )
      continue;	// only test non-axial edges

    // try the six possible slanted axials from this edge
    for ( axis = 0 ; axis < 3 ; axis++ )
    {
      for ( dir = -1 ; dir <= 1 ; dir += 2 )
      {
        // construct a plane
        vec2 = glm::vec3(0);
        vec2[axis] = dir;
        plane = glm::vec4(glm::cross(vec, vec2), plane.w);
        if (glm::length (plane) < 0.5)
          continue;
        plane = glm::normalize(plane);
        plane[3] = glm::dot (w->p[j], glm::vec3(plane));

        // if all the points of the facet winding are
        // behind this plane, it is a proper edge bevel
        for ( l = 0 ; l < w->numpoints ; l++ )
        {
          d = glm::dot (w->p[l], glm::vec3(plane)) - plane[3];
          if (d > 0.1)
            break;	// point in front
        }
        if ( l < w->numpoints )
          continue;

        //if it's the surface plane
        if (plane_equal(&planes[facet->surface_plane], plane, &flipped)) {
          continue;
        }
        // see if the plane is allready present
        for ( i = 0 ; i < facet->n_borders ; i++ ) {
          if (plane_equal(&planes[facet->border_planes[i]], plane, &flipped)) {
            break;
          }
        }

        if ( i == facet->n_borders ) {
          if (facet->n_borders > 4 + 6 + 16) warning_ln("too many bevels");
          facet->border_planes[facet->n_borders] = find_plane2(plane, &flipped);

          for ( k = 0 ; k < facet->n_borders ; k++ ) {
            if (facet->border_planes[facet->n_borders] ==
                facet->border_planes[k]) warning_ln("bevel plane already used");
          }

          facet->border_no_adjust[facet->n_borders] = 0;
          facet->border_inward[facet->n_borders] = flipped;
          //
          w2 = new bsp_winding;
          w2->numpoints = w->numpoints;
          w2->p[0] = w->p[0];
          w2->p[1] = w->p[1];
          w2->p[2] = w->p[2];
          w2->p[3] = w->p[3];
          newplane = planes[facet->border_planes[facet->n_borders]].plane;
          if (!facet->border_inward[facet->n_borders])
          {
            newplane = -newplane;
            newplane[3] = -newplane[3];
          } //end if
          chop_winding_in_place( &w2, newplane, newplane[3], 0.1f );
          if (!w2) {
            warning_ln("invalid bevel");
            continue;
          }
          else {
            // FreeWinding(w2);
          }
          //
          facet->n_borders++;
          //already got a bevel
          //					break;
        }
      }
    }
  }
  // FreeWinding( w );

  //add opposite plane
  facet->border_planes[facet->n_borders] = facet->surface_plane;
  facet->border_no_adjust[facet->n_borders] = 0;
  facet->border_inward[facet->n_borders] = true;
  facet->n_borders++;

}

void bsp::_create_patch_collide_from_grid(bsp_grid *grid, bsp::bsp_patchcollide *pf) {
  int				i, j;
  glm::vec3			p1, p2, p3;
  int				gridPlanes[MAX_GRID_SIZE][MAX_GRID_SIZE][2];
  bsp_facet			*facet;
  int				borders[4];
  int				noAdjust[4];

  numPlanes = 0;
  numFacets = 0;

  // find the planes for each triangle of the grid
  for ( i = 0 ; i < grid->width - 1 ; i++ ) {
    for ( j = 0 ; j < grid->height - 1 ; j++ ) {
      p1 = grid->points[i][j];
      p2 = grid->points[i+1][j];
      p3 = grid->points[i+1][j+1];
      gridPlanes[i][j][0] = find_plane( p1, p2, p3 );

      p1 = grid->points[i+1][j+1];
      p2 = grid->points[i][j+1];
      p3 = grid->points[i][j];
      gridPlanes[i][j][1] = find_plane( p1, p2, p3 );
    }
  }

  // create the borders for each facet
  for ( i = 0 ; i < grid->width - 1 ; i++ ) {
    for ( j = 0 ; j < grid->height - 1 ; j++ ) {

      borders[EN_TOP] = -1;
      if ( j > 0 ) {
        borders[EN_TOP] = gridPlanes[i][j-1][1];
      } else if ( grid->wrap_height ) {
        borders[EN_TOP] = gridPlanes[i][grid->height-2][1];
      } 
      noAdjust[EN_TOP] = ( borders[EN_TOP] == gridPlanes[i][j][0] );
      if ( borders[EN_TOP] == -1 || noAdjust[EN_TOP] ) {
        borders[EN_TOP] = _edge_plane_for_num( grid, gridPlanes, i, j, 0 );
      }

      borders[EN_BOTTOM] = -1;
      if ( j < grid->height - 2 ) {
        borders[EN_BOTTOM] = gridPlanes[i][j+1][0];
      } else if ( grid->wrap_height ) {
        borders[EN_BOTTOM] = gridPlanes[i][0][0];
      }
      noAdjust[EN_BOTTOM] = ( borders[EN_BOTTOM] == gridPlanes[i][j][1] );
      if ( borders[EN_BOTTOM] == -1 || noAdjust[EN_BOTTOM] ) {
        borders[EN_BOTTOM] = _edge_plane_for_num( grid, gridPlanes, i, j, 2 );
      }

      borders[EN_LEFT] = -1;
      if ( i > 0 ) {
        borders[EN_LEFT] = gridPlanes[i-1][j][0];
      } else if ( grid->wrap_width ) {
        borders[EN_LEFT] = gridPlanes[grid->width-2][j][0];
      }
      noAdjust[EN_LEFT] = ( borders[EN_LEFT] == gridPlanes[i][j][1] );
      if ( borders[EN_LEFT] == -1 || noAdjust[EN_LEFT] ) {
        borders[EN_LEFT] = _edge_plane_for_num( grid, gridPlanes, i, j, 3 );
      }

      borders[EN_RIGHT] = -1;
      if ( i < grid->width - 2 ) {
        borders[EN_RIGHT] = gridPlanes[i+1][j][1];
      } else if ( grid->wrap_width ) {
        borders[EN_RIGHT] = gridPlanes[0][j][1];
      }
      noAdjust[EN_RIGHT] = ( borders[EN_RIGHT] == gridPlanes[i][j][0] );
      if ( borders[EN_RIGHT] == -1 || noAdjust[EN_RIGHT] ) {
        borders[EN_RIGHT] = _edge_plane_for_num( grid, gridPlanes, i, j, 1 );
      }

      assertf(numFacets != MAX_FACETS, "boop that snoot rn");
      facet = &facets[numFacets];
      memset( facet, 0, sizeof( *facet ) );

      if ( gridPlanes[i][j][0] == gridPlanes[i][j][1] ) {
        if ( gridPlanes[i][j][0] == -1 ) {
          continue;		// degenrate
        }
        facet->surface_plane = gridPlanes[i][j][0];
        facet->n_borders = 4;
        facet->border_planes[0] = borders[EN_TOP];
        facet->border_no_adjust[0] = noAdjust[EN_TOP];
        facet->border_planes[1] = borders[EN_RIGHT];
        facet->border_no_adjust[1] = noAdjust[EN_RIGHT];
        facet->border_planes[2] = borders[EN_BOTTOM];
        facet->border_no_adjust[2] = noAdjust[EN_BOTTOM];
        facet->border_planes[3] = borders[EN_LEFT];
        facet->border_no_adjust[3] = noAdjust[EN_LEFT];
        set_border_inward( facet, grid, gridPlanes, i, j, -1 );
        if ( validate_facet( facet ) ) {
          add_facet_bevels( facet );
          numFacets++;
        }
      } else {
        // two seperate triangles
        facet->surface_plane = gridPlanes[i][j][0];
        facet->n_borders = 3;
        facet->border_planes[0] = borders[EN_TOP];
        facet->border_no_adjust[0] = noAdjust[EN_TOP];
        facet->border_planes[1] = borders[EN_RIGHT];
        facet->border_no_adjust[1] = noAdjust[EN_RIGHT];
        facet->border_planes[2] = gridPlanes[i][j][1];
        if ( facet->border_planes[2] == -1 ) {
          facet->border_planes[2] = borders[EN_BOTTOM];
          if ( facet->border_planes[2] == -1 ) {
            facet->border_planes[2] = _edge_plane_for_num( grid, gridPlanes, i, j, 4 );
          }
        }
        set_border_inward( facet, grid, gridPlanes, i, j, 0 );
        if ( validate_facet( facet ) ) {
          add_facet_bevels( facet );
          numFacets++;
        }

        if ( numFacets == MAX_FACETS ) {
          Com_Error( ERR_DROP, "MAX_FACETS" );
        }
        facet = &facets[numFacets];
        Com_Memset( facet, 0, sizeof( *facet ) );

        facet->surfacePlane = gridPlanes[i][j][1];
        facet->n_borders = 3;
        facet->borderPlanes[0] = borders[EN_BOTTOM];
        facet->borderNoAdjust[0] = noAdjust[EN_BOTTOM];
        facet->borderPlanes[1] = borders[EN_LEFT];
        facet->borderNoAdjust[1] = noAdjust[EN_LEFT];
        facet->borderPlanes[2] = gridPlanes[i][j][0];
        if ( facet->borderPlanes[2] == -1 ) {
          facet->borderPlanes[2] = borders[EN_TOP];
          if ( facet->borderPlanes[2] == -1 ) {
            facet->borderPlanes[2] = _edge_plane_for_num( grid, gridPlanes, i, j, 5 );
          }
        }
        CM_SetBorderInward( facet, grid, gridPlanes, i, j, 1 );
        if ( validate_facet( facet ) ) {
          CM_AddFacetBevels( facet );
          numFacets++;
        }
      }
    }
  }

  // copy the results out
  pf->numPlanes = numPlanes;
  pf->numFacets = numFacets;
  pf->facets = Hunk_Alloc( numFacets * sizeof( *pf->facets ), h_high );
  Com_Memcpy( pf->facets, facets, numFacets * sizeof( *pf->facets ) );
  pf->planes = Hunk_Alloc( numPlanes * sizeof( *pf->planes ), h_high );
  Com_Memcpy( pf->planes, planes, numPlanes * sizeof( *pf->planes ) );
}

enum class lump {
  entities = 0,
  shaders,
  planes,
  nodes,
  leaves,
  leaffaces, // Q3: leafsurfaces
  leafbrushes,
  models,
  brushes,
  brushsides,
  vertices,
  meshverts,
  fogs,
  faces,
  lightmaps,
  lightgrid,
  visdata
};

enum class face {
  polygon = 1,
  patch,
  mesh,
  billboard
};

struct bsp_direntry {
  int offset;
  int length;
};

struct bsp_header {
  char magic[4];
  int version;
  bsp_direntry direntries[17];
};

struct bsp_raw_vertex {
  glm::vec3 position;
  glm::vec2 decal;
  glm::vec2 lightmap;
  glm::vec3 normal;
  unsigned char color[4];
};

struct bsp_raw_leaf {
  int cluster;
  int area;
  glm::ivec3 mins;
  glm::ivec3 maxs;
  int leafface;
  int n_leaffaces;
  int leafbrush;
  int n_leafbrushes;
};

template <typename T>
void load_lump(std::ifstream &ifs, const bsp_header *header
    , lump lump_type, std::vector<T> &container) {
  int num_elements = header->direntries[(int)lump_type].length / sizeof(T);
  container.reserve(num_elements);
  ifs.seekg(header->direntries[(int)lump_type].offset);
  for (int i = 0; i < num_elements; ++i) {
    T element;
    ifs.read((char*)&element, sizeof(T));
    container.push_back(element);
  }
}

void bsp::_load_file(const char *filename, float world_scale
    , int tesselation_level) {
  std::ifstream ifs(filename, std::ios::binary);

  bsp_header header;
  ifs.read((char*)&header, sizeof(bsp_header));
  std::string magic = std::string(header.magic, 4);
  if (magic != "IBSP") {
    ifs.close();
    die("\"%s\": invalid magic \"%s\"", filename, magic.c_str());
  }
  if (header.version != 0x2e) {
    ifs.close();
    die("\"%s\": invalid version %d = 0x%x", filename, header.version
        , header.version);
  }

  load_lump(ifs, &header, lump::planes, _planes);
  for (bsp_plane &p : _planes) {
    std::swap(p.normal.y, p.normal.z);
    p.normal.z = -p.normal.z;
    p.dist /= world_scale;
  }

  load_lump(ifs, &header, lump::nodes, _nodes);
  for (bsp_node &n : _nodes) {
    n.mins /= world_scale;
    std::swap(n.mins.y, n.mins.z);
    n.mins.z = -n.mins.z;
    n.maxs /= world_scale;
    std::swap(n.maxs.y, n.maxs.z);
    n.maxs.z = -n.maxs.z;
  }

  std::vector<bsp_raw_leaf> raw_leaves;
  load_lump(ifs, &header, lump::leaves, raw_leaves);
  for (const bsp_raw_leaf &l : raw_leaves) {
    bsp_leaf conv_leaf;
    conv_leaf.cluster = l.cluster;
    conv_leaf.area = l.area;
    conv_leaf.mins = glm::vec3(l.mins.x, l.mins.z, -l.mins.y) / world_scale;
    conv_leaf.maxs = glm::vec3(l.maxs.x, l.maxs.z, -l.maxs.y) / world_scale;
    conv_leaf.leafface = l.leafface;
    conv_leaf.n_leaffaces = l.n_leaffaces;
    conv_leaf.leafbrush = l.leafbrush;
    conv_leaf.n_leafbrushes = l.n_leafbrushes;
    _leaves.push_back(conv_leaf);
  }

  load_lump(ifs, &header, lump::leaffaces, _leaffaces);

  load_lump(ifs, &header, lump::leafbrushes, _leafbrushes);

  load_lump(ifs, &header, lump::brushes, _brushes);

  load_lump(ifs, &header, lump::brushsides, _brushsides);

  std::vector<bsp_raw_vertex> raw_vertices;
  load_lump(ifs, &header, lump::vertices, raw_vertices);
  for (const bsp_raw_vertex &v : raw_vertices) {
    bsp_vertex conv_vertex;
    conv_vertex.position = v.position / world_scale;
    std::swap(conv_vertex.position.y, conv_vertex.position.z);
    conv_vertex.position.z = -conv_vertex.position.z;
    conv_vertex.lightmap = v.lightmap;
    _vertices.push_back(conv_vertex);
  }

  load_lump(ifs, &header, lump::meshverts, _meshverts);

  load_lump(ifs, &header, lump::faces, _faces);
  _visible_faces.resize(_faces.size(), 0);
  _patchcollides.resize(_faces.size());
  _create_patch_collides();

  _create_patches(tesselation_level);

  for (size_t i = 0; i < _meshverts.size(); i += 3) {
    if (i + 2 >= _meshverts.size()) {
      warning("meshvert winding swap routine ended too soon");
      break;
    }
    std::swap(_meshverts[i + 1], _meshverts[i + 2]);
  }

  load_lump(ifs, &header, lump::lightmaps, _lightmaps);
  _lightmap_texture_ids.resize(_lightmaps.size());

  ifs.seekg(header.direntries[(int)lump::visdata].offset);
  ifs.read((char*)&_visdata.n_vecs, sizeof(int));
  ifs.read((char*)&_visdata.sz_vecs, sizeof(int));
  int size = _visdata.n_vecs * _visdata.sz_vecs;
  _visdata.vecs.resize(size);
  ifs.read((char*)&_visdata.vecs[0], size * sizeof(unsigned char));

  ifs.close();
}

void bsp::_create_patches(int tesselation_level) {
  int patch_count = 0
    , patch_size = (tesselation_level + 1) * (tesselation_level + 1)
    , patch_index_size = tesselation_level * tesselation_level * 6;
  for (const bsp_face &f : _faces)
    if (f.type == (int)face::patch)
      patch_count += ((f.size[0] - 1) / 2) * ((f.size[0] - 1) / 2);

  size_t vertex_count = _vertices.size(), meshverts_count = _meshverts.size();
  _vertices.resize(_vertices.size() + patch_count * patch_size);
  _meshverts.resize(_meshverts.size() + patch_count * patch_index_size);
  for (size_t i = 0, voff = vertex_count, eoff = meshverts_count; i < _faces.size(); ++i)
    if (_faces[i].type == (int)face::patch) {
      int dim_x = (_faces[i].size[0] - 1) / 2, dim_y = (_faces[i].size[1] - 1) / 2;
      _faces[i].meshvert = eoff;
      for (int x = 0, n = 0; n < dim_x; x = 2 * (++n))
        for (int y = 0, m = 0; m < dim_y; y = 2 * (++m)) {
          _tesselate(tesselation_level
              , _faces[i].vertex + x + _faces[i].size[0] * y, _faces[i].size[0]
              , voff, eoff);
          voff += patch_size;
          eoff += patch_index_size;
        }
      _faces[i].n_meshverts = eoff - _faces[i].meshvert;
    } else
      for (int j = 0; j < _faces[i].n_meshverts; ++j)
        _meshverts[_faces[i].meshvert + j].offset += _faces[i].vertex;
}

void bsp::_tesselate(int tesselation_level, int control_offset
    , int control_width, int voff, int eoff) {
  bsp_vertex controls[9];
  for (int c = 0, c_idx = 0; c < 3; c++) {
    int pos = c * control_width;
    controls[c_idx++] = _vertices[control_offset + pos];
    controls[c_idx++] = _vertices[control_offset + pos + 1];
    controls[c_idx++] = _vertices[control_offset + pos + 2];
  }

  int L1 = tesselation_level + 1;

  for (int j = 0; j <= tesselation_level; ++j) {
    float a = (float)j / tesselation_level, b = 1.f - a;
    _vertices[voff + j] = controls[0] * b * b + controls[3] * 2 * b * a
      + controls[6] * a * a;
  }

  for (int i = 1; i <= tesselation_level; ++i) {
    float a = (float)i / tesselation_level, b = 1.f - a;
    bsp_vertex temp[3];
    for (int j = 0; j < 3; ++j) {
      int k = 3 * j;
      temp[j] = controls[k] * b * b + controls[k + 1] * 2 * b * a + controls[k + 2] * a * a;
    }
    for (int j = 0; j <= tesselation_level; ++j) {
      float n_a = (float)j / tesselation_level, n_b = 1.f - n_a;
      _vertices[voff + i * L1 + j] = temp[0] * n_b * n_b
        + temp[1] * 2 * n_b * n_a + temp[2] * n_a * n_a;
    }
  }

  for (int i = 0; i <= tesselation_level; ++i)
    for (int j = 0; j <= tesselation_level; ++j) {
      int offset = eoff + (i * tesselation_level + j) * 6;
      _meshverts[offset + 0].offset = (i    ) * L1 + (j    ) + voff;
      _meshverts[offset + 1].offset = (i    ) * L1 + (j + 1) + voff;
      _meshverts[offset + 2].offset = (i + 1) * L1 + (j + 1) + voff;
      _meshverts[offset + 3].offset = (i + 1) * L1 + (j + 1) + voff;
      _meshverts[offset + 4].offset = (i + 1) * L1 + (j    ) + voff;
      _meshverts[offset + 5].offset = (i    ) * L1 + (j    ) + voff;
    }
}

void bsp::_create_patch_collides() {
  for (size_t i = 0; i < _faces.size(); ++i) {
    glm::vec3 points[max_patch_verts];
    int width = _faces[i].size[0], height = _faces[i].size[1]
      , c = width * height;
    bsp_vertex *dv_p = &_vertices[_faces[i].vertex];
    for (int j = 0; j < c; j++, dv_p++)
      points[j] = dv_p->position;
    // _patchcollides[i];
    assertf(width > 2 && height > 2, "bad patch size: (%d, %d)", width, height);
    assertf(!(!(width & 1) || !(height & 1)), "even sizes are invalid for"
        " quadratic meshes: (%d, %d)", width, height);
    assertf(width <= MAX_GRID_SIZE && height <= MAX_GRID_SIZE
        , "source is > MAX_GRID_SIZE");

    bsp_grid grid;
    grid.width = width;
    grid.height = height;
    grid.wrap_width = false;
    grid.wrap_height = false;
    for (int k = 0 ; k < width ; k++)
      for (int l = 0 ; l < height ; l++)
        grid.points[k][l] = points[l*width + k];

    // subdivide the grid
    grid.set_wrap_width();
    grid.subdivide_columns();
    grid.remove_degenerate_columns();

    grid.transpose();

    grid.set_wrap_width();
    grid.subdivide_columns();
    grid.remove_degenerate_columns();

    // we now have a grid of points exactly on the curve
    // the aproximate surface defined by these points will be
    // collided against
    _patchcollides[i].clear_bounds();
    for (int k = 0 ; k < grid.width ; k++ )
      for (int l = 0 ; l < grid.height ; l++ )
        _patchcollides[i].add_point_to_bounds(grid.points[k][l]);

    // c_totalPatchBlocks += ( grid.width - 1 ) * ( grid.height - 1 );

    // generate a bsp tree for the surface
    CM_PatchCollideFromGrid( &grid, pf );

    // expand by one unit for epsilon purposes
    pf->bounds[0][0] -= 1;
    pf->bounds[0][1] -= 1;
    pf->bounds[0][2] -= 1;

    pf->bounds[1][0] += 1;
    pf->bounds[1][1] += 1;
    pf->bounds[1][2] += 1;
  }
}

int bsp::_find_leaf(glm::vec3 position) {
  int index = 0;
  while (index >= 0) {
    bsp_node *node = &_nodes[index];
    bsp_plane *plane = &_planes[node->plane];
    if (plane->normal.x * position.x + plane->normal.y * position.y
        + plane->normal.z * position.z > plane->dist)
      index = node->front;
    else
      index = node->back;
  }
  return -(index + 1); // leaf index
}

int bsp::_cluster_visible(int vis_cluster, int test_cluster) {
  if (vis_cluster < 0 || _visdata.vecs.size() == 0)
    return 1;
  int i = vis_cluster * _visdata.sz_vecs + (test_cluster >> 3);
  unsigned char vis_set = _visdata.vecs[i];
  return vis_set & (1 << (test_cluster & 7));
}

void bsp::_set_visible_faces(const glm::vec3 &camera_pos, const frustum &f) {
  int leaf_index = _find_leaf(camera_pos);
  std::fill(_visible_faces.begin(), _visible_faces.end(), 0);
  for (bsp_leaf &l : _leaves)
    if (_cluster_visible(_leaves[leaf_index].cluster, l.cluster))
      if (f.box_in_frustum(l.mins, l.maxs))
        for (int j = 0; j < l.n_leaffaces; j++)
          _visible_faces[_leaffaces[l.leafface + j].face] = 1;
}

bsp::bsp(const char *filename, float world_scale, int tesselation_level)
  : _sp(shaders::map_vert, shaders::map_frag) {
  _load_file(filename, world_scale, tesselation_level);

  _vao.bind();

  _sp.use_this_prog();
  _vertex_pos_attr = _sp.bind_attrib("vertex_pos");
  _lightmap_coord_attr = _sp.bind_attrib("lightmap_coord");
  _mvp_mat_unif = _sp.bind_uniform("mvp");
  glUniform1i(_sp.bind_uniform("lightmap_sampler"), 1);

  _vbo.bind();
  _vbo.upload(sizeof(_vertices[0]) * _vertices.size(), &_vertices[0]);

  _ebo.bind();
  _ebo.upload(sizeof(_meshverts[0]) * _meshverts.size(), &_meshverts[0]);

  glEnableVertexAttribArray(_vertex_pos_attr);
  glEnableVertexAttribArray(_lightmap_coord_attr);

  glActiveTexture(GL_TEXTURE1);
  glGenTextures(_lightmaps.size(), &_lightmap_texture_ids[0]);
  for (size_t i = 0; i < _lightmaps.size(); ++i) {
    glBindTexture(GL_TEXTURE_2D, _lightmap_texture_ids[i]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 128, 128, 0, GL_RGB
        , GL_UNSIGNED_BYTE, _lightmaps[i].map);
    glGenerateMipmap(GL_TEXTURE_2D);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }
  glVertexAttribPointer(_vertex_pos_attr, 3, GL_FLOAT, GL_FALSE
      , sizeof(bsp_vertex), (void*)offsetof(bsp_vertex, position));
  glVertexAttribPointer(_lightmap_coord_attr, 2, GL_FLOAT, GL_FALSE
      , sizeof(bsp_vertex), (void*)offsetof(bsp_vertex, lightmap));

  _sp.dont_use_this_prog();
  _vao.unbind();
  _vbo.unbind();
  _ebo.unbind();
}

bsp::~bsp() {
  glDeleteTextures(_lightmaps.size(), _lightmap_texture_ids.data());
}

void bsp::draw(const glm::vec3 &position, const glm::mat4 &mvp, const frustum &f) {
  _set_visible_faces(position, f);

  _vao.bind();
  _sp.use_this_prog();

  glUniformMatrix4fv(_mvp_mat_unif, 1, GL_FALSE, glm::value_ptr(mvp));

  glActiveTexture(GL_TEXTURE1);
  for (size_t i = 0; i < _faces.size(); i++) {
    if (!_visible_faces[i])
      continue;
    if (_faces[i].type == (int)face::polygon
        || _faces[i].type == (int)face::mesh
        || _faces[i].type == (int)face::patch) {
      glBindTexture(GL_TEXTURE_2D, _lightmap_texture_ids[_faces[i].lm_index]);
      glDrawElements(GL_TRIANGLES, _faces[i].n_meshverts, GL_UNSIGNED_INT
          , (void*)(_faces[i].meshvert * sizeof(GLuint)));
    }
  }

  _vao.unbind();
  _sp.dont_use_this_prog();
}

void bsp::trace_sphere(trace_result *tr, const glm::vec3 &start
      , const glm::vec3 &end, float radius) {
  trace_description t { trace_type::sphere, radius };
  _trace(tr, t, start, end);
}

void bsp::_trace(trace_result *tr, const trace_description &td
    , const glm::vec3 &start, const glm::vec3 &end) {
  tr->clip_plane_normal = glm::vec3(0, 0, 0);
  tr->fraction = 1.0f;
  tr->end = end;
  tr->start_solid = tr->all_solid = false;

  _trace_node(tr, td, 0, 0.0f, 1.0f, start, end);

  if (tr->fraction == 1.0f)
    tr->end = end;
  else
    tr->end = start + tr->fraction * (end - start);
}

// Q3: CM_TraceThroughTree
void bsp::_trace_node(trace_result *tr, const trace_description &td
    , int node_index, float start_fraction, float end_fraction
    , const glm::vec3 &start, const glm::vec3 &end) {
  if (tr->fraction <= start_fraction)
    return;

  if (node_index < 0) {
    _trace_leaf(tr, td, &_leaves[-node_index - 1], start, end);
    return;
  }

  bsp_node *node = &_nodes[node_index];
  bsp_plane *plane = &_planes[node->plane];
  // TODO: plane->type not respected
  float start_dist = glm::dot(start, plane->normal) - plane->dist
    , end_dist = glm::dot(end, plane->normal) - plane->dist, offset = 0;

  switch (td.type) {
    case trace_type::ray:
      offset = 0;
      break;
    case trace_type::sphere:
      offset = td.radius;
      break;
    default:
      offset = 0;
      break;
  }

  if (start_dist >= offset + 1 && end_dist >= offset + 1) {
    _trace_node(tr, td, node->front, start_fraction, end_fraction, start, end);
    return;
  } else if (start_dist < -offset - 1 && end_dist < -offset - 1) {
    _trace_node(tr, td, node->back, start_fraction, end_fraction, start, end);
    return;
  }

  int side;
  float fraction1, fraction2;
  if (start_dist < end_dist) {
    const float inverse_dist = 1.0f / (start_dist - end_dist);
    side = 1; // back
    fraction1 = (start_dist - offset + surface_clip_eps) * inverse_dist;
    fraction2 = (start_dist + offset + surface_clip_eps) * inverse_dist;
  } else if (start_dist > end_dist) {
    const float inverse_dist = 1.0f / (start_dist - end_dist);
    side = 0; // front
    fraction1 = (start_dist + offset + surface_clip_eps) * inverse_dist;
    fraction2 = (start_dist - offset - surface_clip_eps) * inverse_dist;
  } else {
    side = 0;
    fraction1 = 1.0f;
    fraction2 = 0.0f;
  }

  // TODO
  // fraction1 = clamp(fraction1, 0.f, 1.f);
  // fraction2 = clamp(fraction2, 0.f, 1.f);
  if (fraction1 < 0.f)
    fraction1 = 0.f;
  else if (fraction1 > 1.f)
    fraction1 = 1.f;
  if (fraction2 < 0.f)
    fraction2 = 0.f;
  else if (fraction2 > 1.f)
    fraction2 = 1.f;

  float middle_fraction = start_fraction + (end_fraction - start_fraction)
    * fraction1;

  glm::vec3 middle = start + fraction1 * (end - start);

  if (side == 0)
    _trace_node(tr, td, node->front, start_fraction, middle_fraction, start
        , middle);
  else
    _trace_node(tr, td, node->back, start_fraction, middle_fraction, start
        , middle);

  middle_fraction = start_fraction + (end_fraction - start_fraction) * fraction2;
  middle = start + fraction2 * (end - start);

  if (side == 0)
    _trace_node(tr, td, node->back, middle_fraction, end_fraction, middle, end);
  else
    _trace_node(tr, td, node->front, middle_fraction, end_fraction, middle, end);
}

// Q3: CM_TraceThroughBrush
void bsp::_trace_brush(trace_result *tr, const trace_description &td
    , const bsp_brush *b, const glm::vec3 &input_start
    , const glm::vec3 &input_end) {
  bsp_plane *clip_plane = nullptr;
  float start_fraction = -1.0f, end_fraction = 1.0f;
  if (!b->n_brushsides)
    return;
  bool get_out = false, starts_out = false;
  for (int i = 0; i < b->n_brushsides; ++i) {
    bsp_brushside *brushside = &_brushsides[b->brushside + i];
    bsp_plane *plane = &_planes[brushside->plane];
    float start_dist, end_dist;

    switch (td.type) {
      case trace_type::ray:
        start_dist = glm::dot(input_start, plane->normal) - plane->dist;
        end_dist = glm::dot(input_end, plane->normal) - plane->dist;
        break;
      case trace_type::sphere:
        start_dist = glm::dot(input_start, plane->normal) - plane->dist - td.radius;
        end_dist = glm::dot(input_end, plane->normal) - plane->dist - td.radius;
        break;
      default:
        start_dist = glm::dot(input_start, plane->normal) - plane->dist;
        end_dist = glm::dot(input_end, plane->normal) - plane->dist;
        break;
    }

    if (end_dist > 0)
      get_out = true;

    if (start_dist > 0)
      starts_out = true;

    if (start_dist > 0 && (end_dist >= surface_clip_eps || end_dist >= start_dist))
      return;

    if (start_dist <= 0 && end_dist <= 0)
      continue;

    if (start_dist > end_dist) {
      float fraction = (start_dist - surface_clip_eps) / (start_dist - end_dist);
      if (fraction < 0)
        fraction = 0;
      if (fraction > start_fraction) {
        start_fraction = fraction;
        clip_plane = plane;
      }
    } else {
      float fraction = (start_dist + surface_clip_eps) / (start_dist - end_dist);
      if (fraction > 1)
        fraction = 1;
      if (fraction < end_fraction)
        end_fraction = fraction;
    }
  }

  if (!starts_out) {
    tr->start_solid = true;
    if (!get_out) {
      tr->all_solid = true;
      tr->fraction = 0;
    }
    return;
  }

  if (start_fraction < end_fraction)
    if (start_fraction > -1 && start_fraction < tr->fraction) {
      if (start_fraction < 0)
        start_fraction = 0;
      tr->fraction = start_fraction;
      tr->clip_plane_normal = clip_plane->normal;
    }
}

/// `const glm::vec3 &start' and `const glm::vec3 &end'
/// are all passed to only be passed to _trace_brush
void bsp::_trace_leaf(trace_result *tr, const trace_description &td
    , const bsp_leaf *l, const glm::vec3 &start
    , const glm::vec3 &end) {
  // trace line against all brushes in the leaf
  for (int k = 0; k < l->n_leafbrushes; ++k) {
    bsp_brush *b = &_brushes[_leafbrushes[l->leafbrush + k].brush];
    // if (b->checkcount == cm.checkcount)
    //   continue; // already checked this brush in another leaf
    // b->checkcount = cm.checkcount;
    _trace_brush(tr, td, b, start, end);
    if (!tr->fraction)
      return;
  }

  // trace line against all patches in the leaf
  for (int k = 0; k < l->n_leaffaces; ++k) {
    bsp_face *f = &_faces[_leaffaces[l->leafface + k].face];
    if (!f)
      continue;
    // if (patch->checkcount == cm.checkcount)
    //   continue; // already checked this patch in another leaf
    // patch->checkcount = cm.checkcount;
    float old_frac = tr->fraction;
    // CM_TraceThroughPatchCollide(tw, patch->pc);
    if (!tr->fraction)
      return;
  }
}

