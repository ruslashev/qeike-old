#include "bounds.hh"

idBounds::idBounds( void ) {
}

idBounds::idBounds( const idVec3 &mins, const idVec3 &maxs ) {
  b[0] = mins;
  b[1] = maxs;
}

idBounds::idBounds( const idVec3 &point ) {
  b[0] = point;
  b[1] = point;
}

const idVec3 &idBounds::operator[]( const int index ) const {
  return b[index];
}

idVec3 &idBounds::operator[]( const int index ) {
  return b[index];
}

idBounds idBounds::operator+( const idVec3 &t ) const {
  return idBounds( b[0] + t, b[1] + t );
}

idBounds &idBounds::operator+=( const idVec3 &t ) {
  b[0] += t;
  b[1] += t;
  return *this;
}

idBounds idBounds::operator*( const idMat3 &r ) const {
  idBounds bounds;
  bounds.FromTransformedBounds( *this, vec3_origin, r );
  return bounds;
}

idBounds &idBounds::operator*=( const idMat3 &r ) {
  this->FromTransformedBounds( *this, vec3_origin, r );
  return *this;
}

idBounds idBounds::operator+( const idBounds &a ) const {
  idBounds newBounds;
  newBounds = *this;
  newBounds.AddBounds( a );
  return newBounds;
}

idBounds &idBounds::operator+=( const idBounds &a ) {
  idBounds::AddBounds( a );
  return *this;
}

idBounds idBounds::operator-( const idBounds &a ) const {
  assert( b[1][0] - b[0][0] > a.b[1][0] - a.b[0][0] &&
      b[1][1] - b[0][1] > a.b[1][1] - a.b[0][1] &&
      b[1][2] - b[0][2] > a.b[1][2] - a.b[0][2] );
  return idBounds( idVec3( b[0][0] + a.b[1][0], b[0][1] + a.b[1][1], b[0][2] + a.b[1][2] ),
      idVec3( b[1][0] + a.b[0][0], b[1][1] + a.b[0][1], b[1][2] + a.b[0][2] ) );
}

idBounds &idBounds::operator-=( const idBounds &a ) {
  assert( b[1][0] - b[0][0] > a.b[1][0] - a.b[0][0] &&
      b[1][1] - b[0][1] > a.b[1][1] - a.b[0][1] &&
      b[1][2] - b[0][2] > a.b[1][2] - a.b[0][2] );
  b[0] += a.b[1];
  b[1] += a.b[0];
  return *this;
}

bool idBounds::Compare( const idBounds &a ) const {
  return ( b[0].Compare( a.b[0] ) && b[1].Compare( a.b[1] ) );
}

bool idBounds::Compare( const idBounds &a, const float epsilon ) const {
  return ( b[0].Compare( a.b[0], epsilon ) && b[1].Compare( a.b[1], epsilon ) );
}

bool idBounds::operator==( const idBounds &a ) const {
  return Compare( a );
}

bool idBounds::operator!=( const idBounds &a ) const {
  return !Compare( a );
}

void idBounds::Clear( void ) {
  b[0][0] = b[0][1] = b[0][2] = qkmath::C_INFINITY;
  b[1][0] = b[1][1] = b[1][2] = -qkmath::C_INFINITY;
}

void idBounds::Zero( void ) {
  b[0][0] = b[0][1] = b[0][2] =
    b[1][0] = b[1][1] = b[1][2] = 0;
}

idVec3 idBounds::GetCenter( void ) const {
  return idVec3( ( b[1][0] + b[0][0] ) * 0.5f, ( b[1][1] + b[0][1] ) * 0.5f, ( b[1][2] + b[0][2] ) * 0.5f );
}

float idBounds::GetVolume( void ) const {
  if ( b[0][0] >= b[1][0] || b[0][1] >= b[1][1] || b[0][2] >= b[1][2] ) {
    return 0.0f;
  }
  return ( ( b[1][0] - b[0][0] ) * ( b[1][1] - b[0][1] ) * ( b[1][2] - b[0][2] ) );
}

bool idBounds::IsCleared( void ) const {
  return b[0][0] > b[1][0];
}

bool idBounds::AddPoint( const idVec3 &v ) {
  bool expanded = false;
  if ( v[0] < b[0][0]) {
    b[0][0] = v[0];
    expanded = true;
  }
  if ( v[0] > b[1][0]) {
    b[1][0] = v[0];
    expanded = true;
  }
  if ( v[1] < b[0][1] ) {
    b[0][1] = v[1];
    expanded = true;
  }
  if ( v[1] > b[1][1]) {
    b[1][1] = v[1];
    expanded = true;
  }
  if ( v[2] < b[0][2] ) {
    b[0][2] = v[2];
    expanded = true;
  }
  if ( v[2] > b[1][2]) {
    b[1][2] = v[2];
    expanded = true;
  }
  return expanded;
}

bool idBounds::AddBounds( const idBounds &a ) {
  bool expanded = false;
  if ( a.b[0][0] < b[0][0] ) {
    b[0][0] = a.b[0][0];
    expanded = true;
  }
  if ( a.b[0][1] < b[0][1] ) {
    b[0][1] = a.b[0][1];
    expanded = true;
  }
  if ( a.b[0][2] < b[0][2] ) {
    b[0][2] = a.b[0][2];
    expanded = true;
  }
  if ( a.b[1][0] > b[1][0] ) {
    b[1][0] = a.b[1][0];
    expanded = true;
  }
  if ( a.b[1][1] > b[1][1] ) {
    b[1][1] = a.b[1][1];
    expanded = true;
  }
  if ( a.b[1][2] > b[1][2] ) {
    b[1][2] = a.b[1][2];
    expanded = true;
  }
  return expanded;
}

idBounds idBounds::Intersect( const idBounds &a ) const {
  idBounds n;
  n.b[0][0] = ( a.b[0][0] > b[0][0] ) ? a.b[0][0] : b[0][0];
  n.b[0][1] = ( a.b[0][1] > b[0][1] ) ? a.b[0][1] : b[0][1];
  n.b[0][2] = ( a.b[0][2] > b[0][2] ) ? a.b[0][2] : b[0][2];
  n.b[1][0] = ( a.b[1][0] < b[1][0] ) ? a.b[1][0] : b[1][0];
  n.b[1][1] = ( a.b[1][1] < b[1][1] ) ? a.b[1][1] : b[1][1];
  n.b[1][2] = ( a.b[1][2] < b[1][2] ) ? a.b[1][2] : b[1][2];
  return n;
}

idBounds &idBounds::IntersectSelf( const idBounds &a ) {
  if ( a.b[0][0] > b[0][0] ) {
    b[0][0] = a.b[0][0];
  }
  if ( a.b[0][1] > b[0][1] ) {
    b[0][1] = a.b[0][1];
  }
  if ( a.b[0][2] > b[0][2] ) {
    b[0][2] = a.b[0][2];
  }
  if ( a.b[1][0] < b[1][0] ) {
    b[1][0] = a.b[1][0];
  }
  if ( a.b[1][1] < b[1][1] ) {
    b[1][1] = a.b[1][1];
  }
  if ( a.b[1][2] < b[1][2] ) {
    b[1][2] = a.b[1][2];
  }
  return *this;
}

idBounds idBounds::Expand( const float d ) const {
  return idBounds( idVec3( b[0][0] - d, b[0][1] - d, b[0][2] - d ),
      idVec3( b[1][0] + d, b[1][1] + d, b[1][2] + d ) );
}

idBounds &idBounds::ExpandSelf( const float d ) {
  b[0][0] -= d;
  b[0][1] -= d;
  b[0][2] -= d;
  b[1][0] += d;
  b[1][1] += d;
  b[1][2] += d;
  return *this;
}

idBounds idBounds::Translate( const idVec3 &translation ) const {
  return idBounds( b[0] + translation, b[1] + translation );
}

idBounds &idBounds::TranslateSelf( const idVec3 &translation ) {
  b[0] += translation;
  b[1] += translation;
  return *this;
}

idBounds idBounds::Rotate( const idMat3 &rotation ) const {
  idBounds bounds;
  bounds.FromTransformedBounds( *this, vec3_origin, rotation );
  return bounds;
}

idBounds &idBounds::RotateSelf( const idMat3 &rotation ) {
  FromTransformedBounds( *this, vec3_origin, rotation );
  return *this;
}

bool idBounds::ContainsPoint( const idVec3 &p ) const {
  if ( p[0] < b[0][0] || p[1] < b[0][1] || p[2] < b[0][2]
      || p[0] > b[1][0] || p[1] > b[1][1] || p[2] > b[1][2] ) {
    return false;
  }
  return true;
}

bool idBounds::IntersectsBounds( const idBounds &a ) const {
  if ( a.b[1][0] < b[0][0] || a.b[1][1] < b[0][1] || a.b[1][2] < b[0][2]
      || a.b[0][0] > b[1][0] || a.b[0][1] > b[1][1] || a.b[0][2] > b[1][2] ) {
    return false;
  }
  return true;
}

#if 0
idSphere idBounds::ToSphere( void ) const {
  idSphere sphere;
  sphere.SetOrigin( ( b[0] + b[1] ) * 0.5f );
  sphere.SetRadius( ( b[1] - sphere.GetOrigin() ).Length() );
  return sphere;
}
#endif

void idBounds::AxisProjection( const idVec3 &dir, float &min, float &max ) const {
  float d1, d2;
  idVec3 center, extents;

  center = ( b[0] + b[1] ) * 0.5f;
  extents = b[1] - center;

  d1 = dir * center;
  d2 = qkmath::fabs( extents[0] * dir[0] ) + qkmath::fabs( extents[1] * dir[1] ) +
    qkmath::fabs( extents[2] * dir[2] );

  min = d1 - d2;
  max = d1 + d2;
}

void idBounds::AxisProjection( const idVec3 &origin, const idMat3 &axis, const idVec3 &dir, float &min, float &max ) const {
  float d1, d2;
  idVec3 center, extents;

  center = ( b[0] + b[1] ) * 0.5f;
  extents = b[1] - center;
  center = origin + center * axis;

  d1 = dir * center;
  d2 = qkmath::fabs( extents[0] * ( dir * axis[0] ) ) +
    qkmath::fabs( extents[1] * ( dir * axis[1] ) ) +
    qkmath::fabs( extents[2] * ( dir * axis[2] ) );

  min = d1 - d2;
  max = d1 + d2;
}

idBounds bounds_zero( vec3_zero, vec3_zero );

/*
============
idBounds::GetRadius
============
*/
float idBounds::GetRadius( void ) const {
	int		i;
	float	total, b0, b1;

	total = 0.0f;
	for ( i = 0; i < 3; i++ ) {
		b0 = (float)qkmath::fabs( b[0][i] );
		b1 = (float)qkmath::fabs( b[1][i] );
		if ( b0 > b1 ) {
			total += b0 * b0;
		} else {
			total += b1 * b1;
		}
	}
	return qkmath::sqrt( total );
}

/*
============
idBounds::GetRadius
============
*/
float idBounds::GetRadius( const idVec3 &center ) const {
	int		i;
	float	total, b0, b1;

	total = 0.0f;
	for ( i = 0; i < 3; i++ ) {
		b0 = (float)qkmath::fabs( center[i] - b[0][i] );
		b1 = (float)qkmath::fabs( b[1][i] - center[i] );
		if ( b0 > b1 ) {
			total += b0 * b0;
		} else {
			total += b1 * b1;
		}
	}
	return qkmath::sqrt( total );
}

/*
================
idBounds::PlaneDistance
================
*/
float idBounds::PlaneDistance( const idPlane &plane ) const {
	idVec3 center;
	float d1, d2;

	center = ( b[0] + b[1] ) * 0.5f;

	d1 = plane.Distance( center );
	d2 = qkmath::fabs( ( b[1][0] - center[0] ) * plane.Normal()[0] ) +
			qkmath::fabs( ( b[1][1] - center[1] ) * plane.Normal()[1] ) +
				qkmath::fabs( ( b[1][2] - center[2] ) * plane.Normal()[2] );

	if ( d1 - d2 > 0.0f ) {
		return d1 - d2;
	}
	if ( d1 + d2 < 0.0f ) {
		return d1 + d2;
	}
	return 0.0f;
}

/*
================
idBounds::PlaneSide
================
*/
int idBounds::PlaneSide( const idPlane &plane, const float epsilon ) const {
	idVec3 center;
	float d1, d2;

	center = ( b[0] + b[1] ) * 0.5f;

	d1 = plane.Distance( center );
	d2 = qkmath::fabs( ( b[1][0] - center[0] ) * plane.Normal()[0] ) +
			qkmath::fabs( ( b[1][1] - center[1] ) * plane.Normal()[1] ) +
				qkmath::fabs( ( b[1][2] - center[2] ) * plane.Normal()[2] );

	if ( d1 - d2 > epsilon ) {
		return PLANESIDE_FRONT;
	}
	if ( d1 + d2 < -epsilon ) {
		return PLANESIDE_BACK;
	}
	return PLANESIDE_CROSS;
}

/*
============
idBounds::LineIntersection

  Returns true if the line intersects the bounds between the start and end point.
============
*/
bool idBounds::LineIntersection( const idVec3 &start, const idVec3 &end ) const {
	float ld[3];
	idVec3 center = ( b[0] + b[1] ) * 0.5f;
	idVec3 extents = b[1] - center;
	idVec3 lineDir = 0.5f * ( end - start );
	idVec3 lineCenter = start + lineDir;
	idVec3 dir = lineCenter - center;

	ld[0] = qkmath::fabs( lineDir[0] );
	if ( qkmath::fabs( dir[0] ) > extents[0] + ld[0] ) {
		return false;
	}

	ld[1] = qkmath::fabs( lineDir[1] );
	if ( qkmath::fabs( dir[1] ) > extents[1] + ld[1] ) {
		return false;
	}

	ld[2] = qkmath::fabs( lineDir[2] );
	if ( qkmath::fabs( dir[2] ) > extents[2] + ld[2] ) {
		return false;
	}

	idVec3 cross = lineDir.Cross( dir );

	if ( qkmath::fabs( cross[0] ) > extents[1] * ld[2] + extents[2] * ld[1] ) {
		return false;
	}

	if ( qkmath::fabs( cross[1] ) > extents[0] * ld[2] + extents[2] * ld[0] ) {
		return false;
	}

	if ( qkmath::fabs( cross[2] ) > extents[0] * ld[1] + extents[1] * ld[0] ) {
		return false;
	}

	return true;
}

/*
============
idBounds::RayIntersection

  Returns true if the ray intersects the bounds.
  The ray can intersect the bounds in both directions from the start point.
  If start is inside the bounds it is considered an intersection with scale = 0
============
*/
bool idBounds::RayIntersection( const idVec3 &start, const idVec3 &dir, float &scale ) const {
	int i, ax0, ax1, ax2, side, inside;
	float f;
	idVec3 hit;

	ax0 = -1;
	inside = 0;
	for ( i = 0; i < 3; i++ ) {
		if ( start[i] < b[0][i] ) {
			side = 0;
		}
		else if ( start[i] > b[1][i] ) {
			side = 1;
		}
		else {
			inside++;
			continue;
		}
		if ( dir[i] == 0.0f ) {
			continue;
		}
		f = ( start[i] - b[side][i] );
		if ( ax0 < 0 || qkmath::fabs( f ) > qkmath::fabs( scale * dir[i] ) ) {
			scale = - ( f / dir[i] );
			ax0 = i;
		}
	}

	if ( ax0 < 0 ) {
		scale = 0.0f;
		// return true if the start point is inside the bounds
		return ( inside == 3 );
	}

	ax1 = (ax0+1)%3;
	ax2 = (ax0+2)%3;
	hit[ax1] = start[ax1] + scale * dir[ax1];
	hit[ax2] = start[ax2] + scale * dir[ax2];

	return ( hit[ax1] >= b[0][ax1] && hit[ax1] <= b[1][ax1] &&
				hit[ax2] >= b[0][ax2] && hit[ax2] <= b[1][ax2] );
}

/*
============
idBounds::FromTransformedBounds
============
*/
void idBounds::FromTransformedBounds( const idBounds &bounds, const idVec3 &origin, const idMat3 &axis ) {
	int i;
	idVec3 center, extents, rotatedExtents;

	center = (bounds[0] + bounds[1]) * 0.5f;
	extents = bounds[1] - center;

	for ( i = 0; i < 3; i++ ) {
		rotatedExtents[i] = qkmath::fabs( extents[0] * axis[0][i] ) +
							qkmath::fabs( extents[1] * axis[1][i] ) +
							qkmath::fabs( extents[2] * axis[2][i] );
	}

	center = origin + center * axis;
	b[0] = center - rotatedExtents;
	b[1] = center + rotatedExtents;
}

#if 0
/*
============
idBounds::FromPoints

  Most tight bounds for a point set.
============
*/
void idBounds::FromPoints( const idVec3 *points, const int numPoints ) {
	SIMDProcessor->MinMax( b[0], b[1], points, numPoints );
}
#endif

/*
============
idBounds::FromPointTranslation

  Most tight bounds for the translational movement of the given point.
============
*/
void idBounds::FromPointTranslation( const idVec3 &point, const idVec3 &translation ) {
	int i;

	for ( i = 0; i < 3; i++ ) {
		if ( translation[i] < 0.0f ) {
			b[0][i] = point[i] + translation[i];
			b[1][i] = point[i];
		}
		else {
			b[0][i] = point[i];
			b[1][i] = point[i] + translation[i];
		}
	}
}

/*
============
idBounds::FromBoundsTranslation

  Most tight bounds for the translational movement of the given bounds.
============
*/
void idBounds::FromBoundsTranslation( const idBounds &bounds, const idVec3 &origin, const idMat3 &axis, const idVec3 &translation ) {
	int i;

	if ( axis.IsRotated() ) {
		FromTransformedBounds( bounds, origin, axis );
	}
	else {
		b[0] = bounds[0] + origin;
		b[1] = bounds[1] + origin;
	}
	for ( i = 0; i < 3; i++ ) {
		if ( translation[i] < 0.0f ) {
			b[0][i] += translation[i];
		}
		else {
			b[1][i] += translation[i];
		}
	}
}

/*
================
BoundsForPointRotation

  only for rotations < 180 degrees
================
*/
idBounds BoundsForPointRotation( const idVec3 &start, const idRotation &rotation ) {
	int i;
	float radiusSqr;
	idVec3 v1, v2;
	idVec3 origin, axis, end;
	idBounds bounds;

	end = start * rotation;
	axis = rotation.GetVec();
	origin = rotation.GetOrigin() + axis * ( axis * ( start - rotation.GetOrigin() ) );
	radiusSqr = ( start - origin ).LengthSqr();
	v1 = ( start - origin ).Cross( axis );
	v2 = ( end - origin ).Cross( axis );

	for ( i = 0; i < 3; i++ ) {
		// if the derivative changes sign along this axis during the rotation from start to end
		if ( ( v1[i] > 0.0f && v2[i] < 0.0f ) || ( v1[i] < 0.0f && v2[i] > 0.0f ) ) {
			if ( ( 0.5f * (start[i] + end[i]) - origin[i] ) > 0.0f ) {
				bounds[0][i] = qkmath::min( start[i], end[i] );
				bounds[1][i] = origin[i] + qkmath::sqrt( radiusSqr * ( 1.0f - axis[i] * axis[i] ) );
			}
			else {
				bounds[0][i] = origin[i] - qkmath::sqrt( radiusSqr * ( 1.0f - axis[i] * axis[i] ) );
				bounds[1][i] = qkmath::max( start[i], end[i] );
			}
		}
		else if ( start[i] > end[i] ) {
			bounds[0][i] = end[i];
			bounds[1][i] = start[i];
		}
		else {
			bounds[0][i] = start[i];
			bounds[1][i] = end[i];
		}
	}

	return bounds;
}

/*
============
idBounds::FromPointRotation

  Most tight bounds for the rotational movement of the given point.
============
*/
void idBounds::FromPointRotation( const idVec3 &point, const idRotation &rotation ) {
	float radius;

	if ( qkmath::fabs( rotation.GetAngle() ) < 180.0f ) {
		(*this) = BoundsForPointRotation( point, rotation );
	}
	else {

		radius = ( point - rotation.GetOrigin() ).Length();

		// FIXME: these bounds are usually way larger
		b[0].Set( -radius, -radius, -radius );
		b[1].Set( radius, radius, radius );
	}
}

/*
============
idBounds::FromBoundsRotation

  Most tight bounds for the rotational movement of the given bounds.
============
*/
void idBounds::FromBoundsRotation( const idBounds &bounds, const idVec3 &origin, const idMat3 &axis, const idRotation &rotation ) {
	int i;
	float radius;
	idVec3 point;
	idBounds rBounds;

	if ( qkmath::fabs( rotation.GetAngle() ) < 180.0f ) {

		(*this) = BoundsForPointRotation( bounds[0] * axis + origin, rotation );
		for ( i = 1; i < 8; i++ ) {
			point[0] = bounds[(i^(i>>1))&1][0];
			point[1] = bounds[(i>>1)&1][1];
			point[2] = bounds[(i>>2)&1][2];
			(*this) += BoundsForPointRotation( point * axis + origin, rotation );
		}
	}
	else {

		point = (bounds[1] - bounds[0]) * 0.5f;
		radius = (bounds[1] - point).Length() + (point - rotation.GetOrigin()).Length();

		// FIXME: these bounds are usually way larger
		b[0].Set( -radius, -radius, -radius );
		b[1].Set( radius, radius, radius );
	}
}

/*
============
idBounds::ToPoints
============
*/
void idBounds::ToPoints( idVec3 points[8] ) const {
	for ( int i = 0; i < 8; i++ ) {
		points[i][0] = b[(i^(i>>1))&1][0];
		points[i][1] = b[(i>>1)&1][1];
		points[i][2] = b[(i>>2)&1][2];
	}
}

