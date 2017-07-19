#include "pluecker.hh"
#include "plane.hh"

idPluecker::idPluecker( void ) {
}

idPluecker::idPluecker( const float *a ) {
  memcpy( p, a, 6 * sizeof( float ) );
}

idPluecker::idPluecker( const idVec3 &start, const idVec3 &end ) {
  FromLine( start, end );
}

idPluecker::idPluecker( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 ) {
  p[0] = a1;
  p[1] = a2;
  p[2] = a3;
  p[3] = a4;
  p[4] = a5;
  p[5] = a6;
}

idPluecker idPluecker::operator-() const {
  return idPluecker( -p[0], -p[1], -p[2], -p[3], -p[4], -p[5] );
}

float idPluecker::operator[]( const int index ) const {
  return p[index];
}

float &idPluecker::operator[]( const int index ) {
  return p[index];
}

idPluecker idPluecker::operator*( const float a ) const {
  return idPluecker( p[0]*a, p[1]*a, p[2]*a, p[3]*a, p[4]*a, p[5]*a );
}

float idPluecker::operator*( const idPluecker &a ) const {
  return p[0] * a.p[4] + p[1] * a.p[5] + p[2] * a.p[3] + p[4] * a.p[0] + p[5] * a.p[1] + p[3] * a.p[2];
}

idPluecker idPluecker::operator/( const float a ) const {
  float inva;

  assert( a != 0.0f );
  inva = 1.0f / a;
  return idPluecker( p[0]*inva, p[1]*inva, p[2]*inva, p[3]*inva, p[4]*inva, p[5]*inva );
}

idPluecker idPluecker::operator+( const idPluecker &a ) const {
  return idPluecker( p[0] + a[0], p[1] + a[1], p[2] + a[2], p[3] + a[3], p[4] + a[4], p[5] + a[5] );
}

idPluecker idPluecker::operator-( const idPluecker &a ) const {
  return idPluecker( p[0] - a[0], p[1] - a[1], p[2] - a[2], p[3] - a[3], p[4] - a[4], p[5] - a[5] );
}

idPluecker &idPluecker::operator*=( const float a ) {
  p[0] *= a;
  p[1] *= a;
  p[2] *= a;
  p[3] *= a;
  p[4] *= a;
  p[5] *= a;
  return *this;
}

idPluecker &idPluecker::operator/=( const float a ) {
  float inva;

  assert( a != 0.0f );
  inva = 1.0f / a;
  p[0] *= inva;
  p[1] *= inva;
  p[2] *= inva;
  p[3] *= inva;
  p[4] *= inva;
  p[5] *= inva;
  return *this;
}

idPluecker &idPluecker::operator+=( const idPluecker &a ) {
  p[0] += a[0];
  p[1] += a[1];
  p[2] += a[2];
  p[3] += a[3];
  p[4] += a[4];
  p[5] += a[5];
  return *this;
}

idPluecker &idPluecker::operator-=( const idPluecker &a ) {
  p[0] -= a[0];
  p[1] -= a[1];
  p[2] -= a[2];
  p[3] -= a[3];
  p[4] -= a[4];
  p[5] -= a[5];
  return *this;
}

bool idPluecker::Compare( const idPluecker &a ) const {
  return ( ( p[0] == a[0] ) && ( p[1] == a[1] ) && ( p[2] == a[2] ) &&
      ( p[3] == a[3] ) && ( p[4] == a[4] ) && ( p[5] == a[5] ) );
}

bool idPluecker::Compare( const idPluecker &a, const float epsilon ) const {
  if ( qkmath::fabs( p[0] - a[0] ) > epsilon ) {
    return false;
  }

  if ( qkmath::fabs( p[1] - a[1] ) > epsilon ) {
    return false;
  }

  if ( qkmath::fabs( p[2] - a[2] ) > epsilon ) {
    return false;
  }

  if ( qkmath::fabs( p[3] - a[3] ) > epsilon ) {
    return false;
  }

  if ( qkmath::fabs( p[4] - a[4] ) > epsilon ) {
    return false;
  }

  if ( qkmath::fabs( p[5] - a[5] ) > epsilon ) {
    return false;
  }

  return true;
}

bool idPluecker::operator==( const idPluecker &a ) const {
  return Compare( a );
}

bool idPluecker::operator!=( const idPluecker &a ) const {
  return !Compare( a );
}

void idPluecker::Set( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 ) {
  p[0] = a1;
  p[1] = a2;
  p[2] = a3;
  p[3] = a4;
  p[4] = a5;
  p[5] = a6;
}

void idPluecker::Zero( void ) {
  p[0] = p[1] = p[2] = p[3] = p[4] = p[5] = 0.0f;
}

void idPluecker::FromLine( const idVec3 &start, const idVec3 &end ) {
  p[0] = start[0] * end[1] - end[0] * start[1];
  p[1] = start[0] * end[2] - end[0] * start[2];
  p[2] = start[0] - end[0];
  p[3] = start[1] * end[2] - end[1] * start[2];
  p[4] = start[2] - end[2];
  p[5] = end[1] - start[1];
}

void idPluecker::FromRay( const idVec3 &start, const idVec3 &dir ) {
  p[0] = start[0] * dir[1] - dir[0] * start[1];
  p[1] = start[0] * dir[2] - dir[0] * start[2];
  p[2] = -dir[0];
  p[3] = start[1] * dir[2] - dir[1] * start[2];
  p[4] = -dir[2];
  p[5] = dir[1];
}

bool idPluecker::ToLine( idVec3 &start, idVec3 &end ) const {
  idVec3 dir1, dir2;
  float d;

  dir1[0] = p[3];
  dir1[1] = -p[1];
  dir1[2] = p[0];

  dir2[0] = -p[2];
  dir2[1] = p[5];
  dir2[2] = -p[4];

  d = dir2 * dir2;
  if ( d == 0.0f ) {
    return false; // pluecker coordinate does not represent a line
  }

  start = dir2.Cross(dir1) * (1.0f / d);
  end = start + dir2;
  return true;
}

bool idPluecker::ToRay( idVec3 &start, idVec3 &dir ) const {
  idVec3 dir1;
  float d;

  dir1[0] = p[3];
  dir1[1] = -p[1];
  dir1[2] = p[0];

  dir[0] = -p[2];
  dir[1] = p[5];
  dir[2] = -p[4];

  d = dir * dir;
  if ( d == 0.0f ) {
    return false; // pluecker coordinate does not represent a line
  }

  start = dir.Cross(dir1) * (1.0f / d);
  return true;
}

void idPluecker::ToDir( idVec3 &dir ) const {
  dir[0] = -p[2];
  dir[1] = p[5];
  dir[2] = -p[4];
}

float idPluecker::PermutedInnerProduct( const idPluecker &a ) const {
  return p[0] * a.p[4] + p[1] * a.p[5] + p[2] * a.p[3] + p[4] * a.p[0] + p[5] * a.p[1] + p[3] * a.p[2];
}

float idPluecker::Length( void ) const {
  return ( float )qkmath::sqrt( p[5] * p[5] + p[4] * p[4] + p[2] * p[2] );
}

float idPluecker::LengthSqr( void ) const {
  return ( p[5] * p[5] + p[4] * p[4] + p[2] * p[2] );
}

float idPluecker::NormalizeSelf( void ) {
  float l, d;

  l = LengthSqr();
  if ( l == 0.0f ) {
    return l; // pluecker coordinate does not represent a line
  }
  d = qkmath::inv_sqrt( l );
  p[0] *= d;
  p[1] *= d;
  p[2] *= d;
  p[3] *= d;
  p[4] *= d;
  p[5] *= d;
  return d * l;
}

idPluecker idPluecker::Normalize( void ) const {
  float d;

  d = LengthSqr();
  if ( d == 0.0f ) {
    return *this; // pluecker coordinate does not represent a line
  }
  d = qkmath::inv_sqrt( d );
  return idPluecker( p[0]*d, p[1]*d, p[2]*d, p[3]*d, p[4]*d, p[5]*d );
}

int idPluecker::GetDimension( void ) const {
  return 6;
}

const float *idPluecker::ToFloatPtr( void ) const {
  return p;
}

float *idPluecker::ToFloatPtr( void ) {
  return p;
}

idPluecker pluecker_origin( 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f );

/*
================
idPluecker::FromPlanes

  pluecker coordinate for the intersection of two planes
================
*/
#if 0
bool idPluecker::FromPlanes( const idPlane &p1, const idPlane &p2 ) {

	p[0] = -( p1[2] * -p2[3] - p2[2] * -p1[3] );
	p[1] = -( p2[1] * -p1[3] - p1[1] * -p2[3] );
	p[2] = p1[1] * p2[2] - p2[1] * p1[2];

	p[3] = -( p1[0] * -p2[3] - p2[0] * -p1[3] );
	p[4] = p1[0] * p2[1] - p2[0] * p1[1];
	p[5] = p1[0] * p2[2] - p2[0] * p1[2];

	return ( p[2] != 0.0f || p[5] != 0.0f || p[4] != 0.0f );
}
#endif

/*
================
idPluecker::Distance3DSqr

  calculates square of shortest distance between the two
  3D lines represented by their pluecker coordinates
================
*/
float idPluecker::Distance3DSqr( const idPluecker &a ) const {
	float d, s;
	idVec3 dir;

	dir[0] = -a.p[5] *  p[4] -  a.p[4] * -p[5];
	dir[1] =  a.p[4] *  p[2] -  a.p[2] *  p[4];
	dir[2] =  a.p[2] * -p[5] - -a.p[5] *  p[2];
	if ( dir[0] == 0.0f && dir[1] == 0.0f && dir[2] == 0.0f ) {
		return -1.0f;	// FIXME: implement for parallel lines
	}
	d = a.p[4] * ( p[2]*dir[1] - -p[5]*dir[0]) +
		a.p[5] * ( p[2]*dir[2] -  p[4]*dir[0]) +
		a.p[2] * (-p[5]*dir[2] -  p[4]*dir[1]);
	s = PermutedInnerProduct( a ) / d;
	return ( dir * dir ) * ( s * s );
}

