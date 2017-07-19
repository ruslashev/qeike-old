#include "plane.hh"

idPlane::idPlane( void ) {
}

idPlane::idPlane( float a, float b, float c, float d ) {
  this->a = a;
  this->b = b;
  this->c = c;
  this->d = d;
}

idPlane::idPlane( const idVec3 &normal, const float dist ) {
  this->a = normal.x;
  this->b = normal.y;
  this->c = normal.z;
  this->d = -dist;
}

float idPlane::operator[]( int index ) const {
  return ( &a )[ index ];
}

float& idPlane::operator[]( int index ) {
  return ( &a )[ index ];
}

idPlane idPlane::operator-() const {
  return idPlane( -a, -b, -c, -d );
}

idPlane &idPlane::operator=( const idVec3 &v ) {
  a = v.x;
  b = v.y;
  c = v.z;
  d = 0;
  return *this;
}

idPlane idPlane::operator+( const idPlane &p ) const {
  return idPlane( a + p.a, b + p.b, c + p.c, d + p.d );
}

idPlane idPlane::operator-( const idPlane &p ) const {
  return idPlane( a - p.a, b - p.b, c - p.c, d - p.d );
}

idPlane &idPlane::operator*=( const idMat3 &m ) {
  Normal() *= m;
  return *this;
}

bool idPlane::Compare( const idPlane &p ) const {
  return ( a == p.a && b == p.b && c == p.c && d == p.d );
}

bool idPlane::Compare( const idPlane &p, const float epsilon ) const {
  if ( qkmath::fabs( a - p.a ) > epsilon ) {
    return false;
  }

  if ( qkmath::fabs( b - p.b ) > epsilon ) {
    return false;
  }

  if ( qkmath::fabs( c - p.c ) > epsilon ) {
    return false;
  }

  if ( qkmath::fabs( d - p.d ) > epsilon ) {
    return false;
  }

  return true;
}

bool idPlane::Compare( const idPlane &p, const float normalEps, const float distEps ) const {
  if ( qkmath::fabs( d - p.d ) > distEps ) {
    return false;
  }
  if ( !Normal().Compare( p.Normal(), normalEps ) ) {
    return false;
  }
  return true;
}

bool idPlane::operator==( const idPlane &p ) const {
  return Compare( p );
}

bool idPlane::operator!=( const idPlane &p ) const {
  return !Compare( p );
}

void idPlane::Zero( void ) {
  a = b = c = d = 0.0f;
}

void idPlane::SetNormal( const idVec3 &normal ) {
  a = normal.x;
  b = normal.y;
  c = normal.z;
}

const idVec3 &idPlane::Normal( void ) const {
  return *reinterpret_cast<const idVec3 *>(&a);
}

idVec3 &idPlane::Normal( void ) {
  return *reinterpret_cast<idVec3 *>(&a);
}

float idPlane::Normalize( bool fixDegenerate ) {
  float length = reinterpret_cast<idVec3 *>(&a)->Normalize();

  if ( fixDegenerate ) {
    FixDegenerateNormal();
  }
  return length;
}

bool idPlane::FixDegenerateNormal( void ) {
  return Normal().FixDegenerateNormal();
}

#if 0
bool idPlane::FixDegeneracies( float distEpsilon ) {
  bool fixedNormal = FixDegenerateNormal();
  // only fix dist if the normal was degenerate
  if ( fixedNormal ) {
    if ( qkmath::fabs( d - idMath::Rint( d ) ) < distEpsilon ) {
      d = idMath::Rint( d );
    }
  }
  return fixedNormal;
}
#endif

float idPlane::Dist( void ) const {
  return -d;
}

void idPlane::SetDist( const float dist ) {
  d = -dist;
}

bool idPlane::FromPoints( const idVec3 &p1, const idVec3 &p2, const idVec3 &p3, bool fixDegenerate ) {
  Normal() = (p1 - p2).Cross( p3 - p2 );
  if ( Normalize( fixDegenerate ) == 0.0f ) {
    return false;
  }
  d = -( Normal() * p2 );
  return true;
}

bool idPlane::FromVecs( const idVec3 &dir1, const idVec3 &dir2, const idVec3 &p, bool fixDegenerate ) {
  Normal() = dir1.Cross( dir2 );
  if ( Normalize( fixDegenerate ) == 0.0f ) {
    return false;
  }
  d = -( Normal() * p );
  return true;
}

void idPlane::FitThroughPoint( const idVec3 &p ) {
  d = -( Normal() * p );
}

idPlane idPlane::Translate( const idVec3 &translation ) const {
  return idPlane( a, b, c, d - translation * Normal() );
}

idPlane &idPlane::TranslateSelf( const idVec3 &translation ) {
  d -= translation * Normal();
  return *this;
}

idPlane idPlane::Rotate( const idVec3 &origin, const idMat3 &axis ) const {
  idPlane p;
  p.Normal() = Normal() * axis;
  p.d = d + origin * Normal() - origin * p.Normal();
  return p;
}

idPlane &idPlane::RotateSelf( const idVec3 &origin, const idMat3 &axis ) {
  d += origin * Normal();
  Normal() *= axis;
  d -= origin * Normal();
  return *this;
}

float idPlane::Distance( const idVec3 &v ) const {
  return a * v.x + b * v.y + c * v.z + d;
}

int idPlane::Side( const idVec3 &v, const float epsilon ) const {
  float dist = Distance( v );
  if ( dist > epsilon ) {
    return PLANESIDE_FRONT;
  }
  else if ( dist < -epsilon ) {
    return PLANESIDE_BACK;
  }
  else {
    return PLANESIDE_ON;
  }
}

bool idPlane::LineIntersection( const idVec3 &start, const idVec3 &end ) const {
  float d1, d2, fraction;

  d1 = Normal() * start + d;
  d2 = Normal() * end + d;
  if ( d1 == d2 ) {
    return false;
  }
  if ( d1 > 0.0f && d2 > 0.0f ) {
    return false;
  }
  if ( d1 < 0.0f && d2 < 0.0f ) {
    return false;
  }
  fraction = ( d1 / ( d1 - d2 ) );
  return ( fraction >= 0.0f && fraction <= 1.0f );
}

bool idPlane::RayIntersection( const idVec3 &start, const idVec3 &dir, float &scale ) const {
  float d1, d2;

  d1 = Normal() * start + d;
  d2 = Normal() * dir;
  if ( d2 == 0.0f ) {
    return false;
  }
  scale = -( d1 / d2 );
  return true;
}

int idPlane::GetDimension( void ) const {
  return 4;
}

const idVec4 &idPlane::ToVec4( void ) const {
  return *reinterpret_cast<const idVec4 *>(&a);
}

idVec4 &idPlane::ToVec4( void ) {
  return *reinterpret_cast<idVec4 *>(&a);
}

const float *idPlane::ToFloatPtr( void ) const {
  return reinterpret_cast<const float *>(&a);
}

float *idPlane::ToFloatPtr( void ) {
  return reinterpret_cast<float *>(&a);
}

idPlane plane_origin( 0.0f, 0.0f, 0.0f, 0.0f );

/*
================
idPlane::Type
================
*/
int idPlane::Type( void ) const {
	if ( Normal()[0] == 0.0f ) {
		if ( Normal()[1] == 0.0f ) {
			return Normal()[2] > 0.0f ? PLANETYPE_Z : PLANETYPE_NEGZ;
		}
		else if ( Normal()[2] == 0.0f ) {
			return Normal()[1] > 0.0f ? PLANETYPE_Y : PLANETYPE_NEGY;
		}
		else {
			return PLANETYPE_ZEROX;
		}
	}
	else if ( Normal()[1] == 0.0f ) {
		if ( Normal()[2] == 0.0f ) {
			return Normal()[0] > 0.0f ? PLANETYPE_X : PLANETYPE_NEGX;
		}
		else {
			return PLANETYPE_ZEROY;
		}
	}
	else if ( Normal()[2] == 0.0f ) {
		return PLANETYPE_ZEROZ;
	}
	else {
		return PLANETYPE_NONAXIAL;
	}
}

/*
================
idPlane::HeightFit
================
*/
bool idPlane::HeightFit( const idVec3 *points, const int numPoints ) {
	int i;
	float sumXX = 0.0f, sumXY = 0.0f, sumXZ = 0.0f;
	float sumYY = 0.0f, sumYZ = 0.0f;
	idVec3 sum, average, dir;

	if ( numPoints == 1 ) {
		a = 0.0f;
		b = 0.0f;
		c = 1.0f;
		d = -points[0].z;
		return true;
	}
	if ( numPoints == 2 ) {
		dir = points[1] - points[0];
		Normal() = dir.Cross( idVec3( 0, 0, 1 ) ).Cross( dir );
		Normalize();
		d = -( Normal() * points[0] );
		return true;
	}

	sum.Zero();
	for ( i = 0; i < numPoints; i++) {
		sum += points[i];
	}
	average = sum / numPoints;

	for ( i = 0; i < numPoints; i++ ) {
		dir = points[i] - average;
		sumXX += dir.x * dir.x;
		sumXY += dir.x * dir.y;
		sumXZ += dir.x * dir.z;
		sumYY += dir.y * dir.y;
		sumYZ += dir.y * dir.z;
	}

	idMat2 m( sumXX, sumXY, sumXY, sumYY );
	if ( !m.InverseSelf() ) {
		return false;
	}

	a = - sumXZ * m[0][0] - sumYZ * m[0][1];
	b = - sumXZ * m[1][0] - sumYZ * m[1][1];
	c = 1.0f;
	Normalize();
	d = -( a * average.x + b * average.y + c * average.z );
	return true;
}

/*
================
idPlane::PlaneIntersection
================
*/
bool idPlane::PlaneIntersection( const idPlane &plane, idVec3 &start, idVec3 &dir ) const {
	double n00, n01, n11, det, invDet, f0, f1;

	n00 = Normal().LengthSqr();
	n01 = Normal() * plane.Normal();
	n11 = plane.Normal().LengthSqr();
	det = n00 * n11 - n01 * n01;

	if ( qkmath::fabs(det) < 1e-6f ) {
		return false;
	}

	invDet = 1.0f / det;
	f0 = ( n01 * plane.d - n11 * d ) * invDet;
	f1 = ( n01 * d - n00 * plane.d ) * invDet;

	dir = Normal().Cross( plane.Normal() );
	start = f0 * Normal() + f1 * plane.Normal();
	return true;
}

#if 0
/*
=============
idPlane::ToString
=============
*/
const char *idPlane::ToString( int precision ) const {
	return idStr::FloatArrayToString( ToFloatPtr(), GetDimension(), precision );
}
#endif

