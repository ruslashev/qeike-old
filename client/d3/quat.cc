#include "quat.hh"

idQuat::idQuat( void ) {
}

idQuat::idQuat( float x, float y, float z, float w ) {
  this->x = x;
  this->y = y;
  this->z = z;
  this->w = w;
}

float idQuat::operator[]( int index ) const {
  assert( ( index >= 0 ) && ( index < 4 ) );
  return ( &x )[ index ];
}

float& idQuat::operator[]( int index ) {
  assert( ( index >= 0 ) && ( index < 4 ) );
  return ( &x )[ index ];
}

idQuat idQuat::operator-() const {
  return idQuat( -x, -y, -z, -w );
}

idQuat &idQuat::operator=( const idQuat &a ) {
  x = a.x;
  y = a.y;
  z = a.z;
  w = a.w;

  return *this;
}

idQuat idQuat::operator+( const idQuat &a ) const {
  return idQuat( x + a.x, y + a.y, z + a.z, w + a.w );
}

idQuat& idQuat::operator+=( const idQuat &a ) {
  x += a.x;
  y += a.y;
  z += a.z;
  w += a.w;

  return *this;
}

idQuat idQuat::operator-( const idQuat &a ) const {
  return idQuat( x - a.x, y - a.y, z - a.z, w - a.w );
}

idQuat& idQuat::operator-=( const idQuat &a ) {
  x -= a.x;
  y -= a.y;
  z -= a.z;
  w -= a.w;

  return *this;
}

idQuat idQuat::operator*( const idQuat &a ) const {
  return idQuat(	w*a.x + x*a.w + y*a.z - z*a.y,
      w*a.y + y*a.w + z*a.x - x*a.z,
      w*a.z + z*a.w + x*a.y - y*a.x,
      w*a.w - x*a.x - y*a.y - z*a.z );
}

idVec3 idQuat::operator*( const idVec3 &a ) const {
#if 0
  // it's faster to do the conversion to a 3x3 matrix and multiply the vector by this 3x3 matrix
  return ( ToMat3() * a );
#else
  // result = this->Inverse() * idQuat( a.x, a.y, a.z, 0.0f ) * (*this)
  float xxzz = x*x - z*z;
  float wwyy = w*w - y*y;

  float xw2 = x*w*2.0f;
  float xy2 = x*y*2.0f;
  float xz2 = x*z*2.0f;
  float yw2 = y*w*2.0f;
  float yz2 = y*z*2.0f;
  float zw2 = z*w*2.0f;

  return idVec3(
      (xxzz + wwyy)*a.x		+ (xy2 + zw2)*a.y		+ (xz2 - yw2)*a.z,
      (xy2 - zw2)*a.x			+ (y*y+w*w-x*x-z*z)*a.y	+ (yz2 + xw2)*a.z,
      (xz2 + yw2)*a.x			+ (yz2 - xw2)*a.y		+ (wwyy - xxzz)*a.z
      );
#endif
}

idQuat idQuat::operator*( float a ) const {
  return idQuat( x * a, y * a, z * a, w * a );
}

idQuat operator*( const float a, const idQuat &b ) {
  return b * a;
}

idVec3 operator*( const idVec3 &a, const idQuat &b ) {
  return b * a;
}

idQuat& idQuat::operator*=( const idQuat &a ) {
  *this = *this * a;

  return *this;
}

idQuat& idQuat::operator*=( float a ) {
  x *= a;
  y *= a;
  z *= a;
  w *= a;

  return *this;
}

bool idQuat::Compare( const idQuat &a ) const {
  return ( ( x == a.x ) && ( y == a.y ) && ( z == a.z ) && ( w == a.w ) );
}

bool idQuat::Compare( const idQuat &a, const float epsilon ) const {
  if ( qkmath::fabs( x - a.x ) > epsilon ) {
    return false;
  }
  if ( qkmath::fabs( y - a.y ) > epsilon ) {
    return false;
  }
  if ( qkmath::fabs( z - a.z ) > epsilon ) {
    return false;
  }
  if ( qkmath::fabs( w - a.w ) > epsilon ) {
    return false;
  }
  return true;
}

bool idQuat::operator==( const idQuat &a ) const {
  return Compare( a );
}

bool idQuat::operator!=( const idQuat &a ) const {
  return !Compare( a );
}

void idQuat::Set( float x, float y, float z, float w ) {
  this->x = x;
  this->y = y;
  this->z = z;
  this->w = w;
}

idQuat idQuat::Inverse( void ) const {
  return idQuat( -x, -y, -z, w );
}

float idQuat::Length( void ) const {
  float len;

  len = x * x + y * y + z * z + w * w;
  return qkmath::sqrt( len );
}

idQuat& idQuat::Normalize( void ) {
  float len;
  float ilength;

  len = this->Length();
  if ( len ) {
    ilength = 1 / len;
    x *= ilength;
    y *= ilength;
    z *= ilength;
    w *= ilength;
  }
  return *this;
}

float idQuat::CalcW( void ) const {
  // take the absolute value because floating point rounding may cause the dot of x,y,z to be larger than 1
  return sqrt( fabs( 1.0f - ( x * x + y * y + z * z ) ) );
}

int idQuat::GetDimension( void ) const {
  return 4;
}

const float *idQuat::ToFloatPtr( void ) const {
  return &x;
}

float *idQuat::ToFloatPtr( void ) {
  return &x;
}

