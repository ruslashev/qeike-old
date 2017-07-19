#include "rotation.hh"

idRotation::idRotation( void ) {
}

idRotation::idRotation( const idVec3 &rotationOrigin, const idVec3 &rotationVec, const float rotationAngle ) {
  origin = rotationOrigin;
  vec = rotationVec;
  angle = rotationAngle;
  axisValid = false;
}

void idRotation::Set( const idVec3 &rotationOrigin, const idVec3 &rotationVec, const float rotationAngle ) {
  origin = rotationOrigin;
  vec = rotationVec;
  angle = rotationAngle;
  axisValid = false;
}

void idRotation::SetOrigin( const idVec3 &rotationOrigin ) {
  origin = rotationOrigin;
}

void idRotation::SetVec( const idVec3 &rotationVec ) {
  vec = rotationVec;
  axisValid = false;
}

void idRotation::SetVec( float x, float y, float z ) {
  vec[0] = x;
  vec[1] = y;
  vec[2] = z;
  axisValid = false;
}

void idRotation::SetAngle( const float rotationAngle ) {
  angle = rotationAngle;
  axisValid = false;
}

void idRotation::Scale( const float s ) {
  angle *= s;
  axisValid = false;
}

void idRotation::ReCalculateMatrix( void ) {
  axisValid = false;
  ToMat3();
}

const idVec3 &idRotation::GetOrigin( void ) const {
  return origin;
}

const idVec3 &idRotation::GetVec( void ) const  {
  return vec;
}

float idRotation::GetAngle( void ) const  {
  return angle;
}

idRotation idRotation::operator-() const {
  return idRotation( origin, vec, -angle );
}

idRotation idRotation::operator*( const float s ) const {
  return idRotation( origin, vec, angle * s );
}

idRotation idRotation::operator/( const float s ) const {
  assert( s != 0.0f );
  return idRotation( origin, vec, angle / s );
}

idRotation &idRotation::operator*=( const float s ) {
  angle *= s;
  axisValid = false;
  return *this;
}

idRotation &idRotation::operator/=( const float s ) {
  assert( s != 0.0f );
  angle /= s;
  axisValid = false;
  return *this;
}

idVec3 idRotation::operator*( const idVec3 &v ) const {
  if ( !axisValid ) {
    ToMat3();
  }
  return ((v - origin) * axis + origin);
}

idRotation operator*( const float s, const idRotation &r ) {
  return r * s;
}

idVec3 operator*( const idVec3 &v, const idRotation &r ) {
  return r * v;
}

idVec3 &operator*=( idVec3 &v, const idRotation &r ) {
  v = r * v;
  return v;
}

void idRotation::RotatePoint( idVec3 &point ) const {
  if ( !axisValid ) {
    ToMat3();
  }
  point = ((point - origin) * axis + origin);
}



#if 0
/*
   ============
   idRotation::ToAngles
   ============
   */
idAngles idRotation::ToAngles( void ) const {
  return ToMat3().ToAngles();
}

/*
   ============
   idRotation::ToQuat
   ============
   */
idQuat idRotation::ToQuat( void ) const {
  float a, s, c;

  a = angle * ( idMath::M_DEG2RAD * 0.5f );
  idMath::SinCos( a, s, c );
  return idQuat( vec.x * s, vec.y * s, vec.z * s, c );
}
#endif

/*
   ============
   idRotation::toMat3
   ============
   */
const idMat3 &idRotation::ToMat3( void ) const {
  float wx, wy, wz;
  float xx, yy, yz;
  float xy, xz, zz;
  float x2, y2, z2;
  float a, c, s, x, y, z;

  if ( axisValid ) {
    return axis;
  }

  a = angle * ( qkmath::C_DEG2RAD * 0.5f );
  s = sinf(a);
  c = cosf(a);

  x = vec[0] * s;
  y = vec[1] * s;
  z = vec[2] * s;

  x2 = x + x;
  y2 = y + y;
  z2 = z + z;

  xx = x * x2;
  xy = x * y2;
  xz = x * z2;

  yy = y * y2;
  yz = y * z2;
  zz = z * z2;

  wx = c * x2;
  wy = c * y2;
  wz = c * z2;

  axis[ 0 ][ 0 ] = 1.0f - ( yy + zz );
  axis[ 0 ][ 1 ] = xy - wz;
  axis[ 0 ][ 2 ] = xz + wy;

  axis[ 1 ][ 0 ] = xy + wz;
  axis[ 1 ][ 1 ] = 1.0f - ( xx + zz );
  axis[ 1 ][ 2 ] = yz - wx;

  axis[ 2 ][ 0 ] = xz - wy;
  axis[ 2 ][ 1 ] = yz + wx;
  axis[ 2 ][ 2 ] = 1.0f - ( xx + yy );

  axisValid = true;

  return axis;
}

/*
   ============
   idRotation::ToMat4
   ============
   */
idMat4 idRotation::ToMat4( void ) const {
  return ToMat3().ToMat4();
}

/*
   ============
   idRotation::ToAngularVelocity
   ============
   */
idVec3 idRotation::ToAngularVelocity( void ) const {
  return vec * DEG2RAD( angle );
}

/*
   ============
   idRotation::Normalize180
   ============
   */
void idRotation::Normalize180( void ) {
  angle -= floor( angle / 360.0f ) * 360.0f;
  if ( angle > 180.0f ) {
    angle -= 360.0f;
  }
  else if ( angle < -180.0f ) {
    angle += 360.0f;
  }
}

/*
   ============
   idRotation::Normalize360
   ============
   */
void idRotation::Normalize360( void ) {
  angle -= floor( angle / 360.0f ) * 360.0f;
  if ( angle > 360.0f ) {
    angle -= 360.0f;
  }
  else if ( angle < 0.0f ) {
    angle += 360.0f;
  }
}

