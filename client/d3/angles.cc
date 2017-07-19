#include "angles.hh"
#include "quat.hh"
#include "rotation.hh"

idAngles::idAngles( void ) {
}

idAngles::idAngles( float pitch, float yaw, float roll ) {
  this->pitch = pitch;
  this->yaw	= yaw;
  this->roll	= roll;
}

idAngles::idAngles( const idVec3 &v ) {
  this->pitch = v[0];
  this->yaw	= v[1];
  this->roll	= v[2];
}

void idAngles::Set( float pitch, float yaw, float roll ) {
  this->pitch = pitch;
  this->yaw	= yaw;
  this->roll	= roll;
}

idAngles &idAngles::Zero( void ) {
  pitch = yaw = roll = 0.0f;
  return *this;
}

float idAngles::operator[]( int index ) const {
  assert( ( index >= 0 ) && ( index < 3 ) );
  return ( &pitch )[ index ];
}

float &idAngles::operator[]( int index ) {
  assert( ( index >= 0 ) && ( index < 3 ) );
  return ( &pitch )[ index ];
}

idAngles idAngles::operator-() const {
  return idAngles( -pitch, -yaw, -roll );
}

idAngles &idAngles::operator=( const idAngles &a ) {
  pitch	= a.pitch;
  yaw		= a.yaw;
  roll	= a.roll;
  return *this;
}

idAngles idAngles::operator+( const idAngles &a ) const {
  return idAngles( pitch + a.pitch, yaw + a.yaw, roll + a.roll );
}

idAngles& idAngles::operator+=( const idAngles &a ) {
  pitch	+= a.pitch;
  yaw		+= a.yaw;
  roll	+= a.roll;

  return *this;
}

idAngles idAngles::operator-( const idAngles &a ) const {
  return idAngles( pitch - a.pitch, yaw - a.yaw, roll - a.roll );
}

idAngles& idAngles::operator-=( const idAngles &a ) {
  pitch	-= a.pitch;
  yaw		-= a.yaw;
  roll	-= a.roll;

  return *this;
}

idAngles idAngles::operator*( const float a ) const {
  return idAngles( pitch * a, yaw * a, roll * a );
}

idAngles& idAngles::operator*=( float a ) {
  pitch	*= a;
  yaw		*= a;
  roll	*= a;
  return *this;
}

idAngles idAngles::operator/( const float a ) const {
  float inva = 1.0f / a;
  return idAngles( pitch * inva, yaw * inva, roll * inva );
}

idAngles& idAngles::operator/=( float a ) {
  float inva = 1.0f / a;
  pitch	*= inva;
  yaw		*= inva;
  roll	*= inva;
  return *this;
}

idAngles operator*( const float a, const idAngles &b ) {
  return idAngles( a * b.pitch, a * b.yaw, a * b.roll );
}

bool idAngles::Compare( const idAngles &a ) const {
  return ( ( a.pitch == pitch ) && ( a.yaw == yaw ) && ( a.roll == roll ) );
}

bool idAngles::Compare( const idAngles &a, const float epsilon ) const {
  if ( qkmath::fabs( pitch - a.pitch ) > epsilon ) {
    return false;
  }

  if ( qkmath::fabs( yaw - a.yaw ) > epsilon ) {
    return false;
  }

  if ( qkmath::fabs( roll - a.roll ) > epsilon ) {
    return false;
  }

  return true;
}

bool idAngles::operator==( const idAngles &a ) const {
  return Compare( a );
}

bool idAngles::operator!=( const idAngles &a ) const {
  return !Compare( a );
}

void idAngles::Clamp( const idAngles &min, const idAngles &max ) {
  if ( pitch < min.pitch ) {
    pitch = min.pitch;
  } else if ( pitch > max.pitch ) {
    pitch = max.pitch;
  }
  if ( yaw < min.yaw ) {
    yaw = min.yaw;
  } else if ( yaw > max.yaw ) {
    yaw = max.yaw;
  }
  if ( roll < min.roll ) {
    roll = min.roll;
  } else if ( roll > max.roll ) {
    roll = max.roll;
  }
}

int idAngles::GetDimension( void ) const {
  return 3;
}

const float *idAngles::ToFloatPtr( void ) const {
  return &pitch;
}

float *idAngles::ToFloatPtr( void ) {
  return &pitch;
}

idAngles ang_zero( 0.0f, 0.0f, 0.0f );

/*
=================
idAngles::Normalize360

returns angles normalized to the range [0 <= angle < 360]
=================
*/
idAngles& idAngles::Normalize360( void ) {
	int i;

	for ( i = 0; i < 3; i++ ) {
		if ( ( (*this)[i] >= 360.0f ) || ( (*this)[i] < 0.0f ) ) {
			(*this)[i] -= floor( (*this)[i] / 360.0f ) * 360.0f;

			if ( (*this)[i] >= 360.0f ) {
				(*this)[i] -= 360.0f;
			}
			if ( (*this)[i] < 0.0f ) {
				(*this)[i] += 360.0f;
			}
		}
	}

	return *this;
}

/*
=================
idAngles::Normalize180

returns angles normalized to the range [-180 < angle <= 180]
=================
*/
idAngles& idAngles::Normalize180( void ) {
	Normalize360();

	if ( pitch > 180.0f ) {
		pitch -= 360.0f;
	}

	if ( yaw > 180.0f ) {
		yaw -= 360.0f;
	}

	if ( roll > 180.0f ) {
		roll -= 360.0f;
	}
	return *this;
}

/*
=================
idAngles::ToVectors
=================
*/
void idAngles::ToVectors( idVec3 *forward, idVec3 *right, idVec3 *up ) const {
	float sr, sp, sy, cr, cp, cy;

	qkmath::sincos( DEG2RAD( yaw ), sy, cy );
	qkmath::sincos( DEG2RAD( pitch ), sp, cp );
	qkmath::sincos( DEG2RAD( roll ), sr, cr );

	if ( forward ) {
		forward->Set( cp * cy, cp * sy, -sp );
	}

	if ( right ) {
		right->Set( -sr * sp * cy + cr * sy, -sr * sp * sy + -cr * cy, -sr * cp );
	}

	if ( up ) {
		up->Set( cr * sp * cy + -sr * -sy, cr * sp * sy + -sr * cy, cr * cp );
	}
}

/*
=================
idAngles::ToForward
=================
*/
idVec3 idAngles::ToForward( void ) const {
	float sp, sy, cp, cy;

	qkmath::sincos( DEG2RAD( yaw ), sy, cy );
	qkmath::sincos( DEG2RAD( pitch ), sp, cp );

	return idVec3( cp * cy, cp * sy, -sp );
}

/*
=================
idAngles::ToQuat
=================
*/
idQuat idAngles::ToQuat( void ) const {
	float sx, cx, sy, cy, sz, cz;
	float sxcy, cxcy, sxsy, cxsy;

	qkmath::sincos( DEG2RAD( yaw ) * 0.5f, sz, cz );
	qkmath::sincos( DEG2RAD( pitch ) * 0.5f, sy, cy );
	qkmath::sincos( DEG2RAD( roll ) * 0.5f, sx, cx );

	sxcy = sx * cy;
	cxcy = cx * cy;
	sxsy = sx * sy;
	cxsy = cx * sy;

	return idQuat( cxsy*sz - sxcy*cz, -cxsy*cz - sxcy*sz, sxsy*cz - cxcy*sz, cxcy*cz + sxsy*sz );
}

/*
=================
idAngles::ToRotation
=================
*/
idRotation idAngles::ToRotation( void ) const {
	idVec3 vec;
	float angle, w;
	float sx, cx, sy, cy, sz, cz;
	float sxcy, cxcy, sxsy, cxsy;

	if ( pitch == 0.0f ) {
		if ( yaw == 0.0f ) {
			return idRotation( vec3_origin, idVec3( -1.0f, 0.0f, 0.0f ), roll );
		}
		if ( roll == 0.0f ) {
			return idRotation( vec3_origin, idVec3( 0.0f, 0.0f, -1.0f ), yaw );
		}
	} else if ( yaw == 0.0f && roll == 0.0f ) {
		return idRotation( vec3_origin, idVec3( 0.0f, -1.0f, 0.0f ), pitch );
	}

	qkmath::sincos( DEG2RAD( yaw ) * 0.5f, sz, cz );
	qkmath::sincos( DEG2RAD( pitch ) * 0.5f, sy, cy );
	qkmath::sincos( DEG2RAD( roll ) * 0.5f, sx, cx );

	sxcy = sx * cy;
	cxcy = cx * cy;
	sxsy = sx * sy;
	cxsy = cx * sy;

	vec.x =  cxsy * sz - sxcy * cz;
	vec.y = -cxsy * cz - sxcy * sz;
	vec.z =  sxsy * cz - cxcy * sz;
	w =		 cxcy * cz + sxsy * sz;
	angle = qkmath::acos( w );
	if ( angle == 0.0f ) {
		vec.Set( 0.0f, 0.0f, 1.0f );
	} else {
		//vec *= (1.0f / sin( angle ));
		vec.Normalize();
		vec.FixDegenerateNormal();
		angle *= 2.0f * qkmath::C_RAD2DEG;
	}
	return idRotation( vec3_origin, vec, angle );
}

/*
=================
idAngles::ToMat3
=================
*/
idMat3 idAngles::ToMat3( void ) const {
	idMat3 mat;
	float sr, sp, sy, cr, cp, cy;

	qkmath::sincos( DEG2RAD( yaw ), sy, cy );
	qkmath::sincos( DEG2RAD( pitch ), sp, cp );
	qkmath::sincos( DEG2RAD( roll ), sr, cr );

	mat[ 0 ].Set( cp * cy, cp * sy, -sp );
	mat[ 1 ].Set( sr * sp * cy + cr * -sy, sr * sp * sy + cr * cy, sr * cp );
	mat[ 2 ].Set( cr * sp * cy + -sr * -sy, cr * sp * sy + -sr * cy, cr * cp );

	return mat;
}

/*
=================
idAngles::ToMat4
=================
*/
idMat4 idAngles::ToMat4( void ) const {
	return ToMat3().ToMat4();
}

/*
=================
idAngles::ToAngularVelocity
=================
*/
idVec3 idAngles::ToAngularVelocity( void ) const {
	idRotation rotation = idAngles::ToRotation();
	return rotation.GetVec() * DEG2RAD( rotation.GetAngle() );
}

