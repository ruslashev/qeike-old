#pragma once

#include "vector.hh"
#include "matrix.hh"
#include "quat.hh"

class idRotation {
  friend class idAngles;
  friend class idQuat;
  friend class idMat3;

public:
  idRotation( void );
  idRotation( const idVec3 &rotationOrigin, const idVec3 &rotationVec, const float rotationAngle );

  void				Set( const idVec3 &rotationOrigin, const idVec3 &rotationVec, const float rotationAngle );
  void				SetOrigin( const idVec3 &rotationOrigin );
  void				SetVec( const idVec3 &rotationVec );					// has to be normalized
  void				SetVec( const float x, const float y, const float z );	// has to be normalized
  void				SetAngle( const float rotationAngle );
  void				Scale( const float s );
  void				ReCalculateMatrix( void );
  const idVec3 &		GetOrigin( void ) const;
  const idVec3 &		GetVec( void ) const;
  float				GetAngle( void ) const;

  idRotation			operator-() const;										// flips rotation
  idRotation			operator*( const float s ) const;						// scale rotation
  idRotation			operator/( const float s ) const;						// scale rotation
  idRotation &		operator*=( const float s );							// scale rotation
  idRotation &		operator/=( const float s );							// scale rotation
  idVec3				operator*( const idVec3 &v ) const;						// rotate vector

  friend idRotation	operator*( const float s, const idRotation &r );		// scale rotation
  friend idVec3		operator*( const idVec3 &v, const idRotation &r );		// rotate vector
  friend idVec3 &		operator*=( idVec3 &v, const idRotation &r );			// rotate vector

  idAngles			ToAngles( void ) const;
  idQuat				ToQuat( void ) const;
  const idMat3 &		ToMat3( void ) const;
  idMat4				ToMat4( void ) const;
  idVec3				ToAngularVelocity( void ) const;

  void				RotatePoint( idVec3 &point ) const;

  void				Normalize180( void );
  void				Normalize360( void );

private:
  idVec3				origin;			// origin of rotation
  idVec3				vec;			// normalized vector to rotate around
  float				angle;			// angle of rotation in degrees
  mutable idMat3		axis;			// rotation axis
  mutable bool		axisValid;		// true if rotation axis is valid
};

