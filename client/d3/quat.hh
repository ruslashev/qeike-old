#pragma once

#include "vector.hh"
#include "rotation.hh"
#include "matrix.hh"

class idRotation;
class idMat4;

class idQuat {
public:
  float			x;
  float			y;
  float			z;
  float			w;

  idQuat( void );
  idQuat( float x, float y, float z, float w );

  void			Set( float x, float y, float z, float w );

  float			operator[]( int index ) const;
  float &			operator[]( int index );
  idQuat			operator-() const;
  idQuat &		operator=( const idQuat &a );
  idQuat			operator+( const idQuat &a ) const;
  idQuat &		operator+=( const idQuat &a );
  idQuat			operator-( const idQuat &a ) const;
  idQuat &		operator-=( const idQuat &a );
  idQuat			operator*( const idQuat &a ) const;
  idVec3			operator*( const idVec3 &a ) const;
  idQuat			operator*( float a ) const;
  idQuat &		operator*=( const idQuat &a );
  idQuat &		operator*=( float a );

  friend idQuat	operator*( const float a, const idQuat &b );
  friend idVec3	operator*( const idVec3 &a, const idQuat &b );

  bool			Compare( const idQuat &a ) const;						// exact compare, no epsilon
  bool			Compare( const idQuat &a, const float epsilon ) const;	// compare with epsilon
  bool			operator==(	const idQuat &a ) const;					// exact compare, no epsilon
  bool			operator!=(	const idQuat &a ) const;					// exact compare, no epsilon

  idQuat			Inverse( void ) const;
  float			Length( void ) const;
  idQuat &		Normalize( void );

  float			CalcW( void ) const;
  int				GetDimension( void ) const;

  idAngles		ToAngles( void ) const;
  idRotation		ToRotation( void ) const;
  idMat3			ToMat3( void ) const;
  idMat4			ToMat4( void ) const;
  idVec3			ToAngularVelocity( void ) const;
  const float *	ToFloatPtr( void ) const;
  float *			ToFloatPtr( void );
  const char *	ToString( int precision = 2 ) const;

  idQuat &		Slerp( const idQuat &from, const idQuat &to, float t );
};

