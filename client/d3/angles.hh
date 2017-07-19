#pragma once

#include "vector.hh"

// angle indexes
#define	PITCH				0		// up / down
#define	YAW					1		// left / right
#define	ROLL				2		// fall over

class idVec3;
class idQuat;
class idRotation;
class idMat3;
class idMat4;

class idAngles {
public:
	float			pitch;
	float			yaw;
	float			roll;

					idAngles( void );
					idAngles( float pitch, float yaw, float roll );
					explicit idAngles( const idVec3 &v );

	void			Set( float pitch, float yaw, float roll );
	idAngles &		Zero( void );

	float			operator[]( int index ) const;
	float &			operator[]( int index );
	idAngles		operator-() const;			// negate angles, in general not the inverse rotation
	idAngles &		operator=( const idAngles &a );
	idAngles		operator+( const idAngles &a ) const;
	idAngles &		operator+=( const idAngles &a );
	idAngles		operator-( const idAngles &a ) const;
	idAngles &		operator-=( const idAngles &a );
	idAngles		operator*( const float a ) const;
	idAngles &		operator*=( const float a );
	idAngles		operator/( const float a ) const;
	idAngles &		operator/=( const float a );

	friend idAngles	operator*( const float a, const idAngles &b );

	bool			Compare( const idAngles &a ) const;							// exact compare, no epsilon
	bool			Compare( const idAngles &a, const float epsilon ) const;	// compare with epsilon
	bool			operator==(	const idAngles &a ) const;						// exact compare, no epsilon
	bool			operator!=(	const idAngles &a ) const;						// exact compare, no epsilon

	idAngles &		Normalize360( void );	// normalizes 'this'
	idAngles &		Normalize180( void );	// normalizes 'this'

	void			Clamp( const idAngles &min, const idAngles &max );

	int				GetDimension( void ) const;

	void			ToVectors( idVec3 *forward, idVec3 *right = NULL, idVec3 *up = NULL ) const;
	idVec3			ToForward( void ) const;
	idQuat			ToQuat( void ) const;
	idRotation		ToRotation( void ) const;
	idMat3			ToMat3( void ) const;
	idMat4			ToMat4( void ) const;
	idVec3			ToAngularVelocity( void ) const;
	const float *	ToFloatPtr( void ) const;
	float *			ToFloatPtr( void );
	const char *	ToString( int precision = 2 ) const;
};

extern idAngles ang_zero;

