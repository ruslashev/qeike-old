#pragma once

#include "vector.hh"

class idPluecker {
public:
					idPluecker( void );
					explicit idPluecker( const float *a );
					explicit idPluecker( const idVec3 &start, const idVec3 &end );
					explicit idPluecker( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 );

	float			operator[]( const int index ) const;
	float &			operator[]( const int index );
	idPluecker		operator-() const;											// flips the direction
	idPluecker		operator*( const float a ) const;
	idPluecker		operator/( const float a ) const;
	float			operator*( const idPluecker &a ) const;						// permuted inner product
	idPluecker		operator-( const idPluecker &a ) const;
	idPluecker		operator+( const idPluecker &a ) const;
	idPluecker &	operator*=( const float a );
	idPluecker &	operator/=( const float a );
	idPluecker &	operator+=( const idPluecker &a );
	idPluecker &	operator-=( const idPluecker &a );

	bool			Compare( const idPluecker &a ) const;						// exact compare, no epsilon
	bool			Compare( const idPluecker &a, const float epsilon ) const;	// compare with epsilon
	bool			operator==(	const idPluecker &a ) const;					// exact compare, no epsilon
	bool			operator!=(	const idPluecker &a ) const;					// exact compare, no epsilon

	void			Set( const float a1, const float a2, const float a3, const float a4, const float a5, const float a6 );
	void			Zero( void );

	void			FromLine( const idVec3 &start, const idVec3 &end );			// pluecker from line
	void			FromRay( const idVec3 &start, const idVec3 &dir );			// pluecker from ray
	// bool			FromPlanes( const idPlane &p1, const idPlane &p2 );			// pluecker from intersection of planes
	bool			ToLine( idVec3 &start, idVec3 &end ) const;					// pluecker to line
	bool			ToRay( idVec3 &start, idVec3 &dir ) const;					// pluecker to ray
	void			ToDir( idVec3 &dir ) const;									// pluecker to direction
	float			PermutedInnerProduct( const idPluecker &a ) const;			// pluecker permuted inner product
	float			Distance3DSqr( const idPluecker &a ) const;					// pluecker line distance

	float			Length( void ) const;										// pluecker length
	float			LengthSqr( void ) const;									// pluecker squared length
	idPluecker		Normalize( void ) const;									// pluecker normalize
	float			NormalizeSelf( void );										// pluecker normalize

	int				GetDimension( void ) const;

	const float *	ToFloatPtr( void ) const;
	float *			ToFloatPtr( void );
	const char *	ToString( int precision = 2 ) const;

private:
	float			p[6];
};

extern idPluecker pluecker_origin;
#define pluecker_zero pluecker_origin

