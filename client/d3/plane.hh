#pragma once

#include "vector.hh"
#include "matrix.hh"

class idVec3;
class idMat3;

#define	ON_EPSILON					0.1f
#define DEGENERATE_DIST_EPSILON		1e-4f

#define	SIDE_FRONT					0
#define	SIDE_BACK					1
#define	SIDE_ON						2
#define	SIDE_CROSS					3

// plane sides
#define PLANESIDE_FRONT				0
#define PLANESIDE_BACK				1
#define PLANESIDE_ON				2
#define PLANESIDE_CROSS				3

// plane types
#define PLANETYPE_X					0
#define PLANETYPE_Y					1
#define PLANETYPE_Z					2
#define PLANETYPE_NEGX				3
#define PLANETYPE_NEGY				4
#define PLANETYPE_NEGZ				5
#define PLANETYPE_TRUEAXIAL			6	// all types < 6 are true axial planes
#define PLANETYPE_ZEROX				6
#define PLANETYPE_ZEROY				7
#define PLANETYPE_ZEROZ				8
#define PLANETYPE_NONAXIAL			9

class idPlane {
public:
  idPlane( void );
  idPlane( float a, float b, float c, float d );
  idPlane( const idVec3 &normal, const float dist );

  float			operator[]( int index ) const;
  float &			operator[]( int index );
  idPlane			operator-() const;						// flips plane
  idPlane &		operator=( const idVec3 &v );			// sets normal and sets idPlane::d to zero
  idPlane			operator+( const idPlane &p ) const;	// add plane equations
  idPlane			operator-( const idPlane &p ) const;	// subtract plane equations
  idPlane &		operator*=( const idMat3 &m );			// Normal() *= m

  bool			Compare( const idPlane &p ) const;						// exact compare, no epsilon
  bool			Compare( const idPlane &p, const float epsilon ) const;	// compare with epsilon
  bool			Compare( const idPlane &p, const float normalEps, const float distEps ) const;	// compare with epsilon
  bool			operator==(	const idPlane &p ) const;					// exact compare, no epsilon
  bool			operator!=(	const idPlane &p ) const;					// exact compare, no epsilon

  void			Zero( void );							// zero plane
  void			SetNormal( const idVec3 &normal );		// sets the normal
  const idVec3 &	Normal( void ) const;					// reference to const normal
  idVec3 &		Normal( void );							// reference to normal
  float			Normalize( bool fixDegenerate = true );	// only normalizes the plane normal, does not adjust d
  bool			FixDegenerateNormal( void );			// fix degenerate normal
  bool			FixDegeneracies( float distEpsilon );	// fix degenerate normal and dist
  float			Dist( void ) const;						// returns: -d
  void			SetDist( const float dist );			// sets: d = -dist
  int				Type( void ) const;						// returns plane type

  bool			FromPoints( const idVec3 &p1, const idVec3 &p2, const idVec3 &p3, bool fixDegenerate = true );
  bool			FromVecs( const idVec3 &dir1, const idVec3 &dir2, const idVec3 &p, bool fixDegenerate = true );
  void			FitThroughPoint( const idVec3 &p );	// assumes normal is valid
  bool			HeightFit( const idVec3 *points, const int numPoints );
  idPlane			Translate( const idVec3 &translation ) const;
  idPlane &		TranslateSelf( const idVec3 &translation );
  idPlane			Rotate( const idVec3 &origin, const idMat3 &axis ) const;
  idPlane &		RotateSelf( const idVec3 &origin, const idMat3 &axis );

  float			Distance( const idVec3 &v ) const;
  int				Side( const idVec3 &v, const float epsilon = 0.0f ) const;

  bool			LineIntersection( const idVec3 &start, const idVec3 &end ) const;
  // intersection point is start + dir * scale
  bool			RayIntersection( const idVec3 &start, const idVec3 &dir, float &scale ) const;
  bool			PlaneIntersection( const idPlane &plane, idVec3 &start, idVec3 &dir ) const;

  int				GetDimension( void ) const;

  const idVec4 &	ToVec4( void ) const;
  idVec4 &		ToVec4( void );
  const float *	ToFloatPtr( void ) const;
  float *			ToFloatPtr( void );
  const char *	ToString( int precision = 2 ) const;

private:
  float			a;
  float			b;
  float			c;
  float			d;
};

extern idPlane plane_origin;
#define plane_zero plane_origin

