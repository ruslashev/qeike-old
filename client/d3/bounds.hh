#pragma once

#include "vector.hh"
#include "plane.hh"
#include "rotation.hh"

class idBounds {
public:
  idBounds( void );
  explicit idBounds( const idVec3 &mins, const idVec3 &maxs );
  explicit idBounds( const idVec3 &point );

  const idVec3 &	operator[]( const int index ) const;
  idVec3 &		operator[]( const int index );
  idBounds		operator+( const idVec3 &t ) const;				// returns translated bounds
  idBounds &		operator+=( const idVec3 &t );					// translate the bounds
  idBounds		operator*( const idMat3 &r ) const;				// returns rotated bounds
  idBounds &		operator*=( const idMat3 &r );					// rotate the bounds
  idBounds		operator+( const idBounds &a ) const;
  idBounds &		operator+=( const idBounds &a );
  idBounds		operator-( const idBounds &a ) const;
  idBounds &		operator-=( const idBounds &a );

  bool			Compare( const idBounds &a ) const;							// exact compare, no epsilon
  bool			Compare( const idBounds &a, const float epsilon ) const;	// compare with epsilon
  bool			operator==(	const idBounds &a ) const;						// exact compare, no epsilon
  bool			operator!=(	const idBounds &a ) const;						// exact compare, no epsilon

  void			Clear( void );									// inside out bounds
  void			Zero( void );									// single point at origin

  idVec3			GetCenter( void ) const;						// returns center of bounds
  float			GetRadius( void ) const;						// returns the radius relative to the bounds origin
  float			GetRadius( const idVec3 &center ) const;		// returns the radius relative to the given center
  float			GetVolume( void ) const;						// returns the volume of the bounds
  bool			IsCleared( void ) const;						// returns true if bounds are inside out

  bool			AddPoint( const idVec3 &v );					// add the point, returns true if the bounds expanded
  bool			AddBounds( const idBounds &a );					// add the bounds, returns true if the bounds expanded
  idBounds		Intersect( const idBounds &a ) const;			// return intersection of this bounds with the given bounds
  idBounds &		IntersectSelf( const idBounds &a );				// intersect this bounds with the given bounds
  idBounds		Expand( const float d ) const;					// return bounds expanded in all directions with the given value
  idBounds &		ExpandSelf( const float d );					// expand bounds in all directions with the given value
  idBounds		Translate( const idVec3 &translation ) const;	// return translated bounds
  idBounds &		TranslateSelf( const idVec3 &translation );		// translate this bounds
  idBounds		Rotate( const idMat3 &rotation ) const;			// return rotated bounds
  idBounds &		RotateSelf( const idMat3 &rotation );			// rotate this bounds

  float			PlaneDistance( const idPlane &plane ) const;
  int				PlaneSide( const idPlane &plane, const float epsilon = ON_EPSILON ) const;

  bool			ContainsPoint( const idVec3 &p ) const;			// includes touching
  bool			IntersectsBounds( const idBounds &a ) const;	// includes touching
  bool			LineIntersection( const idVec3 &start, const idVec3 &end ) const;
  // intersection point is start + dir * scale
  bool			RayIntersection( const idVec3 &start, const idVec3 &dir, float &scale ) const;

  // most tight bounds for the given transformed bounds
  void			FromTransformedBounds( const idBounds &bounds, const idVec3 &origin, const idMat3 &axis );
  // most tight bounds for a point set
  void			FromPoints( const idVec3 *points, const int numPoints );
  // most tight bounds for a translation
  void			FromPointTranslation( const idVec3 &point, const idVec3 &translation );
  void			FromBoundsTranslation( const idBounds &bounds, const idVec3 &origin, const idMat3 &axis, const idVec3 &translation );
  // most tight bounds for a rotation
  void			FromPointRotation( const idVec3 &point, const idRotation &rotation );
  void			FromBoundsRotation( const idBounds &bounds, const idVec3 &origin, const idMat3 &axis, const idRotation &rotation );

  void			ToPoints( idVec3 points[8] ) const;
  // idSphere		ToSphere( void ) const;

  void			AxisProjection( const idVec3 &dir, float &min, float &max ) const;
  void			AxisProjection( const idVec3 &origin, const idMat3 &axis, const idVec3 &dir, float &min, float &max ) const;

private:
  idVec3			b[2];
};

extern idBounds	bounds_zero;

