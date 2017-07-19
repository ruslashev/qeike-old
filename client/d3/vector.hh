#pragma once

#include <cmath>
#include <cstring>
#include "../math.hh"

/*
   ===============================================================================

   Vector classes

   ===============================================================================
   */

#define VECTOR_EPSILON		0.001f

class idAngles;
class idPolar3;
class idMat3;

//===============================================================
//
//	idVec2 - 2D vector
//
//===============================================================

class idVec2 {
public:
  float			x;
  float			y;

  idVec2(void);
  explicit idVec2( const float x, const float y );

  void			Set(const float x, const float y);
  void			Zero(void);

  float			operator[](int index) const;
  float &			operator[](int index);
  idVec2			operator-() const;
  float			operator*(const idVec2 &a) const;
  idVec2			operator*(const float a) const;
  idVec2			operator/(const float a) const;
  idVec2			operator+(const idVec2 &a) const;
  idVec2			operator-(const idVec2 &a) const;
  idVec2 &		operator+=(const idVec2 &a);
  idVec2 &		operator-=(const idVec2 &a);
  idVec2 &		operator/=(const idVec2 &a);
  idVec2 &		operator/=(const float a);
  idVec2 &		operator*=(const float a);

  friend idVec2	operator*(const float a, const idVec2 b);

  bool			Compare(const idVec2 &a) const;							// exact compare, no epsilon
  bool			Compare(const idVec2 &a, const float epsilon) const;		// compare with epsilon
  bool			operator==(	const idVec2 &a ) const;						// exact compare, no epsilon
  bool			operator!=(	const idVec2 &a ) const;						// exact compare, no epsilon

  float			Length(void) const;
  float			LengthFast(void) const;
  float			LengthSqr(void) const;
  float			Normalize(void);			// returns length
  float			NormalizeFast(void);		// returns length
  idVec2 &		Truncate(float length);	// cap length
  void			Clamp(const idVec2 &min, const idVec2 &max);
  void			Snap(void);				// snap to closest integer value
  void			SnapInt(void);			// snap towards integer (floor)

  int				GetDimension(void) const;

  const float *	ToFloatPtr(void) const;
  float *			ToFloatPtr(void);
  const char *	ToString(int precision = 2) const;

  void			Lerp(const idVec2 &v1, const idVec2 &v2, const float l);
};

extern idVec2 vec2_origin;
#define vec2_zero vec2_origin

