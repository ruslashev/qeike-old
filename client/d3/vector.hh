#pragma once

#include <cmath>
#include <cstring>
#include "../math.hh"

#define MAX_WORLD_COORD			( 128 * 1024 )
#define MIN_WORLD_COORD			( -128 * 1024 )
#define MAX_WORLD_SIZE			( MAX_WORLD_COORD - MIN_WORLD_COORD )

#define FLOATSIGNBITSET(f)		((*(const unsigned int *)&(f)) >> 31)
#define FLOATSIGNBITNOTSET(f)	((~(*(const unsigned int *)&(f))) >> 31)
#define FLOATNOTZERO(f)			((*(const unsigned int *)&(f)) & ~(1<<31) )
#define INTSIGNBITSET(i)		(((const unsigned int)(i)) >> 31)
#define INTSIGNBITNOTSET(i)		((~((const unsigned int)(i))) >> 31)

/*
   ===============================================================================

   Vector classes

   ===============================================================================
   */

#define VECTOR_EPSILON		0.001f

class idAngles;
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

//===============================================================
//
//	idVec3 - 3D vector
//
//===============================================================

class idVec3 {
public:
  float			x;
  float			y;
  float			z;

  idVec3(void);
  explicit idVec3(const float x, const float y, const float z);

  void			Set(const float x, const float y, const float z);
  void			Zero(void);

  float			operator[](const int index) const;
  float &			operator[](const int index);
  idVec3			operator-() const;
  idVec3 &		operator=(const idVec3 &a);		// required because of a msvc 6 & 7 bug
  float			operator*(const idVec3 &a) const;
  idVec3			operator*(const float a) const;
  idVec3			operator/(const float a) const;
  idVec3			operator+(const idVec3 &a) const;
  idVec3			operator-(const idVec3 &a) const;
  idVec3 &		operator+=(const idVec3 &a);
  idVec3 &		operator-=(const idVec3 &a);
  idVec3 &		operator/=(const idVec3 &a);
  idVec3 &		operator/=(const float a);
  idVec3 &		operator*=(const float a);

  friend idVec3	operator*(const float a, const idVec3 b);

  bool			Compare(const idVec3 &a) const;							// exact compare, no epsilon
  bool			Compare(const idVec3 &a, const float epsilon) const;		// compare with epsilon
  bool			operator==(	const idVec3 &a ) const;						// exact compare, no epsilon
  bool			operator!=(	const idVec3 &a ) const;						// exact compare, no epsilon

  bool			FixDegenerateNormal(void);	// fix degenerate axial cases
  bool			FixDenormals(void);			// change tiny numbers to zero

  idVec3			Cross(const idVec3 &a) const;
  idVec3 &		Cross(const idVec3 &a, const idVec3 &b);
  float			Length(void) const;
  float			LengthSqr(void) const;
  float			LengthFast(void) const;
  float			Normalize(void);				// returns length
  float			NormalizeFast(void);			// returns length
  idVec3 &		Truncate(float length);		// cap length
  void			Clamp(const idVec3 &min, const idVec3 &max);
  void			Snap(void);					// snap to closest integer value
  void			SnapInt(void);				// snap towards integer (floor)

  int				GetDimension(void) const;

  float			ToYaw(void) const;
  float			ToPitch(void) const;
  idAngles		ToAngles(void) const;
  // idPolar3		ToPolar(void) const;
  idMat3			ToMat3(void) const;		// vector should be normalized
  const idVec2 &	ToVec2(void) const;
  idVec2 &		ToVec2(void);
  const float *	ToFloatPtr(void) const;
  float *			ToFloatPtr(void);
  const char *	ToString(int precision = 2) const;

  void			NormalVectors(idVec3 &left, idVec3 &down) const;	// vector should be normalized
  void			OrthogonalBasis(idVec3 &left, idVec3 &up) const;

  void			ProjectOntoPlane(const idVec3 &normal, const float overBounce = 1.0f);
  bool			ProjectAlongPlane(const idVec3 &normal, const float epsilon, const float overBounce = 1.0f);
  void			ProjectSelfOntoSphere(const float radius);

  void			Lerp(const idVec3 &v1, const idVec3 &v2, const float l);
  void			SLerp(const idVec3 &v1, const idVec3 &v2, const float l);
};

extern idVec3 vec3_origin;
#define vec3_zero vec3_origin

//===============================================================
//
//	idVec4 - 4D vector
//
//===============================================================

class idVec4 {
public:
  float			x;
  float			y;
  float			z;
  float			w;

  idVec4(void);
  explicit idVec4(const float x, const float y, const float z, const float w);

  void			Set(const float x, const float y, const float z, const float w);
  void			Zero(void);

  float			operator[](const int index) const;
  float &			operator[](const int index);
  idVec4			operator-() const;
  float			operator*(const idVec4 &a) const;
  idVec4			operator*(const float a) const;
  idVec4			operator/(const float a) const;
  idVec4			operator+(const idVec4 &a) const;
  idVec4			operator-(const idVec4 &a) const;
  idVec4 &		operator+=(const idVec4 &a);
  idVec4 &		operator-=(const idVec4 &a);
  idVec4 &		operator/=(const idVec4 &a);
  idVec4 &		operator/=(const float a);
  idVec4 &		operator*=(const float a);

  friend idVec4	operator*(const float a, const idVec4 b);

  bool			Compare(const idVec4 &a) const;							// exact compare, no epsilon
  bool			Compare(const idVec4 &a, const float epsilon) const;		// compare with epsilon
  bool			operator==(	const idVec4 &a ) const;						// exact compare, no epsilon
  bool			operator!=(	const idVec4 &a ) const;						// exact compare, no epsilon

  float			Length(void) const;
  float			LengthSqr(void) const;
  float			Normalize(void);			// returns length
  float			NormalizeFast(void);		// returns length

  int				GetDimension(void) const;

  const idVec2 &	ToVec2(void) const;
  idVec2 &		ToVec2(void);
  const idVec3 &	ToVec3(void) const;
  idVec3 &		ToVec3(void);
  const float *	ToFloatPtr(void) const;
  float *			ToFloatPtr(void);
  const char *	ToString(int precision = 2) const;

  void			Lerp(const idVec4 &v1, const idVec4 &v2, const float l);
};

extern idVec4 vec4_origin;
#define vec4_zero vec4_origin

//===============================================================
//
//	idVec5 - 5D vector
//
//===============================================================

class idVec5 {
public:
  float			x;
  float			y;
  float			z;
  float			s;
  float			t;

  idVec5(void);
  explicit idVec5(const idVec3 &xyz, const idVec2 &st);
  explicit idVec5(const float x, const float y, const float z, const float s, const float t);

  float			operator[](int index) const;
  float &			operator[](int index);
  idVec5 &		operator=(const idVec3 &a);

  int				GetDimension(void) const;

  const idVec3 &	ToVec3(void) const;
  idVec3 &		ToVec3(void);
  const float *	ToFloatPtr(void) const;
  float *			ToFloatPtr(void);
  const char *	ToString(int precision = 2) const;

  void			Lerp(const idVec5 &v1, const idVec5 &v2, const float l);
};

extern idVec5 vec5_origin;
#define vec5_zero vec5_origin

//===============================================================
//
//	idVec6 - 6D vector
//
//===============================================================

class idVec6 {
public:
  idVec6(void);
  explicit idVec6(const float *a);
  explicit idVec6(const float a1, const float a2, const float a3, const float a4, const float a5, const float a6);

  void			Set(const float a1, const float a2, const float a3, const float a4, const float a5, const float a6);
  void			Zero(void);

  float			operator[](const int index) const;
  float &			operator[](const int index);
  idVec6			operator-() const;
  idVec6			operator*(const float a) const;
  idVec6			operator/(const float a) const;
  float			operator*(const idVec6 &a) const;
  idVec6			operator-(const idVec6 &a) const;
  idVec6			operator+(const idVec6 &a) const;
  idVec6 &		operator*=(const float a);
  idVec6 &		operator/=(const float a);
  idVec6 &		operator+=(const idVec6 &a);
  idVec6 &		operator-=(const idVec6 &a);

  friend idVec6	operator*(const float a, const idVec6 b);

  bool			Compare(const idVec6 &a) const;							// exact compare, no epsilon
  bool			Compare(const idVec6 &a, const float epsilon) const;		// compare with epsilon
  bool			operator==(	const idVec6 &a ) const;						// exact compare, no epsilon
  bool			operator!=(	const idVec6 &a ) const;						// exact compare, no epsilon

  float			Length(void) const;
  float			LengthSqr(void) const;
  float			Normalize(void);			// returns length
  float			NormalizeFast(void);		// returns length

  int				GetDimension(void) const;

  const idVec3 &	SubVec3(int index) const;
  idVec3 &		SubVec3(int index);
  const float *	ToFloatPtr(void) const;
  float *			ToFloatPtr(void);
  const char *	ToString(int precision = 2) const;

private:
  float			p[6];
};

extern idVec6 vec6_origin;
#define vec6_zero vec6_origin
extern idVec6 vec6_infinity;

//===============================================================
//
//	idVecX - arbitrary sized vector
//
//  The vector lives on 16 byte aligned and 16 byte padded memory.
//
//	NOTE: due to the temporary memory pool idVecX cannot be used by multiple threads
//
//===============================================================

#define VECX_MAX_TEMP		1024
#define VECX_QUAD(x)		((((x) + 3) & ~3) * sizeof(float))
#define VECX_CLEAREND()		int s = size; while(s < ((s + 3) & ~3)) { p[s++] = 0.0f; }
#define VECX_ALLOCA(n)	((float *) _alloca16(VECX_QUAD(n)))
#define VECX_SIMD

class idVecX {
  friend class idMatX;

public:
  idVecX(void);
  explicit idVecX(int length);
  explicit idVecX(int length, float *data);
  ~idVecX(void);

  float			operator[](const int index) const;
  float &			operator[](const int index);
  idVecX			operator-() const;
  idVecX &		operator=(const idVecX &a);
  idVecX			operator*(const float a) const;
  idVecX			operator/(const float a) const;
  float			operator*(const idVecX &a) const;
  idVecX			operator-(const idVecX &a) const;
  idVecX			operator+(const idVecX &a) const;
  idVecX &		operator*=(const float a);
  idVecX &		operator/=(const float a);
  idVecX &		operator+=(const idVecX &a);
  idVecX &		operator-=(const idVecX &a);

  friend idVecX	operator*(const float a, const idVecX b);

  bool			Compare(const idVecX &a) const;							// exact compare, no epsilon
  bool			Compare(const idVecX &a, const float epsilon) const;		// compare with epsilon
  bool			operator==(	const idVecX &a ) const;						// exact compare, no epsilon
  bool			operator!=(	const idVecX &a ) const;						// exact compare, no epsilon

  void			SetSize(int size);
  void			ChangeSize(int size, bool makeZero = false);
  int				GetSize(void) const { return size; }
  void			SetData(int length, float *data);
  void			Zero(void);
  void			Zero(int length);
  void			Negate(void);
  void			Clamp(float min, float max);
  idVecX &		SwapElements(int e1, int e2);

  float			Length(void) const;
  float			LengthSqr(void) const;
  idVecX			Normalize(void) const;
  float			NormalizeSelf(void);

  int				GetDimension(void) const;

  const idVec3 &	SubVec3(int index) const;
  idVec3 &		SubVec3(int index);
  const idVec6 &	SubVec6(int index) const;
  idVec6 &		SubVec6(int index);
  const float *	ToFloatPtr(void) const;
  float *			ToFloatPtr(void);
  const char *	ToString(int precision = 2) const;

private:
  int				size;					// size of the vector
  int				alloced;				// if -1 p points to data set with SetData
  float *			p;						// memory the vector is stored

  static float	temp[VECX_MAX_TEMP+4];	// used to store intermediate results
  static float *	tempPtr;				// pointer to 16 byte aligned temporary memory
  static int		tempIndex;				// index into memory pool, wraps around

private:
  void			SetTempSize(int size);
};


/*
   ===============================================================================

   Old 3D vector macros, should no longer be used.

   ===============================================================================
   */

#define DotProduct(a, b)			((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])
#define VectorSubtract(a, b, c)	((c)[0]=(a)[0]-(b)[0],(c)[1]=(a)[1]-(b)[1],(c)[2]=(a)[2]-(b)[2])
#define VectorAdd(a, b, c)		((c)[0]=(a)[0]+(b)[0],(c)[1]=(a)[1]+(b)[1],(c)[2]=(a)[2]+(b)[2])
#define	VectorScale(v, s, o)		((o)[0]=(v)[0]*(s),(o)[1]=(v)[1]*(s),(o)[2]=(v)[2]*(s))
#define	VectorMA(v, s, b, o)		((o)[0]=(v)[0]+(b)[0]*(s),(o)[1]=(v)[1]+(b)[1]*(s),(o)[2]=(v)[2]+(b)[2]*(s))
#define VectorCopy(a, b)			((b)[0]=(a)[0],(b)[1]=(a)[1],(b)[2]=(a)[2])
