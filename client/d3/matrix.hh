#pragma once

#include "vector.hh"

#define MATRIX_INVERSE_EPSILON		1e-14
#define MATRIX_EPSILON				1e-6

class idAngles;
class idQuat;
class idCQuat;
class idRotation;
class idMat4;

//===============================================================
//
//	idMat2 - 2x2 matrix
//
//===============================================================

class idMat2 {
public:
  idMat2( void );
  explicit idMat2( const idVec2 &x, const idVec2 &y );
  explicit idMat2( const float xx, const float xy, const float yx, const float yy );
  explicit idMat2( const float src[ 2 ][ 2 ] );

  const idVec2 &	operator[]( int index ) const;
  idVec2 &		operator[]( int index );
  idMat2			operator-() const;
  idMat2			operator*( const float a ) const;
  idVec2			operator*( const idVec2 &vec ) const;
  idMat2			operator*( const idMat2 &a ) const;
  idMat2			operator+( const idMat2 &a ) const;
  idMat2			operator-( const idMat2 &a ) const;
  idMat2 &		operator*=( const float a );
  idMat2 &		operator*=( const idMat2 &a );
  idMat2 &		operator+=( const idMat2 &a );
  idMat2 &		operator-=( const idMat2 &a );

  friend idMat2	operator*( const float a, const idMat2 &mat );
  friend idVec2	operator*( const idVec2 &vec, const idMat2 &mat );
  friend idVec2 &	operator*=( idVec2 &vec, const idMat2 &mat );

  bool			Compare( const idMat2 &a ) const;						// exact compare, no epsilon
  bool			Compare( const idMat2 &a, const float epsilon ) const;	// compare with epsilon
  bool			operator==( const idMat2 &a ) const;					// exact compare, no epsilon
  bool			operator!=( const idMat2 &a ) const;					// exact compare, no epsilon

  void			Zero( void );
  void			Identity( void );
  bool			IsIdentity( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsSymmetric( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsDiagonal( const float epsilon = MATRIX_EPSILON ) const;

  float			Trace( void ) const;
  float			Determinant( void ) const;
  idMat2			Transpose( void ) const;	// returns transpose
  idMat2 &		TransposeSelf( void );
  idMat2			Inverse( void ) const;		// returns the inverse ( m * m.Inverse() = identity )
  bool			InverseSelf( void );		// returns false if determinant is zero
  idMat2			InverseFast( void ) const;	// returns the inverse ( m * m.Inverse() = identity )
  bool			InverseFastSelf( void );	// returns false if determinant is zero

  int				GetDimension( void ) const;

  const float *	ToFloatPtr( void ) const;
  float *			ToFloatPtr( void );
  const char *	ToString( int precision = 2 ) const;

private:
  idVec2			mat[ 2 ];
};

extern idMat2 mat2_zero;
extern idMat2 mat2_identity;
#define mat2_default	mat2_identity

//===============================================================
//
//	idMat3 - 3x3 matrix
//
//	NOTE:	matrix is column-major
//
//===============================================================

class idMat3 {
public:
  idMat3( void );
  explicit idMat3( const idVec3 &x, const idVec3 &y, const idVec3 &z );
  explicit idMat3( const float xx, const float xy, const float xz, const float yx, const float yy, const float yz, const float zx, const float zy, const float zz );
  explicit idMat3( const float src[ 3 ][ 3 ] );

  const idVec3 &	operator[]( int index ) const;
  idVec3 &		operator[]( int index );
  idMat3			operator-() const;
  idMat3			operator*( const float a ) const;
  idVec3			operator*( const idVec3 &vec ) const;
  idMat3			operator*( const idMat3 &a ) const;
  idMat3			operator+( const idMat3 &a ) const;
  idMat3			operator-( const idMat3 &a ) const;
  idMat3 &		operator*=( const float a );
  idMat3 &		operator*=( const idMat3 &a );
  idMat3 &		operator+=( const idMat3 &a );
  idMat3 &		operator-=( const idMat3 &a );

  friend idMat3	operator*( const float a, const idMat3 &mat );
  friend idVec3	operator*( const idVec3 &vec, const idMat3 &mat );
  friend idVec3 &	operator*=( idVec3 &vec, const idMat3 &mat );

  bool			Compare( const idMat3 &a ) const;						// exact compare, no epsilon
  bool			Compare( const idMat3 &a, const float epsilon ) const;	// compare with epsilon
  bool			operator==( const idMat3 &a ) const;					// exact compare, no epsilon
  bool			operator!=( const idMat3 &a ) const;					// exact compare, no epsilon

  void			Zero( void );
  void			Identity( void );
  bool			IsIdentity( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsSymmetric( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsDiagonal( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsRotated( void ) const;

  void			ProjectVector( const idVec3 &src, idVec3 &dst ) const;
  void			UnprojectVector( const idVec3 &src, idVec3 &dst ) const;

  bool			FixDegeneracies( void );	// fix degenerate axial cases
  bool			FixDenormals( void );		// change tiny numbers to zero

  float			Trace( void ) const;
  float			Determinant( void ) const;
  idMat3			OrthoNormalize( void ) const;
  idMat3 &		OrthoNormalizeSelf( void );
  idMat3			Transpose( void ) const;	// returns transpose
  idMat3 &		TransposeSelf( void );
  idMat3			Inverse( void ) const;		// returns the inverse ( m * m.Inverse() = identity )
  bool			InverseSelf( void );		// returns false if determinant is zero
  idMat3			InverseFast( void ) const;	// returns the inverse ( m * m.Inverse() = identity )
  bool			InverseFastSelf( void );	// returns false if determinant is zero
  idMat3			TransposeMultiply( const idMat3 &b ) const;

  idMat3			InertiaTranslate( const float mass, const idVec3 &centerOfMass, const idVec3 &translation ) const;
  idMat3 &		InertiaTranslateSelf( const float mass, const idVec3 &centerOfMass, const idVec3 &translation );
  idMat3			InertiaRotate( const idMat3 &rotation ) const;
  idMat3 &		InertiaRotateSelf( const idMat3 &rotation );

  int				GetDimension( void ) const;

  idAngles		ToAngles( void ) const;
  idQuat			ToQuat( void ) const;
  idCQuat			ToCQuat( void ) const;
  idRotation		ToRotation( void ) const;
  idMat4			ToMat4( void ) const;
  idVec3			ToAngularVelocity( void ) const;
  const float *	ToFloatPtr( void ) const;
  float *			ToFloatPtr( void );
  const char *	ToString( int precision = 2 ) const;

  friend void		TransposeMultiply( const idMat3 &inv, const idMat3 &b, idMat3 &dst );
  friend idMat3	SkewSymmetric( idVec3 const &src );

private:
  idVec3			mat[ 3 ];
};

extern idMat3 mat3_zero;
extern idMat3 mat3_identity;
#define mat3_default	mat3_identity


//===============================================================
//
//	idMat4 - 4x4 matrix
//
//===============================================================

class idMat4 {
public:
  idMat4( void );
  explicit idMat4( const idVec4 &x, const idVec4 &y, const idVec4 &z, const idVec4 &w );
  explicit idMat4(const float xx, const float xy, const float xz, const float xw,
      const float yx, const float yy, const float yz, const float yw,
      const float zx, const float zy, const float zz, const float zw,
      const float wx, const float wy, const float wz, const float ww );
  explicit idMat4( const idMat3 &rotation, const idVec3 &translation );
  explicit idMat4( const float src[ 4 ][ 4 ] );

  const idVec4 &	operator[]( int index ) const;
  idVec4 &		operator[]( int index );
  idMat4			operator*( const float a ) const;
  idVec4			operator*( const idVec4 &vec ) const;
  idVec3			operator*( const idVec3 &vec ) const;
  idMat4			operator*( const idMat4 &a ) const;
  idMat4			operator+( const idMat4 &a ) const;
  idMat4			operator-( const idMat4 &a ) const;
  idMat4 &		operator*=( const float a );
  idMat4 &		operator*=( const idMat4 &a );
  idMat4 &		operator+=( const idMat4 &a );
  idMat4 &		operator-=( const idMat4 &a );

  friend idMat4	operator*( const float a, const idMat4 &mat );
  friend idVec4	operator*( const idVec4 &vec, const idMat4 &mat );
  friend idVec3	operator*( const idVec3 &vec, const idMat4 &mat );
  friend idVec4 &	operator*=( idVec4 &vec, const idMat4 &mat );
  friend idVec3 &	operator*=( idVec3 &vec, const idMat4 &mat );

  bool			Compare( const idMat4 &a ) const;						// exact compare, no epsilon
  bool			Compare( const idMat4 &a, const float epsilon ) const;	// compare with epsilon
  bool			operator==( const idMat4 &a ) const;					// exact compare, no epsilon
  bool			operator!=( const idMat4 &a ) const;					// exact compare, no epsilon

  void			Zero( void );
  void			Identity( void );
  bool			IsIdentity( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsSymmetric( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsDiagonal( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsRotated( void ) const;

  void			ProjectVector( const idVec4 &src, idVec4 &dst ) const;
  void			UnprojectVector( const idVec4 &src, idVec4 &dst ) const;

  float			Trace( void ) const;
  float			Determinant( void ) const;
  idMat4			Transpose( void ) const;	// returns transpose
  idMat4 &		TransposeSelf( void );
  idMat4			Inverse( void ) const;		// returns the inverse ( m * m.Inverse() = identity )
  bool			InverseSelf( void );		// returns false if determinant is zero
  idMat4			InverseFast( void ) const;	// returns the inverse ( m * m.Inverse() = identity )
  bool			InverseFastSelf( void );	// returns false if determinant is zero
  idMat4			TransposeMultiply( const idMat4 &b ) const;

  int				GetDimension( void ) const;

  const float *	ToFloatPtr( void ) const;
  float *			ToFloatPtr( void );
  const char *	ToString( int precision = 2 ) const;

private:
  idVec4			mat[ 4 ];
};

extern idMat4 mat4_zero;
extern idMat4 mat4_identity;
#define mat4_default	mat4_identity

//===============================================================
//
//	idMat5 - 5x5 matrix
//
//===============================================================

class idMat5 {
public:
  idMat5( void );
  explicit idMat5( const idVec5 &v0, const idVec5 &v1, const idVec5 &v2, const idVec5 &v3, const idVec5 &v4 );
  explicit idMat5( const float src[ 5 ][ 5 ] );

  const idVec5 &	operator[]( int index ) const;
  idVec5 &		operator[]( int index );
  idMat5			operator*( const float a ) const;
  idVec5			operator*( const idVec5 &vec ) const;
  idMat5			operator*( const idMat5 &a ) const;
  idMat5			operator+( const idMat5 &a ) const;
  idMat5			operator-( const idMat5 &a ) const;
  idMat5 &		operator*=( const float a );
  idMat5 &		operator*=( const idMat5 &a );
  idMat5 &		operator+=( const idMat5 &a );
  idMat5 &		operator-=( const idMat5 &a );

  friend idMat5	operator*( const float a, const idMat5 &mat );
  friend idVec5	operator*( const idVec5 &vec, const idMat5 &mat );
  friend idVec5 &	operator*=( idVec5 &vec, const idMat5 &mat );

  bool			Compare( const idMat5 &a ) const;						// exact compare, no epsilon
  bool			Compare( const idMat5 &a, const float epsilon ) const;	// compare with epsilon
  bool			operator==( const idMat5 &a ) const;					// exact compare, no epsilon
  bool			operator!=( const idMat5 &a ) const;					// exact compare, no epsilon

  void			Zero( void );
  void			Identity( void );
  bool			IsIdentity( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsSymmetric( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsDiagonal( const float epsilon = MATRIX_EPSILON ) const;

  float			Trace( void ) const;
  float			Determinant( void ) const;
  idMat5			Transpose( void ) const;	// returns transpose
  idMat5 &		TransposeSelf( void );
  idMat5			Inverse( void ) const;		// returns the inverse ( m * m.Inverse() = identity )
  bool			InverseSelf( void );		// returns false if determinant is zero
  idMat5			InverseFast( void ) const;	// returns the inverse ( m * m.Inverse() = identity )
  bool			InverseFastSelf( void );	// returns false if determinant is zero

  int				GetDimension( void ) const;

  const float *	ToFloatPtr( void ) const;
  float *			ToFloatPtr( void );
  const char *	ToString( int precision = 2 ) const;

private:
  idVec5			mat[ 5 ];
};

extern idMat5 mat5_zero;
extern idMat5 mat5_identity;
#define mat5_default	mat5_identity


//===============================================================
//
//	idMat6 - 6x6 matrix
//
//===============================================================

class idMat6 {
public:
  idMat6( void );
  explicit idMat6( const idVec6 &v0, const idVec6 &v1, const idVec6 &v2, const idVec6 &v3, const idVec6 &v4, const idVec6 &v5 );
  explicit idMat6( const idMat3 &m0, const idMat3 &m1, const idMat3 &m2, const idMat3 &m3 );
  explicit idMat6( const float src[ 6 ][ 6 ] );

  const idVec6 &	operator[]( int index ) const;
  idVec6 &		operator[]( int index );
  idMat6			operator*( const float a ) const;
  idVec6			operator*( const idVec6 &vec ) const;
  idMat6			operator*( const idMat6 &a ) const;
  idMat6			operator+( const idMat6 &a ) const;
  idMat6			operator-( const idMat6 &a ) const;
  idMat6 &		operator*=( const float a );
  idMat6 &		operator*=( const idMat6 &a );
  idMat6 &		operator+=( const idMat6 &a );
  idMat6 &		operator-=( const idMat6 &a );

  friend idMat6	operator*( const float a, const idMat6 &mat );
  friend idVec6	operator*( const idVec6 &vec, const idMat6 &mat );
  friend idVec6 &	operator*=( idVec6 &vec, const idMat6 &mat );

  bool			Compare( const idMat6 &a ) const;						// exact compare, no epsilon
  bool			Compare( const idMat6 &a, const float epsilon ) const;	// compare with epsilon
  bool			operator==( const idMat6 &a ) const;					// exact compare, no epsilon
  bool			operator!=( const idMat6 &a ) const;					// exact compare, no epsilon

  void			Zero( void );
  void			Identity( void );
  bool			IsIdentity( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsSymmetric( const float epsilon = MATRIX_EPSILON ) const;
  bool			IsDiagonal( const float epsilon = MATRIX_EPSILON ) const;

  idMat3			SubMat3( int n ) const;
  float			Trace( void ) const;
  float			Determinant( void ) const;
  idMat6			Transpose( void ) const;	// returns transpose
  idMat6 &		TransposeSelf( void );
  idMat6			Inverse( void ) const;		// returns the inverse ( m * m.Inverse() = identity )
  bool			InverseSelf( void );		// returns false if determinant is zero
  idMat6			InverseFast( void ) const;	// returns the inverse ( m * m.Inverse() = identity )
  bool			InverseFastSelf( void );	// returns false if determinant is zero

  int				GetDimension( void ) const;

  const float *	ToFloatPtr( void ) const;
  float *			ToFloatPtr( void );
  const char *	ToString( int precision = 2 ) const;

private:
  idVec6			mat[ 6 ];
};

extern idMat6 mat6_zero;
extern idMat6 mat6_identity;
#define mat6_default	mat6_identity

