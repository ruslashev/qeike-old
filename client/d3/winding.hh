#pragma once

#include "vector.hh"
#include "plane.hh"
#include "bounds.hh"

class idWinding {

public:
					idWinding( void );
					explicit idWinding( const int n );								// allocate for n points
					explicit idWinding( const idVec3 *verts, const int n );			// winding from points
					explicit idWinding( const idVec3 &normal, const float dist );	// base winding for plane
					explicit idWinding( const idPlane &plane );						// base winding for plane
					explicit idWinding( const idWinding &winding );
	virtual			~idWinding( void );

	idWinding &		operator=( const idWinding &winding );
	const idVec5 &	operator[]( const int index ) const;
	idVec5 &		operator[]( const int index );

					// add a point to the end of the winding point array
	idWinding &		operator+=( const idVec3 &v );
	idWinding &		operator+=( const idVec5 &v );
	void			AddPoint( const idVec3 &v );
	void			AddPoint( const idVec5 &v );

					// number of points on winding
	int				GetNumPoints( void ) const;
	void			SetNumPoints( int n );
	virtual void	Clear( void );

					// huge winding for plane, the points go counter clockwise when facing the front of the plane
	void			BaseForPlane( const idVec3 &normal, const float dist );
	void			BaseForPlane( const idPlane &plane );

					// splits the winding into a front and back winding, the winding itself stays unchanged
					// returns a SIDE_?
	int				Split( const idPlane &plane, const float epsilon, idWinding **front, idWinding **back ) const;
					// returns the winding fragment at the front of the clipping plane,
					// if there is nothing at the front the winding itself is destroyed and NULL is returned
	idWinding *		Clip( const idPlane &plane, const float epsilon = ON_EPSILON, const bool keepOn = false );
					// cuts off the part at the back side of the plane, returns true if some part was at the front
					// if there is nothing at the front the number of points is set to zero
	bool			ClipInPlace( const idPlane &plane, const float epsilon = ON_EPSILON, const bool keepOn = false );

					// returns a copy of the winding
	idWinding *		Copy( void ) const;
	idWinding *		Reverse( void ) const;
	void			ReverseSelf( void );
	void			RemoveEqualPoints( const float epsilon = ON_EPSILON );
	void			RemoveColinearPoints( const idVec3 &normal, const float epsilon = ON_EPSILON );
	void			RemovePoint( int point );
	void			InsertPoint( const idVec3 &point, int spot );
	bool			InsertPointIfOnEdge( const idVec3 &point, const idPlane &plane, const float epsilon = ON_EPSILON );
					// add a winding to the convex hull
	void			AddToConvexHull( const idWinding *winding, const idVec3 &normal, const float epsilon = ON_EPSILON );
					// add a point to the convex hull
	void			AddToConvexHull( const idVec3 &point, const idVec3 &normal, const float epsilon = ON_EPSILON );
					// tries to merge 'this' with the given winding, returns NULL if merge fails, both 'this' and 'w' stay intact
					// 'keep' tells if the contacting points should stay even if they create colinear edges
	idWinding *		TryMerge( const idWinding &w, const idVec3 &normal, int keep = false ) const;
					// check whether the winding is valid or not
	bool			Check( bool print = true ) const;

	float			GetArea( void ) const;
	idVec3			GetCenter( void ) const;
	float			GetRadius( const idVec3 &center ) const;
	void			GetPlane( idVec3 &normal, float &dist ) const;
	void			GetPlane( idPlane &plane ) const;
	void			GetBounds( idBounds &bounds ) const;

	bool			IsTiny( void ) const;
	bool			IsHuge( void ) const;	// base winding for a plane is typically huge
	void			Print( void ) const;

	float			PlaneDistance( const idPlane &plane ) const;
	int				PlaneSide( const idPlane &plane, const float epsilon = ON_EPSILON ) const;

	bool			PlanesConcave( const idWinding &w2, const idVec3 &normal1, const idVec3 &normal2, float dist1, float dist2 ) const;

	bool			PointInside( const idVec3 &normal, const idVec3 &point, const float epsilon ) const;
					// returns true if the line or ray intersects the winding
	bool			LineIntersection( const idPlane &windingPlane, const idVec3 &start, const idVec3 &end, bool backFaceCull = false ) const;
					// intersection point is start + dir * scale
	bool			RayIntersection( const idPlane &windingPlane, const idVec3 &start, const idVec3 &dir, float &scale, bool backFaceCull = false ) const;

	static float	TriangleArea( const idVec3 &a, const idVec3 &b, const idVec3 &c );

protected:
	int				numPoints;				// number of points
	idVec5 *		p;						// pointer to point data
	int				allocedSize;

	bool			EnsureAlloced( int n, bool keep = false );
	virtual bool	ReAllocate( int n, bool keep = false );
};

/*
===============================================================================

	idFixedWinding is a fixed buffer size winding not using
	memory allocations.

	When an operation would overflow the fixed buffer a warning
	is printed and the operation is safely cancelled.

===============================================================================
*/

#define	MAX_POINTS_ON_WINDING	64

class idFixedWinding : public idWinding {

public:
					idFixedWinding( void );
					explicit idFixedWinding( const int n );
					explicit idFixedWinding( const idVec3 *verts, const int n );
					explicit idFixedWinding( const idVec3 &normal, const float dist );
					explicit idFixedWinding( const idPlane &plane );
					explicit idFixedWinding( const idWinding &winding );
					explicit idFixedWinding( const idFixedWinding &winding );
	virtual			~idFixedWinding( void );

	idFixedWinding &operator=( const idWinding &winding );

	virtual void	Clear( void );

					// splits the winding in a back and front part, 'this' becomes the front part
					// returns a SIDE_?
	int				Split( idFixedWinding *back, const idPlane &plane, const float epsilon = ON_EPSILON );

protected:
	idVec5			data[MAX_POINTS_ON_WINDING];	// point data

	virtual bool	ReAllocate( int n, bool keep = false );
};

