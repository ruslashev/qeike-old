#include "winding.hh"

idWinding::idWinding( void ) {
  numPoints = allocedSize = 0;
  p = NULL;
}

idWinding::idWinding( int n ) {
  numPoints = allocedSize = 0;
  p = NULL;
  EnsureAlloced( n );
}

idWinding::idWinding( const idVec3 *verts, const int n ) {
  int i;

  numPoints = allocedSize = 0;
  p = NULL;
  if ( !EnsureAlloced( n ) ) {
    numPoints = 0;
    return;
  }
  for ( i = 0; i < n; i++ ) {
    p[i].ToVec3() = verts[i];
    p[i].s = p[i].t = 0.0f;
  }
  numPoints = n;
}

idWinding::idWinding( const idVec3 &normal, const float dist ) {
  numPoints = allocedSize = 0;
  p = NULL;
  BaseForPlane( normal, dist );
}

idWinding::idWinding( const idPlane &plane ) {
  numPoints = allocedSize = 0;
  p = NULL;
  BaseForPlane( plane );
}

idWinding::idWinding( const idWinding &winding ) {
  int i;
  if ( !EnsureAlloced( winding.GetNumPoints() ) ) {
    numPoints = 0;
    return;
  }
  for ( i = 0; i < winding.GetNumPoints(); i++ ) {
    p[i] = winding[i];
  }
  numPoints = winding.GetNumPoints();
}

idWinding::~idWinding( void ) {
  delete[] p;
  p = NULL;
}

idWinding &idWinding::operator=( const idWinding &winding ) {
  int i;

  if ( !EnsureAlloced( winding.numPoints ) ) {
    numPoints = 0;
    return *this;
  }
  for ( i = 0; i < winding.numPoints; i++ ) {
    p[i] = winding.p[i];
  }
  numPoints = winding.numPoints;
  return *this;
}

const idVec5 &idWinding::operator[]( const int index ) const {
  //assert( index >= 0 && index < numPoints );
  return p[ index ];
}

idVec5 &idWinding::operator[]( const int index ) {
  //assert( index >= 0 && index < numPoints );
  return p[ index ];
}

idWinding &idWinding::operator+=( const idVec3 &v ) {
  AddPoint( v );
  return *this;
}

idWinding &idWinding::operator+=( const idVec5 &v ) {
  AddPoint( v );
  return *this;
}

void idWinding::AddPoint( const idVec3 &v ) {
  if ( !EnsureAlloced(numPoints+1, true) ) {
    return;
  }
  p[numPoints] = v;
  numPoints++;
}

void idWinding::AddPoint( const idVec5 &v ) {
  if ( !EnsureAlloced(numPoints+1, true) ) {
    return;
  }
  p[numPoints] = v;
  numPoints++;
}

int idWinding::GetNumPoints( void ) const {
  return numPoints;
}

void idWinding::SetNumPoints( int n ) {
  if ( !EnsureAlloced( n, true ) ) {
    return;
  }
  numPoints = n;
}

void idWinding::Clear( void ) {
  numPoints = 0;
  delete[] p;
  p = NULL;
}

void idWinding::BaseForPlane( const idPlane &plane ) {
  BaseForPlane( plane.Normal(), plane.Dist() );
}

bool idWinding::EnsureAlloced( int n, bool keep ) {
  if ( n > allocedSize ) {
    return ReAllocate( n, keep );
  }
  return true;
}

