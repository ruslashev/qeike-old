#pragma once

#include "vector.hh"

class idDrawVert {
public:
  idVec3 xyz;
  idVec2 st;
  idVec3 normal;
  idVec3 tangents[2];
  unsigned char color[4];

  float			operator[](const int index) const;
  float &			operator[](const int index);

  void			Clear(void);

  void			Lerp(const idDrawVert &a, const idDrawVert &b, const float f);
  void			LerpAll(const idDrawVert &a, const idDrawVert &b, const float f);

  void			Normalize(void);

  void			SetColor(unsigned int color);
  unsigned int			GetColor(void) const;
};

