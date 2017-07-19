#include "drawvert.hh"

float idDrawVert::operator[](const int index) const {
  assert(index >= 0 && index < 5);
  return ((float *)(&xyz))[index];
}
float	&idDrawVert::operator[](const int index) {
  assert(index >= 0 && index < 5);
  return ((float *)(&xyz))[index];
}

void idDrawVert::Clear(void) {
  xyz.Zero();
  st.Zero();
  normal.Zero();
  tangents[0].Zero();
  tangents[1].Zero();
  color[0] = color[1] = color[2] = color[3] = 0;
}

void idDrawVert::Lerp(const idDrawVert &a, const idDrawVert &b, const float f) {
  xyz = a.xyz + f * (b.xyz - a.xyz);
  st = a.st + f * (b.st - a.st);
}

void idDrawVert::LerpAll(const idDrawVert &a, const idDrawVert &b, const float f) {
  xyz = a.xyz + f * (b.xyz - a.xyz);
  st = a.st + f * (b.st - a.st);
  normal = a.normal + f * (b.normal - a.normal);
  tangents[0] = a.tangents[0] + f * (b.tangents[0] - a.tangents[0]);
  tangents[1] = a.tangents[1] + f * (b.tangents[1] - a.tangents[1]);
  color[0] = (unsigned char)(a.color[0] + f * (b.color[0] - a.color[0]));
  color[1] = (unsigned char)(a.color[1] + f * (b.color[1] - a.color[1]));
  color[2] = (unsigned char)(a.color[2] + f * (b.color[2] - a.color[2]));
  color[3] = (unsigned char)(a.color[3] + f * (b.color[3] - a.color[3]));
}

void idDrawVert::SetColor(unsigned int color) {
  *reinterpret_cast<unsigned int *>(this->color) = color;
}

unsigned int idDrawVert::GetColor(void) const {
  return *reinterpret_cast<const unsigned int *>(this->color);
}

/*
   =============
   idDrawVert::Normalize
   =============
   */
void idDrawVert::Normalize(void) {
  normal.Normalize();
  tangents[1].Cross(normal, tangents[0]);
  tangents[1].Normalize();
  tangents[0].Cross(tangents[1], normal);
  tangents[0].Normalize();
}

