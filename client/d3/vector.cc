#include "vector.hh"


idVec2::idVec2(void) {
}

idVec2::idVec2(const float x, const float y) {
  this->x = x;
  this->y = y;
}

void idVec2::Set(const float x, const float y) {
  this->x = x;
  this->y = y;
}

void idVec2::Zero(void) {
  x = y = 0.0f;
}

bool idVec2::Compare(const idVec2 &a) const {
  return ((x == a.x) && (y == a.y));
}

bool idVec2::Compare(const idVec2 &a, const float epsilon) const {
  if (qkmath::fabs(x - a.x) > epsilon) {
    return false;
  }

  if (qkmath::fabs(y - a.y) > epsilon) {
    return false;
  }

  return true;
}

bool idVec2::operator==(const idVec2 &a) const {
  return Compare(a);
}

bool idVec2::operator!=(const idVec2 &a) const {
  return !Compare(a);
}

float idVec2::operator[](int index) const {
  return (&x)[ index ];
}

float& idVec2::operator[](int index) {
  return (&x)[ index ];
}

float idVec2::Length(void) const {
  return (float)sqrtf(x * x + y * y);
}

float idVec2::LengthFast(void) const {
  float sqrLength;

  sqrLength = x * x + y * y;
  return sqrLength * qkmath::fastinvsqrt(sqrLength);
}

float idVec2::LengthSqr(void) const {
  return (x * x + y * y);
}

float idVec2::Normalize(void) {
  float sqrLength, invLength;

  sqrLength = x * x + y * y;
  invLength = qkmath::inv_sqrt(sqrLength);
  x *= invLength;
  y *= invLength;
  return invLength * sqrLength;
}

float idVec2::NormalizeFast(void) {
  float lengthSqr, invLength;

  lengthSqr = x * x + y * y;
  invLength = qkmath::fastinvsqrt(lengthSqr);
  x *= invLength;
  y *= invLength;
  return invLength * lengthSqr;
}

idVec2 &idVec2::Truncate(float length) {
  float length2;
  float ilength;

  if (!length) {
    Zero();
  }
  else {
    length2 = LengthSqr();
    if (length2 > length * length) {
      ilength = length * qkmath::inv_sqrt(length2);
      x *= ilength;
      y *= ilength;
    }
  }

  return *this;
}

void idVec2::Clamp(const idVec2 &min, const idVec2 &max) {
  if (x < min.x) {
    x = min.x;
  } else if (x > max.x) {
    x = max.x;
  }
  if (y < min.y) {
    y = min.y;
  } else if (y > max.y) {
    y = max.y;
  }
}

void idVec2::Snap(void) {
  x = floor(x + 0.5f);
  y = floor(y + 0.5f);
}

void idVec2::SnapInt(void) {
  x = float(int(x));
  y = float(int(y));
}

idVec2 idVec2::operator-() const {
  return idVec2(-x, -y);
}

idVec2 idVec2::operator-(const idVec2 &a) const {
  return idVec2(x - a.x, y - a.y);
}

float idVec2::operator*(const idVec2 &a) const {
  return x * a.x + y * a.y;
}

idVec2 idVec2::operator*(const float a) const {
  return idVec2(x * a, y * a);
}

idVec2 idVec2::operator/(const float a) const {
  float inva = 1.0f / a;
  return idVec2(x * inva, y * inva);
}

idVec2 operator*(const float a, const idVec2 b) {
  return idVec2(b.x * a, b.y * a);
}

idVec2 idVec2::operator+(const idVec2 &a) const {
  return idVec2(x + a.x, y + a.y);
}

idVec2 &idVec2::operator+=(const idVec2 &a) {
  x += a.x;
  y += a.y;

  return *this;
}

idVec2 &idVec2::operator/=(const idVec2 &a) {
  x /= a.x;
  y /= a.y;

  return *this;
}

idVec2 &idVec2::operator/=(const float a) {
  float inva = 1.0f / a;
  x *= inva;
  y *= inva;

  return *this;
}

idVec2 &idVec2::operator-=(const idVec2 &a) {
  x -= a.x;
  y -= a.y;

  return *this;
}

idVec2 &idVec2::operator*=(const float a) {
  x *= a;
  y *= a;

  return *this;
}

int idVec2::GetDimension(void) const {
  return 2;
}

const float *idVec2::ToFloatPtr(void) const {
  return &x;
}

float *idVec2::ToFloatPtr(void) {
  return &x;
}


idVec3::idVec3(void) {
}

idVec3::idVec3(const float x, const float y, const float z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

float idVec3::operator[](const int index) const {
  return (&x)[ index ];
}

float &idVec3::operator[](const int index) {
  return (&x)[ index ];
}

void idVec3::Set(const float x, const float y, const float z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

void idVec3::Zero(void) {
  x = y = z = 0.0f;
}

idVec3 idVec3::operator-() const {
  return idVec3(-x, -y, -z);
}

idVec3 &idVec3::operator=(const idVec3 &a) {
  x = a.x;
  y = a.y;
  z = a.z;
  return *this;
}

idVec3 idVec3::operator-(const idVec3 &a) const {
  return idVec3(x - a.x, y - a.y, z - a.z);
}

float idVec3::operator*(const idVec3 &a) const {
  return x * a.x + y * a.y + z * a.z;
}

idVec3 idVec3::operator*(const float a) const {
  return idVec3(x * a, y * a, z * a);
}

idVec3 idVec3::operator/(const float a) const {
  float inva = 1.0f / a;
  return idVec3(x * inva, y * inva, z * inva);
}

idVec3 operator*(const float a, const idVec3 b) {
  return idVec3(b.x * a, b.y * a, b.z * a);
}

idVec3 idVec3::operator+(const idVec3 &a) const {
  return idVec3(x + a.x, y + a.y, z + a.z);
}

idVec3 &idVec3::operator+=(const idVec3 &a) {
  x += a.x;
  y += a.y;
  z += a.z;

  return *this;
}

idVec3 &idVec3::operator/=(const idVec3 &a) {
  x /= a.x;
  y /= a.y;
  z /= a.z;

  return *this;
}

idVec3 &idVec3::operator/=(const float a) {
  float inva = 1.0f / a;
  x *= inva;
  y *= inva;
  z *= inva;

  return *this;
}

idVec3 &idVec3::operator-=(const idVec3 &a) {
  x -= a.x;
  y -= a.y;
  z -= a.z;

  return *this;
}

idVec3 &idVec3::operator*=(const float a) {
  x *= a;
  y *= a;
  z *= a;

  return *this;
}

bool idVec3::Compare(const idVec3 &a) const {
  return ((x == a.x) && (y == a.y) && (z == a.z));
}

bool idVec3::Compare(const idVec3 &a, const float epsilon) const {
  if (qkmath::fabs(x - a.x) > epsilon) {
    return false;
  }

  if (qkmath::fabs(y - a.y) > epsilon) {
    return false;
  }

  if (qkmath::fabs(z - a.z) > epsilon) {
    return false;
  }

  return true;
}

bool idVec3::operator==(const idVec3 &a) const {
  return Compare(a);
}

bool idVec3::operator!=(const idVec3 &a) const {
  return !Compare(a);
}

float idVec3::NormalizeFast(void) {
  float sqrLength, invLength;

  sqrLength = x * x + y * y + z * z;
  invLength = qkmath::fastinvsqrt(sqrLength);
  x *= invLength;
  y *= invLength;
  z *= invLength;
  return invLength * sqrLength;
}

bool idVec3::FixDegenerateNormal(void) {
  if (x == 0.0f) {
    if (y == 0.0f) {
      if (z > 0.0f) {
        if (z != 1.0f) {
          z = 1.0f;
          return true;
        }
      } else {
        if (z != -1.0f) {
          z = -1.0f;
          return true;
        }
      }
      return false;
    } else if (z == 0.0f) {
      if (y > 0.0f) {
        if (y != 1.0f) {
          y = 1.0f;
          return true;
        }
      } else {
        if (y != -1.0f) {
          y = -1.0f;
          return true;
        }
      }
      return false;
    }
  } else if (y == 0.0f) {
    if (z == 0.0f) {
      if (x > 0.0f) {
        if (x != 1.0f) {
          x = 1.0f;
          return true;
        }
      } else {
        if (x != -1.0f) {
          x = -1.0f;
          return true;
        }
      }
      return false;
    }
  }
  if (qkmath::fabs(x) == 1.0f) {
    if (y != 0.0f || z != 0.0f) {
      y = z = 0.0f;
      return true;
    }
    return false;
  } else if (qkmath::fabs(y) == 1.0f) {
    if (x != 0.0f || z != 0.0f) {
      x = z = 0.0f;
      return true;
    }
    return false;
  } else if (qkmath::fabs(z) == 1.0f) {
    if (x != 0.0f || y != 0.0f) {
      x = y = 0.0f;
      return true;
    }
    return false;
  }
  return false;
}

bool idVec3::FixDenormals(void) {
  bool denormal = false;
  if (qkmath::fabs(x) < 1e-30f) {
    x = 0.0f;
    denormal = true;
  }
  if (qkmath::fabs(y) < 1e-30f) {
    y = 0.0f;
    denormal = true;
  }
  if (qkmath::fabs(z) < 1e-30f) {
    z = 0.0f;
    denormal = true;
  }
  return denormal;
}

idVec3 idVec3::Cross(const idVec3 &a) const {
  return idVec3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x);
}

idVec3 &idVec3::Cross(const idVec3 &a, const idVec3 &b) {
  x = a.y * b.z - a.z * b.y;
  y = a.z * b.x - a.x * b.z;
  z = a.x * b.y - a.y * b.x;

  return *this;
}

float idVec3::Length(void) const {
  return (float)sqrtf(x * x + y * y + z * z);
}

float idVec3::LengthSqr(void) const {
  return (x * x + y * y + z * z);
}

float idVec3::LengthFast(void) const {
  float sqrLength;

  sqrLength = x * x + y * y + z * z;
  return sqrLength * qkmath::fastinvsqrt(sqrLength);
}

float idVec3::Normalize(void) {
  float sqrLength, invLength;

  sqrLength = x * x + y * y + z * z;
  invLength = qkmath::inv_sqrt(sqrLength);
  x *= invLength;
  y *= invLength;
  z *= invLength;
  return invLength * sqrLength;
}

idVec3 &idVec3::Truncate(float length) {
  float length2;
  float ilength;

  if (!length) {
    Zero();
  }
  else {
    length2 = LengthSqr();
    if (length2 > length * length) {
      ilength = length * qkmath::inv_sqrt(length2);
      x *= ilength;
      y *= ilength;
      z *= ilength;
    }
  }

  return *this;
}

void idVec3::Clamp(const idVec3 &min, const idVec3 &max) {
  if (x < min.x) {
    x = min.x;
  } else if (x > max.x) {
    x = max.x;
  }
  if (y < min.y) {
    y = min.y;
  } else if (y > max.y) {
    y = max.y;
  }
  if (z < min.z) {
    z = min.z;
  } else if (z > max.z) {
    z = max.z;
  }
}

void idVec3::Snap(void) {
  x = floor(x + 0.5f);
  y = floor(y + 0.5f);
  z = floor(z + 0.5f);
}

void idVec3::SnapInt(void) {
  x = float(int(x));
  y = float(int(y));
  z = float(int(z));
}

int idVec3::GetDimension(void) const {
  return 3;
}

const idVec2 &idVec3::ToVec2(void) const {
  return *reinterpret_cast<const idVec2 *>(this);
}

idVec2 &idVec3::ToVec2(void) {
  return *reinterpret_cast<idVec2 *>(this);
}

const float *idVec3::ToFloatPtr(void) const {
  return &x;
}

float *idVec3::ToFloatPtr(void) {
  return &x;
}

void idVec3::NormalVectors(idVec3 &left, idVec3 &down) const {
  float d;

  d = x * x + y * y;
  if (!d) {
    left[0] = 1;
    left[1] = 0;
    left[2] = 0;
  } else {
    d = qkmath::inv_sqrt(d);
    left[0] = -y * d;
    left[1] = x * d;
    left[2] = 0;
  }
  down = left.Cross(*this);
}

void idVec3::OrthogonalBasis(idVec3 &left, idVec3 &up) const {
  float l, s;

  if (qkmath::fabs(z) > 0.7f) {
    l = y * y + z * z;
    s = qkmath::inv_sqrt(l);
    up[0] = 0;
    up[1] = z * s;
    up[2] = -y * s;
    left[0] = l * s;
    left[1] = -x * up[2];
    left[2] = x * up[1];
  }
  else {
    l = x * x + y * y;
    s = qkmath::inv_sqrt(l);
    left[0] = -y * s;
    left[1] = x * s;
    left[2] = 0;
    up[0] = -z * left[1];
    up[1] = z * left[0];
    up[2] = l * s;
  }
}

void idVec3::ProjectOntoPlane(const idVec3 &normal, const float overBounce) {
  float backoff;

  backoff = *this * normal;

  if (overBounce != 1.0) {
    if (backoff < 0) {
      backoff *= overBounce;
    } else {
      backoff /= overBounce;
    }
  }

  *this -= backoff * normal;
}

bool idVec3::ProjectAlongPlane(const idVec3 &normal, const float epsilon, const float overBounce) {
  idVec3 cross;
  float len;

  cross = this->Cross(normal).Cross((*this));
  // normalize so a fixed epsilon can be used
  cross.Normalize();
  len = normal * cross;
  if (qkmath::fabs(len) < epsilon) {
    return false;
  }
  cross *= overBounce * (normal * (*this)) / len;
  (*this) -= cross;
  return true;
}

idVec4::idVec4(void) {
}

idVec4::idVec4(const float x, const float y, const float z, const float w) {
  this->x = x;
  this->y = y;
  this->z = z;
  this->w = w;
}

void idVec4::Set(const float x, const float y, const float z, const float w) {
  this->x = x;
  this->y = y;
  this->z = z;
  this->w = w;
}

void idVec4::Zero(void) {
  x = y = z = w = 0.0f;
}

float idVec4::operator[](int index) const {
  return (&x)[ index ];
}

float& idVec4::operator[](int index) {
  return (&x)[ index ];
}

idVec4 idVec4::operator-() const {
  return idVec4(-x, -y, -z, -w);
}

idVec4 idVec4::operator-(const idVec4 &a) const {
  return idVec4(x - a.x, y - a.y, z - a.z, w - a.w);
}

float idVec4::operator*(const idVec4 &a) const {
  return x * a.x + y * a.y + z * a.z + w * a.w;
}

idVec4 idVec4::operator*(const float a) const {
  return idVec4(x * a, y * a, z * a, w * a);
}

idVec4 idVec4::operator/(const float a) const {
  float inva = 1.0f / a;
  return idVec4(x * inva, y * inva, z * inva, w * inva);
}

idVec4 operator*(const float a, const idVec4 b) {
  return idVec4(b.x * a, b.y * a, b.z * a, b.w * a);
}

idVec4 idVec4::operator+(const idVec4 &a) const {
  return idVec4(x + a.x, y + a.y, z + a.z, w + a.w);
}

idVec4 &idVec4::operator+=(const idVec4 &a) {
  x += a.x;
  y += a.y;
  z += a.z;
  w += a.w;

  return *this;
}

idVec4 &idVec4::operator/=(const idVec4 &a) {
  x /= a.x;
  y /= a.y;
  z /= a.z;
  w /= a.w;

  return *this;
}

idVec4 &idVec4::operator/=(const float a) {
  float inva = 1.0f / a;
  x *= inva;
  y *= inva;
  z *= inva;
  w *= inva;

  return *this;
}

idVec4 &idVec4::operator-=(const idVec4 &a) {
  x -= a.x;
  y -= a.y;
  z -= a.z;
  w -= a.w;

  return *this;
}

idVec4 &idVec4::operator*=(const float a) {
  x *= a;
  y *= a;
  z *= a;
  w *= a;

  return *this;
}

bool idVec4::Compare(const idVec4 &a) const {
  return ((x == a.x) && (y == a.y) && (z == a.z) && w == a.w);
}

bool idVec4::Compare(const idVec4 &a, const float epsilon) const {
  if (qkmath::fabs(x - a.x) > epsilon) {
    return false;
  }

  if (qkmath::fabs(y - a.y) > epsilon) {
    return false;
  }

  if (qkmath::fabs(z - a.z) > epsilon) {
    return false;
  }

  if (qkmath::fabs(w - a.w) > epsilon) {
    return false;
  }

  return true;
}

bool idVec4::operator==(const idVec4 &a) const {
  return Compare(a);
}

bool idVec4::operator!=(const idVec4 &a) const {
  return !Compare(a);
}

float idVec4::Length(void) const {
  return (float)sqrtf(x * x + y * y + z * z + w * w);
}

float idVec4::LengthSqr(void) const {
  return (x * x + y * y + z * z + w * w);
}

float idVec4::Normalize(void) {
  float sqrLength, invLength;

  sqrLength = x * x + y * y + z * z + w * w;
  invLength = qkmath::inv_sqrt(sqrLength);
  x *= invLength;
  y *= invLength;
  z *= invLength;
  w *= invLength;
  return invLength * sqrLength;
}

float idVec4::NormalizeFast(void) {
  float sqrLength, invLength;

  sqrLength = x * x + y * y + z * z + w * w;
  invLength = qkmath::fastinvsqrt(sqrLength);
  x *= invLength;
  y *= invLength;
  z *= invLength;
  w *= invLength;
  return invLength * sqrLength;
}

int idVec4::GetDimension(void) const {
  return 4;
}

const idVec2 &idVec4::ToVec2(void) const {
  return *reinterpret_cast<const idVec2 *>(this);
}

idVec2 &idVec4::ToVec2(void) {
  return *reinterpret_cast<idVec2 *>(this);
}

const idVec3 &idVec4::ToVec3(void) const {
  return *reinterpret_cast<const idVec3 *>(this);
}

idVec3 &idVec4::ToVec3(void) {
  return *reinterpret_cast<idVec3 *>(this);
}

const float *idVec4::ToFloatPtr(void) const {
  return &x;
}

float *idVec4::ToFloatPtr(void) {
  return &x;
}


idVec5::idVec5(void) {
}

idVec5::idVec5(const idVec3 &xyz, const idVec2 &st) {
  x = xyz.x;
  y = xyz.y;
  z = xyz.z;
  s = st[0];
  t = st[1];
}

idVec5::idVec5(const float x, const float y, const float z, const float s, const float t) {
  this->x = x;
  this->y = y;
  this->z = z;
  this->s = s;
  this->t = t;
}

float idVec5::operator[](int index) const {
  return (&x)[ index ];
}

float& idVec5::operator[](int index) {
  return (&x)[ index ];
}

idVec5 &idVec5::operator=(const idVec3 &a) {
  x = a.x;
  y = a.y;
  z = a.z;
  s = t = 0;
  return *this;
}

int idVec5::GetDimension(void) const {
  return 5;
}

const idVec3 &idVec5::ToVec3(void) const {
  return *reinterpret_cast<const idVec3 *>(this);
}

idVec3 &idVec5::ToVec3(void) {
  return *reinterpret_cast<idVec3 *>(this);
}

const float *idVec5::ToFloatPtr(void) const {
  return &x;
}

float *idVec5::ToFloatPtr(void) {
  return &x;
}


idVec6::idVec6(void) {
}

idVec6::idVec6(const float *a) {
  memcpy(p, a, 6 * sizeof(float));
}

idVec6::idVec6(const float a1, const float a2, const float a3, const float a4, const float a5, const float a6) {
  p[0] = a1;
  p[1] = a2;
  p[2] = a3;
  p[3] = a4;
  p[4] = a5;
  p[5] = a6;
}

idVec6 idVec6::operator-() const {
  return idVec6(-p[0], -p[1], -p[2], -p[3], -p[4], -p[5]);
}

float idVec6::operator[](const int index) const {
  return p[index];
}

float &idVec6::operator[](const int index) {
  return p[index];
}

idVec6 idVec6::operator*(const float a) const {
  return idVec6(p[0]*a, p[1]*a, p[2]*a, p[3]*a, p[4]*a, p[5]*a);
}

float idVec6::operator*(const idVec6 &a) const {
  return p[0] * a[0] + p[1] * a[1] + p[2] * a[2] + p[3] * a[3] + p[4] * a[4] + p[5] * a[5];
}

idVec6 idVec6::operator/(const float a) const {
  float inva;

  assert(a != 0.0f);
  inva = 1.0f / a;
  return idVec6(p[0]*inva, p[1]*inva, p[2]*inva, p[3]*inva, p[4]*inva, p[5]*inva);
}

idVec6 idVec6::operator+(const idVec6 &a) const {
  return idVec6(p[0] + a[0], p[1] + a[1], p[2] + a[2], p[3] + a[3], p[4] + a[4], p[5] + a[5]);
}

idVec6 idVec6::operator-(const idVec6 &a) const {
  return idVec6(p[0] - a[0], p[1] - a[1], p[2] - a[2], p[3] - a[3], p[4] - a[4], p[5] - a[5]);
}

idVec6 &idVec6::operator*=(const float a) {
  p[0] *= a;
  p[1] *= a;
  p[2] *= a;
  p[3] *= a;
  p[4] *= a;
  p[5] *= a;
  return *this;
}

idVec6 &idVec6::operator/=(const float a) {
  float inva;

  assert(a != 0.0f);
  inva = 1.0f / a;
  p[0] *= inva;
  p[1] *= inva;
  p[2] *= inva;
  p[3] *= inva;
  p[4] *= inva;
  p[5] *= inva;
  return *this;
}

idVec6 &idVec6::operator+=(const idVec6 &a) {
  p[0] += a[0];
  p[1] += a[1];
  p[2] += a[2];
  p[3] += a[3];
  p[4] += a[4];
  p[5] += a[5];
  return *this;
}

idVec6 &idVec6::operator-=(const idVec6 &a) {
  p[0] -= a[0];
  p[1] -= a[1];
  p[2] -= a[2];
  p[3] -= a[3];
  p[4] -= a[4];
  p[5] -= a[5];
  return *this;
}

idVec6 operator*(const float a, const idVec6 b) {
  return b * a;
}

bool idVec6::Compare(const idVec6 &a) const {
  return ((p[0] == a[0]) && (p[1] == a[1]) && (p[2] == a[2]) &&
      (p[3] == a[3]) && (p[4] == a[4]) && (p[5] == a[5]));
}

bool idVec6::Compare(const idVec6 &a, const float epsilon) const {
  if (qkmath::fabs(p[0] - a[0]) > epsilon) {
    return false;
  }

  if (qkmath::fabs(p[1] - a[1]) > epsilon) {
    return false;
  }

  if (qkmath::fabs(p[2] - a[2]) > epsilon) {
    return false;
  }

  if (qkmath::fabs(p[3] - a[3]) > epsilon) {
    return false;
  }

  if (qkmath::fabs(p[4] - a[4]) > epsilon) {
    return false;
  }

  if (qkmath::fabs(p[5] - a[5]) > epsilon) {
    return false;
  }

  return true;
}

bool idVec6::operator==(const idVec6 &a) const {
  return Compare(a);
}

bool idVec6::operator!=(const idVec6 &a) const {
  return !Compare(a);
}

void idVec6::Set(const float a1, const float a2, const float a3, const float a4, const float a5, const float a6) {
  p[0] = a1;
  p[1] = a2;
  p[2] = a3;
  p[3] = a4;
  p[4] = a5;
  p[5] = a6;
}

void idVec6::Zero(void) {
  p[0] = p[1] = p[2] = p[3] = p[4] = p[5] = 0.0f;
}

float idVec6::Length(void) const {
  return (float)sqrtf(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + p[4] * p[4] + p[5] * p[5]);
}

float idVec6::LengthSqr(void) const {
  return (p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + p[4] * p[4] + p[5] * p[5]);
}

float idVec6::Normalize(void) {
  float sqrLength, invLength;

  sqrLength = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + p[4] * p[4] + p[5] * p[5];
  invLength = qkmath::inv_sqrt(sqrLength);
  p[0] *= invLength;
  p[1] *= invLength;
  p[2] *= invLength;
  p[3] *= invLength;
  p[4] *= invLength;
  p[5] *= invLength;
  return invLength * sqrLength;
}

float idVec6::NormalizeFast(void) {
  float sqrLength, invLength;

  sqrLength = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + p[4] * p[4] + p[5] * p[5];
  invLength = qkmath::fastinvsqrt(sqrLength);
  p[0] *= invLength;
  p[1] *= invLength;
  p[2] *= invLength;
  p[3] *= invLength;
  p[4] *= invLength;
  p[5] *= invLength;
  return invLength * sqrLength;
}

int idVec6::GetDimension(void) const {
  return 6;
}

const idVec3 &idVec6::SubVec3(int index) const {
  return *reinterpret_cast<const idVec3 *>(p + index * 3);
}

idVec3 &idVec6::SubVec3(int index) {
  return *reinterpret_cast<idVec3 *>(p + index * 3);
}

const float *idVec6::ToFloatPtr(void) const {
  return p;
}

float *idVec6::ToFloatPtr(void) {
  return p;
}


idVecX::idVecX(void) {
  size = alloced = 0;
  p = NULL;
}

idVecX::idVecX(int length) {
  size = alloced = 0;
  p = NULL;
  SetSize(length);
}

idVecX::idVecX(int length, float *data) {
  size = alloced = 0;
  p = NULL;
  SetData(length, data);
}

idVecX::~idVecX(void) {
  // if not temp memory
  if (p && (p < idVecX::tempPtr || p >= idVecX::tempPtr + VECX_MAX_TEMP) && alloced != -1) {
    free(p);
  }
}

float idVecX::operator[](const int index) const {
  assert(index >= 0 && index < size);
  return p[index];
}

float &idVecX::operator[](const int index) {
  assert(index >= 0 && index < size);
  return p[index];
}

idVecX idVecX::operator-() const {
  int i;
  idVecX m;

  m.SetTempSize(size);
  for (i = 0; i < size; i++) {
    m.p[i] = -p[i];
  }
  return m;
}

idVecX &idVecX::operator=(const idVecX &a) {
  SetSize(a.size);
  memcpy(p, a.p, a.size * sizeof(float));
  idVecX::tempIndex = 0;
  return *this;
}

idVecX idVecX::operator+(const idVecX &a) const {
  idVecX m;

  assert(size == a.size);
  m.SetTempSize(size);
  int i;
  for (i = 0; i < size; i++) {
    m.p[i] = p[i] + a.p[i];
  }
  return m;
}

idVecX idVecX::operator-(const idVecX &a) const {
  idVecX m;

  assert(size == a.size);
  m.SetTempSize(size);
  int i;
  for (i = 0; i < size; i++) {
    m.p[i] = p[i] - a.p[i];
  }
  return m;
}

idVecX &idVecX::operator+=(const idVecX &a) {
  assert(size == a.size);
  int i;
  for (i = 0; i < size; i++) {
    p[i] += a.p[i];
  }
  idVecX::tempIndex = 0;
  return *this;
}

idVecX &idVecX::operator-=(const idVecX &a) {
  assert(size == a.size);
  int i;
  for (i = 0; i < size; i++) {
    p[i] -= a.p[i];
  }
  idVecX::tempIndex = 0;
  return *this;
}

idVecX idVecX::operator*(const float a) const {
  idVecX m;

  m.SetTempSize(size);
  int i;
  for (i = 0; i < size; i++) {
    m.p[i] = p[i] * a;
  }
  return m;
}

idVecX &idVecX::operator*=(const float a) {
  int i;
  for (i = 0; i < size; i++) {
    p[i] *= a;
  }
  return *this;
}

idVecX idVecX::operator/(const float a) const {
  assert(a != 0.0f);
  return (*this) * (1.0f / a);
}

idVecX &idVecX::operator/=(const float a) {
  assert(a != 0.0f);
  (*this) *= (1.0f / a);
  return *this;
}

idVecX operator*(const float a, const idVecX b) {
  return b * a;
}

float idVecX::operator*(const idVecX &a) const {
  int i;
  float sum = 0.0f;

  assert(size == a.size);
  for (i = 0; i < size; i++) {
    sum += p[i] * a.p[i];
  }
  return sum;
}

bool idVecX::Compare(const idVecX &a) const {
  int i;

  assert(size == a.size);
  for (i = 0; i < size; i++) {
    if (p[i] != a.p[i]) {
      return false;
    }
  }
  return true;
}

bool idVecX::Compare(const idVecX &a, const float epsilon) const {
  int i;

  assert(size == a.size);
  for (i = 0; i < size; i++) {
    if (qkmath::fabs(p[i] - a.p[i]) > epsilon) {
      return false;
    }
  }
  return true;
}

bool idVecX::operator==(const idVecX &a) const {
  return Compare(a);
}

bool idVecX::operator!=(const idVecX &a) const {
  return !Compare(a);
}

void idVecX::SetSize(int newSize) {
  int alloc = (newSize + 3) & ~3;
  if (alloc > alloced && alloced != -1) {
    if (p) {
      free(p);
    }
    p = (float *) malloc(alloc * sizeof(float));
    alloced = alloc;
  }
  size = newSize;
  VECX_CLEAREND();
}

void idVecX::ChangeSize(int newSize, bool makeZero) {
  int alloc = (newSize + 3) & ~3;
  if (alloc > alloced && alloced != -1) {
    float *oldVec = p;
    p = (float *) malloc(alloc * sizeof(float));
    alloced = alloc;
    if (oldVec) {
      for (int i = 0; i < size; i++) {
        p[i] = oldVec[i];
      }
      free(oldVec);
    }
    if (makeZero) {
      // zero any new elements
      for (int i = size; i < newSize; i++) {
        p[i] = 0.0f;
      }
    }
  }
  size = newSize;
  VECX_CLEAREND();
}

void idVecX::SetTempSize(int newSize) {

  size = newSize;
  alloced = (newSize + 3) & ~3;
  assert(alloced < VECX_MAX_TEMP);
  if (idVecX::tempIndex + alloced > VECX_MAX_TEMP) {
    idVecX::tempIndex = 0;
  }
  p = idVecX::tempPtr + idVecX::tempIndex;
  idVecX::tempIndex += alloced;
  VECX_CLEAREND();
}

void idVecX::SetData(int length, float *data) {
  if (p && (p < idVecX::tempPtr || p >= idVecX::tempPtr + VECX_MAX_TEMP) && alloced != -1) {
    free(p);
  }
  assert((((uintptr_t) data) & 15) == 0); // data must be 16 byte aligned
  p = data;
  size = length;
  alloced = -1;
  VECX_CLEAREND();
}

void idVecX::Zero(void) {
  memset(p, 0, size * sizeof(float));
}

void idVecX::Zero(int length) {
  SetSize(length);
  memset(p, 0, size * sizeof(float));
}

void idVecX::Negate(void) {
  int i;
  for (i = 0; i < size; i++) {
    p[i] = -p[i];
  }
}

void idVecX::Clamp(float min, float max) {
  int i;
  for (i = 0; i < size; i++) {
    if (p[i] < min) {
      p[i] = min;
    } else if (p[i] > max) {
      p[i] = max;
    }
  }
}

idVecX &idVecX::SwapElements(int e1, int e2) {
  float tmp;
  tmp = p[e1];
  p[e1] = p[e2];
  p[e2] = tmp;
  return *this;
}

float idVecX::Length(void) const {
  int i;
  float sum = 0.0f;

  for (i = 0; i < size; i++) {
    sum += p[i] * p[i];
  }
  return sqrtf(sum);
}

float idVecX::LengthSqr(void) const {
  int i;
  float sum = 0.0f;

  for (i = 0; i < size; i++) {
    sum += p[i] * p[i];
  }
  return sum;
}

idVecX idVecX::Normalize(void) const {
  int i;
  idVecX m;
  float invSqrt, sum = 0.0f;

  m.SetTempSize(size);
  for (i = 0; i < size; i++) {
    sum += p[i] * p[i];
  }
  invSqrt = qkmath::inv_sqrt(sum);
  for (i = 0; i < size; i++) {
    m.p[i] = p[i] * invSqrt;
  }
  return m;
}

float idVecX::NormalizeSelf(void) {
  float invSqrt, sum = 0.0f;
  int i;
  for (i = 0; i < size; i++) {
    sum += p[i] * p[i];
  }
  invSqrt = qkmath::inv_sqrt(sum);
  for (i = 0; i < size; i++) {
    p[i] *= invSqrt;
  }
  return invSqrt * sum;
}

int idVecX::GetDimension(void) const {
  return size;
}

idVec3 &idVecX::SubVec3(int index) {
  assert(index >= 0 && index * 3 + 3 <= size);
  return *reinterpret_cast<idVec3 *>(p + index * 3);
}

const idVec3 &idVecX::SubVec3(int index) const {
  assert(index >= 0 && index * 3 + 3 <= size);
  return *reinterpret_cast<const idVec3 *>(p + index * 3);
}

idVec6 &idVecX::SubVec6(int index) {
  assert(index >= 0 && index * 6 + 6 <= size);
  return *reinterpret_cast<idVec6 *>(p + index * 6);
}

const idVec6 &idVecX::SubVec6(int index) const {
  assert(index >= 0 && index * 6 + 6 <= size);
  return *reinterpret_cast<const idVec6 *>(p + index * 6);
}

const float *idVecX::ToFloatPtr(void) const {
  return p;
}

float *idVecX::ToFloatPtr(void) {
  return p;
}


//===============================================================
//
//	idPolar3
//
//===============================================================

class idPolar3 {
public:
  float			radius, theta, phi;

  idPolar3(void);
  explicit idPolar3(const float radius, const float theta, const float phi);

  void			Set(const float radius, const float theta, const float phi);

  float			operator[](const int index) const;
  float &			operator[](const int index);
  idPolar3		operator-() const;
  idPolar3 &		operator=(const idPolar3 &a);

  idVec3			ToVec3(void) const;
};

idPolar3::idPolar3(void) {
}

idPolar3::idPolar3(const float radius, const float theta, const float phi) {
  assert(radius > 0);
  this->radius = radius;
  this->theta = theta;
  this->phi = phi;
}

void idPolar3::Set(const float radius, const float theta, const float phi) {
  assert(radius > 0);
  this->radius = radius;
  this->theta = theta;
  this->phi = phi;
}

float idPolar3::operator[](const int index) const {
  return (&radius)[ index ];
}

float &idPolar3::operator[](const int index) {
  return (&radius)[ index ];
}

idPolar3 idPolar3::operator-() const {
  return idPolar3(radius, -theta, -phi);
}

idPolar3 &idPolar3::operator=(const idPolar3 &a) {
  radius = a.radius;
  theta = a.theta;
  phi = a.phi;
  return *this;
}

idVec3 idPolar3::ToVec3(void) const {
  float sp, cp, st, ct;
  sp = sinf(phi);
  cp = cosf(phi);
  st = sinf(theta);
  ct = cosf(theta);
  return idVec3(cp * radius * ct, cp * radius * st, radius * sp);
}


idVec2 vec2_origin( 0.0f, 0.0f );
idVec3 vec3_origin( 0.0f, 0.0f, 0.0f );
idVec4 vec4_origin( 0.0f, 0.0f, 0.0f, 0.0f );
idVec5 vec5_origin( 0.0f, 0.0f, 0.0f, 0.0f, 0.0f );
idVec6 vec6_origin( 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f );
idVec6 vec6_infinity( qkmath::C_INFINITY, qkmath::C_INFINITY, qkmath::C_INFINITY, qkmath::C_INFINITY, qkmath::C_INFINITY, qkmath::C_INFINITY );

//===============================================================
//
//	idVec2
//
//===============================================================


/*
=============
Lerp

Linearly inperpolates one vector to another.
=============
*/
void idVec2::Lerp( const idVec2 &v1, const idVec2 &v2, const float l ) {
	if ( l <= 0.0f ) {
		(*this) = v1;
	} else if ( l >= 1.0f ) {
		(*this) = v2;
	} else {
		(*this) = v1 + l * ( v2 - v1 );
	}
}


//===============================================================
//
//	idVec3
//
//===============================================================

/*
=============
idVec3::ToYaw
=============
*/
float idVec3::ToYaw( void ) const {
	float yaw;

	if ( ( y == 0.0f ) && ( x == 0.0f ) ) {
		yaw = 0.0f;
	} else {
		yaw = RAD2DEG( atan2( y, x ) );
		if ( yaw < 0.0f ) {
			yaw += 360.0f;
		}
	}

	return yaw;
}

/*
=============
idVec3::ToPitch
=============
*/
float idVec3::ToPitch( void ) const {
	float	forward;
	float	pitch;

	if ( ( x == 0.0f ) && ( y == 0.0f ) ) {
		if ( z > 0.0f ) {
			pitch = 90.0f;
		} else {
			pitch = 270.0f;
		}
	} else {
		forward = ( float )qkmath::sqrt( x * x + y * y );
		pitch = RAD2DEG( atan2( z, forward ) );
		if ( pitch < 0.0f ) {
			pitch += 360.0f;
		}
	}

	return pitch;
}

#if 0
/*
=============
idVec3::ToAngles
=============
*/
idAngles idVec3::ToAngles( void ) const {
	float forward;
	float yaw;
	float pitch;

	if ( ( x == 0.0f ) && ( y == 0.0f ) ) {
		yaw = 0.0f;
		if ( z > 0.0f ) {
			pitch = 90.0f;
		} else {
			pitch = 270.0f;
		}
	} else {
		yaw = RAD2DEG( atan2( y, x ) );
		if ( yaw < 0.0f ) {
			yaw += 360.0f;
		}

		forward = ( float )qkmath::sqrt( x * x + y * y );
		pitch = RAD2DEG( atan2( z, forward ) );
		if ( pitch < 0.0f ) {
			pitch += 360.0f;
		}
	}

	return idAngles( -pitch, yaw, 0.0f );
}
#endif

/*
=============
idVec3::ToPolar
=============
*/
idPolar3 idVec3::ToPolar( void ) const {
	float forward;
	float yaw;
	float pitch;

	if ( ( x == 0.0f ) && ( y == 0.0f ) ) {
		yaw = 0.0f;
		if ( z > 0.0f ) {
			pitch = 90.0f;
		} else {
			pitch = 270.0f;
		}
	} else {
		yaw = RAD2DEG( atan2( y, x ) );
		if ( yaw < 0.0f ) {
			yaw += 360.0f;
		}

		forward = ( float )qkmath::sqrt( x * x + y * y );
		pitch = RAD2DEG( atan2( z, forward ) );
		if ( pitch < 0.0f ) {
			pitch += 360.0f;
		}
	}
	return idPolar3( qkmath::sqrt( x * x + y * y + z * z ), yaw, -pitch );
}

#if 0
/*
=============
idVec3::ToMat3
=============
*/
idMat3 idVec3::ToMat3( void ) const {
	idMat3	mat;
	float	d;

	mat[0] = *this;
	d = x * x + y * y;
	if ( !d ) {
		mat[1][0] = 1.0f;
		mat[1][1] = 0.0f;
		mat[1][2] = 0.0f;
	} else {
		d = qkmath::inv_sqrt( d );
		mat[1][0] = -y * d;
		mat[1][1] = x * d;
		mat[1][2] = 0.0f;
	}
	mat[2] = Cross( mat[1] );

	return mat;
}
#endif

/*
=============
Lerp

Linearly inperpolates one vector to another.
=============
*/
void idVec3::Lerp( const idVec3 &v1, const idVec3 &v2, const float l ) {
	if ( l <= 0.0f ) {
		(*this) = v1;
	} else if ( l >= 1.0f ) {
		(*this) = v2;
	} else {
		(*this) = v1 + l * ( v2 - v1 );
	}
}

/*
=============
SLerp

Spherical linear interpolation from v1 to v2.
Vectors are expected to be normalized.
=============
*/
#define LERP_DELTA 1e-6

void idVec3::SLerp( const idVec3 &v1, const idVec3 &v2, const float t ) {
	float omega, cosom, sinom, scale0, scale1;

	if ( t <= 0.0f ) {
		(*this) = v1;
		return;
	} else if ( t >= 1.0f ) {
		(*this) = v2;
		return;
	}

	cosom = v1 * v2;
	if ( ( 1.0f - cosom ) > LERP_DELTA ) {
		omega = acos( cosom );
		sinom = sin( omega );
		scale0 = sin( ( 1.0f - t ) * omega ) / sinom;
		scale1 = sin( t * omega ) / sinom;
	} else {
		scale0 = 1.0f - t;
		scale1 = t;
	}

	(*this) = ( v1 * scale0 + v2 * scale1 );
}

/*
=============
ProjectSelfOntoSphere

Projects the z component onto a sphere.
=============
*/
void idVec3::ProjectSelfOntoSphere( const float radius ) {
	float rsqr = radius * radius;
	float len = Length();
	if ( len  < rsqr * 0.5f ) {
		z = sqrt( rsqr - len );
	} else {
		z = rsqr / ( 2.0f * sqrt( len ) );
	}
}



//===============================================================
//
//	idVec4
//
//===============================================================


/*
=============
Lerp

Linearly inperpolates one vector to another.
=============
*/
void idVec4::Lerp( const idVec4 &v1, const idVec4 &v2, const float l ) {
	if ( l <= 0.0f ) {
		(*this) = v1;
	} else if ( l >= 1.0f ) {
		(*this) = v2;
	} else {
		(*this) = v1 + l * ( v2 - v1 );
	}
}


//===============================================================
//
//	idVec5
//
//===============================================================

/*
=============
idVec5::Lerp
=============
*/
void idVec5::Lerp( const idVec5 &v1, const idVec5 &v2, const float l ) {
	if ( l <= 0.0f ) {
		(*this) = v1;
	} else if ( l >= 1.0f ) {
		(*this) = v2;
	} else {
		x = v1.x + l * ( v2.x - v1.x );
		y = v1.y + l * ( v2.y - v1.y );
		z = v1.z + l * ( v2.z - v1.z );
		s = v1.s + l * ( v2.s - v1.s );
		t = v1.t + l * ( v2.t - v1.t );
	}
}


//===============================================================
//
//	idVec6
//
//===============================================================


//===============================================================
//
//	idVecX
//
//===============================================================

float	idVecX::temp[VECX_MAX_TEMP+4];
float *	idVecX::tempPtr = (float *) ( ( (intptr_t) idVecX::temp + 15 ) & ~15 );
int		idVecX::tempIndex = 0;

