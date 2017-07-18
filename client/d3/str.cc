#include "str.hh"

void idStr::EnsureAlloced(int amount, bool keepold) {
  if (amount > alloced) {
    ReAllocate(amount, keepold);
  }
}

void idStr::Init(void) {
  len = 0;
  alloced = STR_ALLOC_BASE;
  data = baseBuffer;
  data[ 0 ] = '\0';
#ifdef ID_DEBUG_UNINITIALIZED_MEMORY
  memset(baseBuffer, 0, sizeof(baseBuffer));
#endif
}

idStr::idStr(void) {
  Init();
}

idStr::idStr(const idStr &text) {
  int l;

  Init();
  l = text.Length();
  EnsureAlloced(l + 1);
  strcpy(data, text.data);
  len = l;
}

idStr::idStr(const idStr &text, int start, int end) {
  int i;
  int l;

  Init();
  if (end > text.Length()) {
    end = text.Length();
  }
  if (start > text.Length()) {
    start = text.Length();
  } else if (start < 0) {
    start = 0;
  }

  l = end - start;
  if (l < 0) {
    l = 0;
  }

  EnsureAlloced(l + 1);

  for (i = 0; i < l; i++) {
    data[ i ] = text[ start + i ];
  }

  data[ l ] = '\0';
  len = l;
}

idStr::idStr(const char *text) {
  int l;

  Init();
  if (text) {
    l = strlen(text);
    EnsureAlloced(l + 1);
    strcpy(data, text);
    len = l;
  }
}

idStr::idStr(const char *text, int start, int end) {
  int i;
  int l = strlen(text);

  Init();
  if (end > l) {
    end = l;
  }
  if (start > l) {
    start = l;
  } else if (start < 0) {
    start = 0;
  }

  l = end - start;
  if (l < 0) {
    l = 0;
  }

  EnsureAlloced(l + 1);

  for (i = 0; i < l; i++) {
    data[ i ] = text[ start + i ];
  }

  data[ l ] = '\0';
  len = l;
}

idStr::idStr(const bool b) {
  Init();
  EnsureAlloced(2);
  data[ 0 ] = b ? '1' : '0';
  data[ 1 ] = '\0';
  len = 1;
}

idStr::idStr(const char c) {
  Init();
  EnsureAlloced(2);
  data[ 0 ] = c;
  data[ 1 ] = '\0';
  len = 1;
}

idStr::idStr(const int i) {
  char text[ 64 ];
  int l;

  Init();
  l = sprintf(text, "%d", i);
  EnsureAlloced(l + 1);
  strcpy(data, text);
  len = l;
}

idStr::idStr(const unsigned u) {
  char text[ 64 ];
  int l;

  Init();
  l = sprintf(text, "%u", u);
  EnsureAlloced(l + 1);
  strcpy(data, text);
  len = l;
}

idStr::idStr(const float f) {
  char text[ 64 ];
  int l;

  Init();
  l = idStr::snPrintf(text, sizeof(text), "%f", (double)f);
  while(l > 0 && text[l-1] == '0') text[--l] = '\0';
  while(l > 0 && text[l-1] == '.') text[--l] = '\0';
  EnsureAlloced(l + 1);
  strcpy(data, text);
  len = l;
}

idStr::~idStr(void) {
  FreeData();
}

size_t idStr::Size(void) const {
  return sizeof(*this) + Allocated();
}

const char *idStr::c_str(void) const {
  return data;
}

idStr::operator const char *(void) {
  return c_str();
}

idStr::operator const char *(void) const {
  return c_str();
}

#pragma GCC diagnostic push
// shut up GCC's stupid "warning: assuming signed overflow does not occur when assuming that
// (X - c) > X is always false [-Wstrict-overflow]"
#pragma GCC diagnostic ignored "-Wstrict-overflow"
char idStr::operator[](int index) const {
  assert((index >= 0) && (index <= len));
  return data[ index ];
}

char &idStr::operator[](int index) {
  assert((index >= 0) && (index <= len));
  return data[ index ];
}
#pragma GCC diagnostic pop

void idStr::operator=(const idStr &text) {
  if (&text == this) {
    return;
  }

  int l;

  l = text.Length();
  EnsureAlloced(l + 1, false);
  memcpy(data, text.data, l);
  data[l] = '\0';
  len = l;
}

idStr operator+(const idStr &a, const idStr &b) {
  idStr result(a);
  result.Append(b);
  return result;
}

idStr operator+(const idStr &a, const char *b) {
  idStr result(a);
  result.Append(b);
  return result;
}

idStr operator+(const char *a, const idStr &b) {
  idStr result(a);
  result.Append(b);
  return result;
}

idStr operator+(const idStr &a, const bool b) {
  idStr result(a);
  result.Append(b ? "true" : "false");
  return result;
}

idStr operator+(const idStr &a, const char b) {
  idStr result(a);
  result.Append(b);
  return result;
}

idStr operator+(const idStr &a, const float b) {
  char	text[ 64 ];
  idStr	result(a);

  sprintf(text, "%f", (double)b);
  result.Append(text);

  return result;
}

idStr operator+(const idStr &a, const int b) {
  char	text[ 64 ];
  idStr	result(a);

  sprintf(text, "%d", b);
  result.Append(text);

  return result;
}

idStr operator+(const idStr &a, const unsigned b) {
  char	text[ 64 ];
  idStr	result(a);

  sprintf(text, "%u", b);
  result.Append(text);

  return result;
}

idStr &idStr::operator+=(const float a) {
  char text[ 64 ];

  sprintf(text, "%f", (double)a);
  Append(text);

  return *this;
}

idStr &idStr::operator+=(const int a) {
  char text[ 64 ];

  sprintf(text, "%d", a);
  Append(text);

  return *this;
}

idStr &idStr::operator+=(const unsigned a) {
  char text[ 64 ];

  sprintf(text, "%u", a);
  Append(text);

  return *this;
}

idStr &idStr::operator+=(const idStr &a) {
  Append(a);
  return *this;
}

idStr &idStr::operator+=(const char *a) {
  Append(a);
  return *this;
}

idStr &idStr::operator+=(const char a) {
  Append(a);
  return *this;
}

idStr &idStr::operator+=(const bool a) {
  Append(a ? "true" : "false");
  return *this;
}

bool operator==(const idStr &a, const idStr &b) {
  return (!idStr::Cmp(a.data, b.data));
}

bool operator==(const idStr &a, const char *b) {
  assert(b);
  return (!idStr::Cmp(a.data, b));
}

bool operator==(const char *a, const idStr &b) {
  assert(a);
  return (!idStr::Cmp(a, b.data));
}

bool operator!=(const idStr &a, const idStr &b) {
  return !(a == b);
}

bool operator!=(const idStr &a, const char *b) {
  return !(a == b);
}

bool operator!=(const char *a, const idStr &b) {
  return !(a == b);
}

int idStr::Cmp(const char *text) const {
  assert(text);
  return idStr::Cmp(data, text);
}

int idStr::Cmpn(const char *text, int n) const {
  assert(text);
  return idStr::Cmpn(data, text, n);
}

int idStr::CmpPrefix(const char *text) const {
  assert(text);
  return idStr::Cmpn(data, text, strlen(text));
}

int idStr::Icmp(const char *text) const {
  assert(text);
  return idStr::Icmp(data, text);
}

int idStr::Icmpn(const char *text, int n) const {
  assert(text);
  return idStr::Icmpn(data, text, n);
}

int idStr::IcmpPrefix(const char *text) const {
  assert(text);
  return idStr::Icmpn(data, text, strlen(text));
}

int idStr::IcmpPath(const char *text) const {
  assert(text);
  return idStr::IcmpPath(data, text);
}

int idStr::IcmpnPath(const char *text, int n) const {
  assert(text);
  return idStr::IcmpnPath(data, text, n);
}

int idStr::IcmpPrefixPath(const char *text) const {
  assert(text);
  return idStr::IcmpnPath(data, text, strlen(text));
}

int idStr::Length(void) const {
  return len;
}

int idStr::Allocated(void) const {
  if (data != baseBuffer) {
    return alloced;
  } else {
    return 0;
  }
}

void idStr::Empty(void) {
  EnsureAlloced(1);
  data[ 0 ] = '\0';
  len = 0;
}

bool idStr::IsEmpty(void) const {
  return (idStr::Cmp(data, "") == 0);
}

void idStr::Clear(void) {
  FreeData();
  Init();
}

void idStr::Append(const char a) {
  EnsureAlloced(len + 2);
  data[ len ] = a;
  len++;
  data[ len ] = '\0';
}

void idStr::Append(const idStr &text) {
  int newLen;
  int i;

  newLen = len + text.Length();
  EnsureAlloced(newLen + 1);
  for (i = 0; i < text.len; i++) {
    data[ len + i ] = text[ i ];
  }
  len = newLen;
  data[ len ] = '\0';
}

void idStr::Append(const char *text) {
  int newLen;
  int i;

  if (text) {
    newLen = len + strlen(text);
    EnsureAlloced(newLen + 1);
    for (i = 0; text[ i ]; i++) {
      data[ len + i ] = text[ i ];
    }
    len = newLen;
    data[ len ] = '\0';
  }
}

void idStr::Append(const char *text, int l) {
  int newLen;
  int i;

  if (text && l) {
    newLen = len + l;
    EnsureAlloced(newLen + 1);
    for (i = 0; text[ i ] && i < l; i++) {
      data[ len + i ] = text[ i ];
    }
    len = newLen;
    data[ len ] = '\0';
  }
}

void idStr::Insert(const char a, int index) {
  int i, l;

  if (index < 0) {
    index = 0;
  } else if (index > len) {
    index = len;
  }

  l = 1;
  EnsureAlloced(len + l + 1);
  for (i = len; i >= index; i--) {
    data[i+l] = data[i];
  }
  data[index] = a;
  len++;
}

void idStr::Insert(const char *text, int index) {
  int i, l;

  if (index < 0) {
    index = 0;
  } else if (index > len) {
    index = len;
  }

  l = strlen(text);
  EnsureAlloced(len + l + 1);
  for (i = len; i >= index; i--) {
    data[i+l] = data[i];
  }
  for (i = 0; i < l; i++) {
    data[index+i] = text[i];
  }
  len += l;
}

void idStr::ToLower(void) {
  for (int i = 0; data[i]; i++ ) {
    if (CharIsUpper(data[i])) {
      data[i] += ('a' - 'A');
    }
  }
}

void idStr::ToUpper(void) {
  for (int i = 0; data[i]; i++ ) {
    if (CharIsLower(data[i])) {
      data[i] -= ('a' - 'A');
    }
  }
}

bool idStr::IsNumeric(void) const {
  return idStr::IsNumeric(data);
}

bool idStr::HasLower(void) const {
  return idStr::HasLower(data);
}

bool idStr::HasUpper(void) const {
  return idStr::HasUpper(data);
}

void idStr::CapLength(int newlen) {
  if (len <= newlen) {
    return;
  }
  data[ newlen ] = 0;
  len = newlen;
}

void idStr::Fill(const char ch, int newlen) {
  EnsureAlloced(newlen + 1);
  len = newlen;
  memset(data, ch, len);
  data[ len ] = 0;
}

int idStr::Find(const char c, int start, int end) const {
  if (end == -1) {
    end = len;
  }
  return idStr::FindChar(data, c, start, end);
}

int idStr::Find(const char *text, bool casesensitive, int start, int end) const {
  if (end == -1) {
    end = len;
  }
  return idStr::FindText(data, text, casesensitive, start, end);
}

bool idStr::Filter(const char *filter, bool casesensitive) const {
  return idStr::Filter(filter, data, casesensitive);
}

const char *idStr::Left(int len, idStr &result) const {
  return Mid(0, len, result);
}

const char *idStr::Right(int len, idStr &result) const {
  if (len >= Length()) {
    result = *this;
    return result;
  }
  return Mid(Length() - len, len, result);
}

idStr idStr::Left(int len) const {
  return Mid(0, len);
}

#pragma GCC diagnostic push
// shut up GCC's stupid "warning: assuming signed overflow does not occur when assuming that
// (X - c) > X is always false [-Wstrict-overflow]"
#pragma GCC diagnostic ignored "-Wstrict-overflow"
idStr idStr::Right(int len) const {
  if (len >= Length()) {
    return *this;
  }
  return Mid(Length() - len, len);
}
#pragma GCC diagnostic pop

void idStr::Strip(const char c) {
  StripLeading(c);
  StripTrailing(c);
}

void idStr::Strip(const char *string) {
  StripLeading(string);
  StripTrailing(string);
}

bool idStr::CheckExtension(const char *ext) {
  return idStr::CheckExtension(data, ext);
}

int idStr::Length(const char *s) {
  int i;
  for (i = 0; s[i]; i++) {}
  return i;
}

char *idStr::ToLower(char *s) {
  for (int i = 0; s[i]; i++) {
    if (CharIsUpper(s[i])) {
      s[i] += ('a' - 'A');
    }
  }
  return s;
}

char *idStr::ToUpper(char *s) {
  for (int i = 0; s[i]; i++) {
    if (CharIsLower(s[i])) {
      s[i] -= ('a' - 'A');
    }
  }
  return s;
}

int idStr::Hash(const char *string) {
  int i, hash = 0;
  for (i = 0; *string != '\0'; i++) {
    hash += (*string++) * (i + 119);
  }
  return hash;
}

int idStr::Hash(const char *string, int length) {
  int i, hash = 0;
  for (i = 0; i < length; i++) {
    hash += (*string++) * (i + 119);
  }
  return hash;
}

int idStr::IHash(const char *string) {
  int i, hash = 0;
  for(i = 0; *string != '\0'; i++) {
    hash += ToLower(*string++) * (i + 119);
  }
  return hash;
}

int idStr::IHash(const char *string, int length) {
  int i, hash = 0;
  for (i = 0; i < length; i++) {
    hash += ToLower(*string++) * (i + 119);
  }
  return hash;
}

char idStr::ToLower(char c) {
  if (c <= 'Z' && c >= 'A') {
    return (c + ('a' - 'A'));
  }
  return c;
}

char idStr::ToUpper(char c) {
  if (c >= 'a' && c <= 'z') {
    return (c - ('a' - 'A'));
  }
  return c;
}

bool idStr::CharIsPrintable(int c) {
  // test for regular ascii and western European high-ascii chars
  return (c >= 0x20 && c <= 0x7E) || (c >= 0xA1 && c <= 0xFF);
}

bool idStr::CharIsLower(int c) {
  // test for regular ascii and western European high-ascii chars
  return (c >= 'a' && c <= 'z') || (c >= 0xE0 && c <= 0xFF);
}

bool idStr::CharIsUpper(int c) {
  // test for regular ascii and western European high-ascii chars
  return (c <= 'Z' && c >= 'A') || (c >= 0xC0 && c <= 0xDF);
}

bool idStr::CharIsAlpha(int c) {
  // test for regular ascii and western European high-ascii chars
  return ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
      (c >= 0xC0 && c <= 0xFF));
}

bool idStr::CharIsNumeric(int c) {
  return (c <= '9' && c >= '0');
}

bool idStr::CharIsNewLine(char c) {
  return (c == '\n' || c == '\r' || c == '\v');
}

bool idStr::CharIsTab(char c) {
  return (c == '\t');
}

int idStr::DynamicMemoryUsed() const {
  return (data == baseBuffer) ? 0 : alloced;
}

const char *units[2][4] =
{
  { "B", "KB", "MB", "GB" },
  { "B/s", "KB/s", "MB/s", "GB/s" }
};

/*
   ============
   idStr::ReAllocate
   ============
   */
void idStr::ReAllocate(int amount, bool keepold) {
  char	*newbuffer;
  int		newsize;
  int		mod;

  //assert(data);
  assert(amount > 0);

  mod = amount % STR_ALLOC_GRAN;
  if (!mod) {
    newsize = amount;
  }
  else {
    newsize = amount + STR_ALLOC_GRAN - mod;
  }
  alloced = newsize;

#ifdef USE_STRING_DATA_ALLOCATOR
  newbuffer = stringDataAllocator.Alloc(alloced);
#else
  newbuffer = new char[ alloced ];
#endif
  if (keepold && data) {
    data[ len ] = '\0';
    strcpy(newbuffer, data);
  }

  if (data && data != baseBuffer) {
#ifdef USE_STRING_DATA_ALLOCATOR
    stringDataAllocator.Free(data);
#else
    delete [] data;
#endif
  }

  data = newbuffer;
}

/*
   ============
   idStr::FreeData
   ============
   */
void idStr::FreeData(void) {
  if (data && data != baseBuffer) {
#ifdef USE_STRING_DATA_ALLOCATOR
    stringDataAllocator.Free(data);
#else
    delete[] data;
#endif
    data = baseBuffer;
  }
}

/*
   ============
   idStr::operator=
   ============
   */
void idStr::operator=(const char *text) {
  int l;
  int diff;
  int i;

  if (!text) {
    // safe behaviour if NULL
    EnsureAlloced(1, false);
    data[ 0 ] = '\0';
    len = 0;
    return;
  }

  if (text == data) {
    return; // copying same thing
  }

  // check if we're aliasing
  if (text >= data && text <= data + len) {
    diff = text - data;

    assert(strlen(text) < (unsigned)len);

    for (i = 0; text[ i ]; i++) {
      data[ i ] = text[ i ];
    }

    data[ i ] = '\0';

    len -= diff;

    return;
  }

  l = strlen(text);
  EnsureAlloced(l + 1, false);
  strcpy(data, text);
  len = l;
}

/*
   ============
   idStr::FindChar

   returns -1 if not found otherwise the index of the char
   ============
   */
int idStr::FindChar(const char *str, const char c, int start, int end) {
  int i;

  if (end == -1) {
    end = strlen(str) - 1;
  }
  for (i = start; i <= end; i++) {
    if (str[i] == c) {
      return i;
    }
  }
  return -1;
}

/*
   ============
   idStr::FindText

   returns -1 if not found otherwise the index of the text
   ============
   */
int idStr::FindText(const char *str, const char *text, bool casesensitive, int start, int end) {
  int l, i, j;

  if (end == -1) {
    end = strlen(str);
  }
  l = end - strlen(text);
  for (i = start; i <= l; i++) {
    if (casesensitive) {
      for (j = 0; text[j]; j++) {
        if (str[i+j] != text[j]) {
          break;
        }
      }
    } else {
      for (j = 0; text[j]; j++) {
        if (::toupper(str[i+j]) != ::toupper(text[j])) {
          break;
        }
      }
    }
    if (!text[j]) {
      return i;
    }
  }
  return -1;
}

/*
   ============
   idStr::Filter

   Returns true if the string conforms the given filter.
   Several metacharacter may be used in the filter.

 *          match any string of zero or more characters
 ?          match any single character
 [abc...]   match any of the enclosed characters; a hyphen can
 be used to specify a range (e.g. a-z, A-Z, 0-9)

 ============
 */
bool idStr::Filter(const char *filter, const char *name, bool casesensitive) {
  idStr buf;
  int i, found, index;

  while(*filter) {
    if (*filter == '*') {
      filter++;
      buf.Empty();
      for (i = 0; *filter; i++) {
        if (*filter == '*' || *filter == '?' || (*filter == '[' && *(filter+1) != '[')) {
          break;
        }
        buf += *filter;
        if (*filter == '[') {
          filter++;
        }
        filter++;
      }
      if (buf.Length()) {
        index = idStr(name).Find(buf.c_str(), casesensitive);
        if (index == -1) {
          return false;
        }
        name += index + strlen(buf);
      }
    }
    else if (*filter == '?') {
      filter++;
      name++;
    }
    else if (*filter == '[') {
      if (*(filter+1) == '[') {
        if (*name != '[') {
          return false;
        }
        filter += 2;
        name++;
      }
      else {
        filter++;
        found = false;
        while(*filter && !found) {
          if (*filter == ']' && *(filter+1) != ']') {
            break;
          }
          if (*(filter+1) == '-' && *(filter+2) && (*(filter+2) != ']' || *(filter+3) == ']')) {
            if (casesensitive) {
              if (*name >= *filter && *name <= *(filter+2)) {
                found = true;
              }
            }
            else {
              if (::toupper(*name) >= ::toupper(*filter) && ::toupper(*name) <= ::toupper(*(filter+2))) {
                found = true;
              }
            }
            filter += 3;
          }
          else {
            if (casesensitive) {
              if (*filter == *name) {
                found = true;
              }
            }
            else {
              if (::toupper(*filter) == ::toupper(*name)) {
                found = true;
              }
            }
            filter++;
          }
        }
        if (!found) {
          return false;
        }
        while(*filter) {
          if (*filter == ']' && *(filter+1) != ']') {
            break;
          }
          filter++;
        }
        filter++;
        name++;
      }
    }
    else {
      if (casesensitive) {
        if (*filter != *name) {
          return false;
        }
      }
      else {
        if (::toupper(*filter) != ::toupper(*name)) {
          return false;
        }
      }
      filter++;
      name++;
    }
  }
  return true;
}

/*
   =============
   idStr::StripMediaName

   makes the string lower case, replaces backslashes with forward slashes, and removes extension
   =============
   */
void idStr::StripMediaName(const char *name, idStr &mediaName) {
  char c;

  mediaName.Empty();

  for (c = *name; c; c = *(++name)) {
    // truncate at an extension
    if (c == '.') {
      break;
    }
    // convert backslashes to forward slashes
    if (c == '\\') {
      mediaName.Append('/');
    } else {
      mediaName.Append(idStr::ToLower(c));
    }
  }
}

/*
   =============
   idStr::CheckExtension
   =============
   */
bool idStr::CheckExtension(const char *name, const char *ext) {
  const char *s1 = name + Length(name) - 1;
  const char *s2 = ext + Length(ext) - 1;
  int c1, c2, d;

  do {
    c1 = *s1--;
    c2 = *s2--;

    d = c1 - c2;
    while(d) {
      if (c1 <= 'Z' && c1 >= 'A') {
        d += ('a' - 'A');
        if (!d) {
          break;
        }
      }
      if (c2 <= 'Z' && c2 >= 'A') {
        d -= ('a' - 'A');
        if (!d) {
          break;
        }
      }
      return false;
    }
  } while(s1 > name && s2 > ext);

  return (s1 >= name);
}

/*
   =============
   idStr::FloatArrayToString
   =============
   */
const char *idStr::FloatArrayToString(const float *array, const int length, const int precision) {
  static int index = 0;
  static char str[4][16384];	// in case called by nested functions
  int i, n;
  char format[16], *s;

  // use an array of string so that multiple calls won't collide
  s = str[ index ];
  index = (index + 1) & 3;

  idStr::snPrintf(format, sizeof(format), "%%.%df", precision);
  n = idStr::snPrintf(s, sizeof(str[0]), format, (double)array[0]);
  if (precision > 0) {
    while(n > 0 && s[n-1] == '0') s[--n] = '\0';
    while(n > 0 && s[n-1] == '.') s[--n] = '\0';
  }
  idStr::snPrintf(format, sizeof(format), " %%.%df", precision);
  for (i = 1; i < length; i++) {
    n += idStr::snPrintf(s + n, sizeof(str[0]) - n, format, (double)array[i]);
    if (precision > 0) {
      while(n > 0 && s[n-1] == '0') s[--n] = '\0';
      while(n > 0 && s[n-1] == '.') s[--n] = '\0';
    }
  }
  return s;
}

/*
   ============
   idStr::Last

   returns -1 if not found otherwise the index of the char
   ============
   */
int idStr::Last(const char c) const {
  int i;

  for(i = Length(); i > 0; i--) {
    if (data[ i - 1 ] == c) {
      return i - 1;
    }
  }

  return -1;
}

/*
   ============
   idStr::StripLeading
   ============
   */
void idStr::StripLeading(const char c) {
  while(data[ 0 ] == c) {
    memmove(&data[ 0 ], &data[ 1 ], len);
    len--;
  }
}

/*
   ============
   idStr::StripLeading
   ============
   */
void idStr::StripLeading(const char *string) {
  int l;

  l = strlen(string);
  if (l > 0) {
    while (!Cmpn(string, l)) {
      memmove(data, data + l, len - l + 1);
      len -= l;
    }
  }
}

/*
   ============
   idStr::StripLeadingOnce
   ============
   */
bool idStr::StripLeadingOnce(const char *string) {
  int l;

  l = strlen(string);
  if ((l > 0) && !Cmpn(string, l)) {
    memmove(data, data + l, len - l + 1);
    len -= l;
    return true;
  }
  return false;
}

/*
   ============
   idStr::StripTrailing
   ============
   */
void idStr::StripTrailing(const char c) {
  int i;

  for(i = Length(); i > 0 && data[ i - 1 ] == c; i--) {
    data[ i - 1 ] = '\0';
    len--;
  }
}

/*
   ============
   idStr::StripLeading
   ============
   */
void idStr::StripTrailing(const char *string) {
  int l;

  l = strlen(string);
  if (l > 0) {
    while ((len >= l) && !Cmpn(string, data + len - l, l)) {
      len -= l;
      data[len] = '\0';
    }
  }
}

/*
   ============
   idStr::StripTrailingOnce
   ============
   */
bool idStr::StripTrailingOnce(const char *string) {
  int l;

  l = strlen(string);
  if ((l > 0) && (len >= l) && !Cmpn(string, data + len - l, l)) {
    len -= l;
    data[len] = '\0';
    return true;
  }
  return false;
}

/*
   ============
   idStr::Replace
   ============
   */
void idStr::Replace(const char *old, const char *nw) {
  int		oldLen, newLen, i, j, count;
  idStr	oldString(data);

  oldLen = strlen(old);
  newLen = strlen(nw);

  // Work out how big the new string will be
  count = 0;
  for(i = 0; i < oldString.Length(); i++) {
    if(!idStr::Cmpn(&oldString[i], old, oldLen)) {
      count++;
      i += oldLen - 1;
    }
  }

  if(count) {
    EnsureAlloced(len + ((newLen - oldLen) * count) + 2, false);

    // Replace the old data with the new data
    for(i = 0, j = 0; i < oldString.Length(); i++) {
      if(!idStr::Cmpn(&oldString[i], old, oldLen)) {
        memcpy(data + j, nw, newLen);
        i += oldLen - 1;
        j += newLen;
      } else {
        data[j] = oldString[i];
        j++;
      }
    }
    data[j] = 0;
    len = strlen(data);
  }
}

/*
   ============
   idStr::Mid
   ============
   */
const char *idStr::Mid(int start, int len, idStr &result) const {
  int i;

  result.Empty();

  i = Length();
  if (i == 0 || len <= 0 || start >= i) {
    return NULL;
  }

  if (start + len >= i) {
    len = i - start;
  }

  result.Append(&data[ start ], len);
  return result;
}

/*
   ============
   idStr::Mid
   ============
   */
idStr idStr::Mid(int start, int len) const {
  int i;
  idStr result;

  i = Length();
  if (i == 0 || len <= 0 || start >= i) {
    return result;
  }

  if (start + len >= i) {
    len = i - start;
  }

  result.Append(&data[ start ], len);
  return result;
}

/*
   ============
   idStr::StripTrailingWhitespace
   ============
   */
void idStr::StripTrailingWhitespace(void) {
  int i;

  // cast to unsigned char to prevent stripping off high-ASCII characters
  for(i = Length(); i > 0 && (unsigned char)(data[ i - 1 ]) <= ' '; i--) {
    data[ i - 1 ] = '\0';
    len--;
  }
}

/*
   ============
   idStr::StripQuotes

   Removes the quotes from the beginning and end of the string
   ============
   */
idStr& idStr::StripQuotes (void)
{
  if (data[0] != '\"')
  {
    return *this;
  }

  // Remove the trailing quote first
  if (data[len-1] == '\"')
  {
    data[len-1] = '\0';
    len--;
  }

  // Strip the leading quote now
  len--;
  memmove(&data[ 0 ], &data[ 1 ], len);
  data[len] = '\0';

  return *this;
}

/*
   =====================================================================

   filename methods

   =====================================================================
   */

/*
   ============
   idStr::FileNameHash
   ============
   */
int idStr::FileNameHash(void) const {
  int		i;
  int		hash;
  char	letter;

  hash = 0;
  i = 0;
  while(data[i] != '\0') {
    letter = idStr::ToLower(data[i]);
    if (letter == '.') {
      break;				// don't include extension
    }
    if (letter =='\\') {
      letter = '/';
    }
    hash += (int)(letter)*(i+119);
    i++;
  }
  hash &= (FILE_HASH_SIZE-1);
  return hash;
}

/*
   ============
   idStr::BackSlashesToSlashes
   ============
   */
idStr &idStr::BackSlashesToSlashes(void) {
  int i;

  for (i = 0; i < len; i++) {
    if (data[ i ] == '\\') {
      data[ i ] = '/';
    }
  }
  return *this;
}

/*
   ============
   idStr::SetFileExtension
   ============
   */
idStr &idStr::SetFileExtension(const char *extension) {
  StripFileExtension();
  if (*extension != '.') {
    Append('.');
  }
  Append(extension);
  return *this;
}

/*
   ============
   idStr::StripFileExtension
   ============
   */
idStr &idStr::StripFileExtension(void) {
  int i;

  for (i = len-1; i >= 0; i--) {
    if (data[i] == '.') {
      data[i] = '\0';
      len = i;
      break;
    }
  }
  return *this;
}

/*
   ============
   idStr::StripAbsoluteFileExtension
   ============
   */
idStr &idStr::StripAbsoluteFileExtension(void) {
  int i;

  for (i = 0; i < len; i++) {
    if (data[i] == '.') {
      data[i] = '\0';
      len = i;
      break;
    }
  }

  return *this;
}

/*
   ==================
   idStr::DefaultFileExtension
   ==================
   */
idStr &idStr::DefaultFileExtension(const char *extension) {
  int i;

  // do nothing if the string already has an extension
  for (i = len-1; i >= 0; i--) {
    if (data[i] == '.') {
      return *this;
    }
  }
  if (*extension != '.') {
    Append('.');
  }
  Append(extension);
  return *this;
}

/*
   ==================
   idStr::DefaultPath
   ==================
   */
idStr &idStr::DefaultPath(const char *basepath) {
  if (((*this)[ 0 ] == '/') || ((*this)[ 0 ] == '\\')) {
    // absolute path location
    return *this;
  }

  *this = basepath + *this;
  return *this;
}

/*
   ====================
   idStr::AppendPath
   ====================
   */
void idStr::AppendPath(const char *text) {
  int pos;
  int i = 0;

  if (text && text[i]) {
    pos = len;
    EnsureAlloced(len + strlen(text) + 2);

    if (pos) {
      if (data[ pos-1 ] != '/') {
        data[ pos++ ] = '/';
      }
    }
    if (text[i] == '/') {
      i++;
    }

    for (; text[ i ]; i++) {
      if (text[ i ] == '\\') {
        data[ pos++ ] = '/';
      } else {
        data[ pos++ ] = text[ i ];
      }
    }
    len = pos;
    data[ pos ] = '\0';
  }
}

/*
   ==================
   idStr::StripFilename
   ==================
   */
idStr &idStr::StripFilename(void) {
  int pos;

  pos = Length() - 1;
  while((pos > 0) && ((*this)[ pos ] != '/') && ((*this)[ pos ] != '\\')) {
    pos--;
  }

  if (pos < 0) {
    pos = 0;
  }

  CapLength(pos);
  return *this;
}

/*
   ==================
   idStr::StripPath
   ==================
   */
idStr &idStr::StripPath(void) {
  int pos;

  pos = Length();
  while((pos > 0) && ((*this)[ pos - 1 ] != '/') && ((*this)[ pos - 1 ] != '\\')) {
    pos--;
  }

  *this = Right(Length() - pos);
  return *this;
}

/*
   ====================
   idStr::ExtractFilePath
   ====================
   */
void idStr::ExtractFilePath(idStr &dest) const {
  int pos;

  //
  // back up until a \ or the start
  //
  pos = Length();
  while((pos > 0) && ((*this)[ pos - 1 ] != '/') && ((*this)[ pos - 1 ] != '\\')) {
    pos--;
  }

  Left(pos, dest);
}

/*
   ====================
   idStr::ExtractFileName
   ====================
   */
void idStr::ExtractFileName(idStr &dest) const {
  int pos;

  //
  // back up until a \ or the start
  //
  pos = Length() - 1;
  while((pos > 0) && ((*this)[ pos - 1 ] != '/') && ((*this)[ pos - 1 ] != '\\')) {
    pos--;
  }

  Right(Length() - pos, dest);
}

/*
   ====================
   idStr::ExtractFileBase
   ====================
   */
void idStr::ExtractFileBase(idStr &dest) const {
  int pos;
  int start;

  //
  // back up until a \ or the start
  //
  pos = Length() - 1;
  while((pos > 0) && ((*this)[ pos - 1 ] != '/') && ((*this)[ pos - 1 ] != '\\')) {
    pos--;
  }

  start = pos;
  while((pos < Length()) && ((*this)[ pos ] != '.')) {
    pos++;
  }

  Mid(start, pos - start, dest);
}

/*
   ====================
   idStr::ExtractFileExtension
   ====================
   */
void idStr::ExtractFileExtension(idStr &dest) const {
  int pos;

  //
  // back up until a . or the start
  //
  pos = Length() - 1;
  while((pos > 0) && ((*this)[ pos - 1 ] != '.')) {
    pos--;
  }

  if (!pos) {
    // no extension
    dest.Empty();
  } else {
    Right(Length() - pos, dest);
  }
}


/*
   =====================================================================

   char * methods to replace library functions

   =====================================================================
   */

/*
   ============
   idStr::IsNumeric

   Checks a string to see if it contains only numerical values.
   ============
   */
bool idStr::IsNumeric(const char *s) {
  int		i;
  bool	dot;

  if (*s == '-') {
    s++;
  }

  dot = false;
  for (i = 0; s[i]; i++) {
    if (!isdigit(s[i])) {
      if ((s[ i ] == '.') && !dot) {
        dot = true;
        continue;
      }
      return false;
    }
  }

  return true;
}

/*
   ============
   idStr::HasLower

   Checks if a string has any lowercase chars
   ============
   */
bool idStr::HasLower(const char *s) {
  if (!s) {
    return false;
  }

  while (*s) {
    if (CharIsLower(*s)) {
      return true;
    }
    s++;
  }

  return false;
}

/*
   ============
   idStr::HasUpper

   Checks if a string has any uppercase chars
   ============
   */
bool idStr::HasUpper(const char *s) {
  if (!s) {
    return false;
  }

  while (*s) {
    if (CharIsUpper(*s)) {
      return true;
    }
    s++;
  }

  return false;
}

/*
   ================
   idStr::Cmp
   ================
   */
int idStr::Cmp(const char *s1, const char *s2) {
  int c1, c2, d;

  do {
    c1 = *s1++;
    c2 = *s2++;

    d = c1 - c2;
    if (d) {
      return (INTSIGNBITNOTSET(d) << 1) - 1;
    }
  } while(c1);

  return 0;		// strings are equal
}

/*
   ================
   idStr::Cmpn
   ================
   */
int idStr::Cmpn(const char *s1, const char *s2, int n) {
  int c1, c2, d;

  assert(n >= 0);

  do {
    c1 = *s1++;
    c2 = *s2++;

    if (!n--) {
      return 0;		// strings are equal until end point
    }

    d = c1 - c2;
    if (d) {
      return (INTSIGNBITNOTSET(d) << 1) - 1;
    }
  } while(c1);

  return 0;		// strings are equal
}

/*
   ================
   idStr::Icmp
   ================
   */
int idStr::Icmp(const char *s1, const char *s2) {
  int c1, c2, d;

  do {
    c1 = *s1++;
    c2 = *s2++;

    d = c1 - c2;
    while(d) {
      if (c1 <= 'Z' && c1 >= 'A') {
        d += ('a' - 'A');
        if (!d) {
          break;
        }
      }
      if (c2 <= 'Z' && c2 >= 'A') {
        d -= ('a' - 'A');
        if (!d) {
          break;
        }
      }
      return (INTSIGNBITNOTSET(d) << 1) - 1;
    }
  } while(c1);

  return 0;		// strings are equal
}

/*
   ================
   idStr::Icmpn
   ================
   */
int idStr::Icmpn(const char *s1, const char *s2, int n) {
  int c1, c2, d;

  assert(n >= 0);

  do {
    c1 = *s1++;
    c2 = *s2++;

    if (!n--) {
      return 0;		// strings are equal until end point
    }

    d = c1 - c2;
    while(d) {
      if (c1 <= 'Z' && c1 >= 'A') {
        d += ('a' - 'A');
        if (!d) {
          break;
        }
      }
      if (c2 <= 'Z' && c2 >= 'A') {
        d -= ('a' - 'A');
        if (!d) {
          break;
        }
      }
      return (INTSIGNBITNOTSET(d) << 1) - 1;
    }
  } while(c1);

  return 0;		// strings are equal
}

/*
   ================
   idStr::IcmpPath
   ================
   */
int idStr::IcmpPath(const char *s1, const char *s2) {
  int c1, c2, d;

#if 0
  //#if !defined(_WIN32)
  idLib::common->Printf("WARNING: IcmpPath used on a case-sensitive filesystem?\n");
#endif

  do {
    c1 = *s1++;
    c2 = *s2++;

    d = c1 - c2;
    while(d) {
      if (c1 <= 'Z' && c1 >= 'A') {
        d += ('a' - 'A');
        if (!d) {
          break;
        }
      }
      if (c1 == '\\') {
        d += ('/' - '\\');
        if (!d) {
          break;
        }
      }
      if (c2 <= 'Z' && c2 >= 'A') {
        d -= ('a' - 'A');
        if (!d) {
          break;
        }
      }
      if (c2 == '\\') {
        d -= ('/' - '\\');
        if (!d) {
          break;
        }
      }
      // make sure folders come first
      while(c1) {
        if (c1 == '/' || c1 == '\\') {
          break;
        }
        c1 = *s1++;
      }
      while(c2) {
        if (c2 == '/' || c2 == '\\') {
          break;
        }
        c2 = *s2++;
      }
      if (c1 && !c2) {
        return -1;
      } else if (!c1 && c2) {
        return 1;
      }
      // same folder depth so use the regular compare
      return (INTSIGNBITNOTSET(d) << 1) - 1;
    }
  } while(c1);

  return 0;
}

/*
   ================
   idStr::IcmpnPath
   ================
   */
int idStr::IcmpnPath(const char *s1, const char *s2, int n) {
  int c1, c2, d;

#if 0
  //#if !defined(_WIN32)
  idLib::common->Printf("WARNING: IcmpPath used on a case-sensitive filesystem?\n");
#endif

  assert(n >= 0);

  do {
    c1 = *s1++;
    c2 = *s2++;

    if (!n--) {
      return 0;		// strings are equal until end point
    }

    d = c1 - c2;
    while(d) {
      if (c1 <= 'Z' && c1 >= 'A') {
        d += ('a' - 'A');
        if (!d) {
          break;
        }
      }
      if (c1 == '\\') {
        d += ('/' - '\\');
        if (!d) {
          break;
        }
      }
      if (c2 <= 'Z' && c2 >= 'A') {
        d -= ('a' - 'A');
        if (!d) {
          break;
        }
      }
      if (c2 == '\\') {
        d -= ('/' - '\\');
        if (!d) {
          break;
        }
      }
      // make sure folders come first
      while(c1) {
        if (c1 == '/' || c1 == '\\') {
          break;
        }
        c1 = *s1++;
      }
      while(c2) {
        if (c2 == '/' || c2 == '\\') {
          break;
        }
        c2 = *s2++;
      }
      if (c1 && !c2) {
        return -1;
      } else if (!c1 && c2) {
        return 1;
      }
      // same folder depth so use the regular compare
      return (INTSIGNBITNOTSET(d) << 1) - 1;
    }
  } while(c1);

  return 0;
}

/*
   =============
   idStr::Copynz

   Safe strncpy that ensures a trailing zero
   =============
   */
void idStr::Copynz(char *dest, const char *src, int destsize) {
  if (!src) {
    warning_ln("idStr::Copynz: NULL src");
    return;
  }
  if (destsize < 1) {
    warning_ln("idStr::Copynz: destsize < 1");
    return;
  }

  strncpy(dest, src, destsize-1);
  dest[destsize-1] = 0;
}

/*
   ================
   idStr::Append

   never goes past bounds or leaves without a terminating 0
   ================
   */
void idStr::Append(char *dest, int size, const char *src) {
  int		l1;

  l1 = strlen(dest);
  if (l1 >= size) {
    die("idStr::Append: already overflowed");
  }
  idStr::Copynz(dest + l1, src, size - l1);
}

/*
   ================
   idStr::snPrintf
   ================
   */
int idStr::snPrintf(char *dest, int size, const char *fmt, ...) {
  int len;
  va_list argptr;
  char buffer[32000];	// big, but small enough to fit in PPC stack

  va_start(argptr, fmt);
  len = vsprintf(buffer, fmt, argptr);
  va_end(argptr);
  if ((unsigned int)len >= sizeof(buffer)) {
    die("idStr::snPrintf: overflowed buffer");
  }
  if (len >= size) {
    warning_ln("idStr::snPrintf: overflow of %i in %i\n", len, size);
    len = size;
  }
  idStr::Copynz(dest, buffer, size);
  return len;
}

/*
   ============
   idStr::vsnPrintf

   vsnprintf portability:

   C99 standard: vsnprintf returns the number of characters (excluding the trailing
   '\0') which would have been written to the final string if enough space had been available
   snprintf and vsnprintf do not write more than size bytes (including the trailing '\0')

win32: _vsnprintf returns the number of characters written, not including the terminating null character,
or a negative value if an output error occurs. If the number of characters to write exceeds count, then count
characters are written and -1 is returned and no trailing '\0' is added.

idStr::vsnPrintf: always appends a trailing '\0', returns number of characters written (not including terminal \0)
or returns -1 on failure or if the buffer would be overflowed.
============
*/
int idStr::vsnPrintf(char *dest, int size, const char *fmt, va_list argptr) {
  int ret;

#ifdef _WIN32
#undef _vsnprintf
  ret = _vsnprintf(dest, size-1, fmt, argptr);
#define _vsnprintf	use_idStr_vsnPrintf
#else
#undef vsnprintf
  ret = vsnprintf(dest, size, fmt, argptr);
#define vsnprintf	use_idStr_vsnPrintf
#endif
  dest[size-1] = '\0';
  if (ret < 0 || ret >= size) {
    return -1;
  }
  return ret;
}

/*
   ============
   sprintf

   Sets the value of the string using a printf interface.
   ============
   */
int sprintf(idStr &string, const char *fmt, ...) {
  int l;
  va_list argptr;
  char buffer[32000];

  va_start(argptr, fmt);
  l = idStr::vsnPrintf(buffer, sizeof(buffer)-1, fmt, argptr);
  va_end(argptr);
  buffer[sizeof(buffer)-1] = '\0';

  string = buffer;
  return l;
}

/*
   ============
   vsprintf

   Sets the value of the string using a vprintf interface.
   ============
   */
int vsprintf(idStr &string, const char *fmt, va_list argptr) {
  int l;
  char buffer[32000];

  l = idStr::vsnPrintf(buffer, sizeof(buffer)-1, fmt, argptr);
  buffer[sizeof(buffer)-1] = '\0';

  string = buffer;
  return l;
}

/*
   ============
   va

   does a varargs printf into a temp buffer
NOTE: not thread safe
============
*/
char *va(const char *fmt, ...) {
  va_list argptr;
  static int index = 0;
  static char string[4][16384];	// in case called by nested functions
  char *buf;

  buf = string[index];
  index = (index + 1) & 3;

  va_start(argptr, fmt);
  vsprintf(buf, fmt, argptr);
  va_end(argptr);

  return buf;
}

/*
   ================
   idStr::InitMemory
   ================
   */
void idStr::InitMemory(void) {
#ifdef USE_STRING_DATA_ALLOCATOR
  stringDataAllocator.Init();
#endif
}

/*
   ================
   idStr::ShutdownMemory
   ================
   */
void idStr::ShutdownMemory(void) {
#ifdef USE_STRING_DATA_ALLOCATOR
  stringDataAllocator.Shutdown();
#endif
}

/*
   ================
   idStr::PurgeMemory
   ================
   */
void idStr::PurgeMemory(void) {
#ifdef USE_STRING_DATA_ALLOCATOR
  stringDataAllocator.FreeEmptyBaseBlocks();
#endif
}

/*
   ================
   idStr::FormatNumber
   ================
   */
struct formatList_t {
  int			gran;
  int			count;
};

// elements of list need to decend in size
formatList_t formatList[] = {
  { 1000000000, 0 },
  { 1000000, 0 },
  { 1000, 0 }
};

int numFormatList = sizeof(formatList) / sizeof(formatList[0]);


idStr idStr::FormatNumber(int number) {
  idStr string;
  bool hit;

  // reset
  for (int i = 0; i < numFormatList; i++) {
    formatList_t *li = formatList + i;
    li->count = 0;
  }

  // main loop
  do {
    hit = false;

    for (int i = 0; i < numFormatList; i++) {
      formatList_t *li = formatList + i;

      if (number >= li->gran) {
        li->count++;
        number -= li->gran;
        hit = true;
        break;
      }
    }
  } while (hit);

  // print out
  bool found = false;

  for (int i = 0; i < numFormatList; i++) {
    formatList_t *li = formatList + i;

    if (li->count) {
      if (!found) {
        string += va("%i,", li->count);
      } else {
        string += va("%3.3i,", li->count);
      }
      found = true;
    }
    else if (found) {
      string += va("%3.3i,", li->count);
    }
  }

  if (found) {
    string += va("%3.3i", number);
  }
  else {
    string += va("%i", number);
  }

  // pad to proper size
  int count = 11 - string.Length();

  for (int i = 0; i < count; i++) {
    string.Insert(" ", 0);
  }

  return string;
}

