#pragma once

#include <cstring>
#include <cstdarg>
#include <cmath>
#include "../utils.hh"

#define PATHSEPERATOR_STR			"/"
#define PATHSEPERATOR_CHAR			'/'

const int STR_ALLOC_BASE   = 20;
const int STR_ALLOC_GRAN   = 32;

#define FILE_HASH_SIZE  1024

#define assert(X) if (X) {} else assertf(0, "fail")

#define INTSIGNBITNOTSET(i)  ((~((const unsigned int)(i))) >> 31)

class idStr {

public:
  idStr( void );
  idStr( const idStr &text );
  idStr( const idStr &text, int start, int end );
  idStr( const char *text );
  idStr( const char *text, int start, int end );
  explicit idStr( const bool b );
  explicit idStr( const char c );
  explicit idStr( const int i );
  explicit idStr( const unsigned u );
  explicit idStr( const float f );
  ~idStr( void );

  size_t    Size( void ) const;
  const char *  c_str( void ) const;
  operator   const char *( void ) const;
  operator   const char *( void );

  char    operator[]( int index ) const;
  char &    operator[]( int index );

  void    operator=( const idStr &text );
  void    operator=( const char *text );

  friend idStr  operator+( const idStr &a, const idStr &b );
  friend idStr  operator+( const idStr &a, const char *b );
  friend idStr  operator+( const char *a, const idStr &b );

  friend idStr  operator+( const idStr &a, const float b );
  friend idStr  operator+( const idStr &a, const int b );
  friend idStr  operator+( const idStr &a, const unsigned b );
  friend idStr  operator+( const idStr &a, const bool b );
  friend idStr  operator+( const idStr &a, const char b );

  idStr &    operator+=( const idStr &a );
  idStr &    operator+=( const char *a );
  idStr &    operator+=( const float a );
  idStr &    operator+=( const char a );
  idStr &    operator+=( const int a );
  idStr &    operator+=( const unsigned a );
  idStr &    operator+=( const bool a );

  // case sensitive compare
  friend bool   operator==( const idStr &a, const idStr &b );
  friend bool   operator==( const idStr &a, const char *b );
  friend bool   operator==( const char *a, const idStr &b );

  // case sensitive compare
  friend bool   operator!=( const idStr &a, const idStr &b );
  friend bool   operator!=( const idStr &a, const char *b );
  friend bool   operator!=( const char *a, const idStr &b );

  // case sensitive compare
  int     Cmp( const char *text ) const;
  int     Cmpn( const char *text, int n ) const;
  int     CmpPrefix( const char *text ) const;

  // case insensitive compare
  int     Icmp( const char *text ) const;
  int     Icmpn( const char *text, int n ) const;
  int     IcmpPrefix( const char *text ) const;

  // compares paths and makes sure folders come first
  int     IcmpPath( const char *text ) const;
  int     IcmpnPath( const char *text, int n ) const;
  int     IcmpPrefixPath( const char *text ) const;

  int     Length( void ) const;
  int     Allocated( void ) const;
  void    Empty( void );
  bool    IsEmpty( void ) const;
  void    Clear( void );
  void    Append( const char a );
  void    Append( const idStr &text );
  void    Append( const char *text );
  void    Append( const char *text, int len );
  void    Insert( const char a, int index );
  void    Insert( const char *text, int index );
  void    ToLower( void );
  void    ToUpper( void );
  bool    IsNumeric( void ) const;
  bool    HasLower( void ) const;
  bool    HasUpper( void ) const;
  void    CapLength( int );
  void    Fill( const char ch, int newlen );

  int     Find( const char c, int start = 0, int end = -1 ) const;
  int     Find( const char *text, bool casesensitive = true, int start = 0, int end = -1 ) const;
  bool    Filter( const char *filter, bool casesensitive ) const;
  int     Last( const char c ) const;      // return the index to the last occurance of 'c', returns -1 if not found
  const char *  Left( int len, idStr &result ) const;   // store the leftmost 'len' characters in the result
  const char *  Right( int len, idStr &result ) const;   // store the rightmost 'len' characters in the result
  const char *  Mid( int start, int len, idStr &result ) const; // store 'len' characters starting at 'start' in result
  idStr    Left( int len ) const;       // return the leftmost 'len' characters
  idStr    Right( int len ) const;       // return the rightmost 'len' characters
  idStr    Mid( int start, int len ) const;    // return 'len' characters starting at 'start'
  void    StripLeading( const char c );     // strip char from front as many times as the char occurs
  void    StripLeading( const char *string );    // strip string from front as many times as the string occurs
  bool    StripLeadingOnce( const char *string );   // strip string from front just once if it occurs
  void    StripTrailing( const char c );     // strip char from end as many times as the char occurs
  void    StripTrailing( const char *string );   // strip string from end as many times as the string occurs
  bool    StripTrailingOnce( const char *string );  // strip string from end just once if it occurs
  void    Strip( const char c );       // strip char from front and end as many times as the char occurs
  void    Strip( const char *string );     // strip string from front and end as many times as the string occurs
  void    StripTrailingWhitespace( void );    // strip trailing white space characters
  idStr &    StripQuotes( void );       // strip quotes around string
  void    Replace( const char *old, const char *nw );

  // file name methods
  int     FileNameHash( void ) const;      // hash key for the filename (skips extension)
  idStr &    BackSlashesToSlashes( void );     // convert slashes
  idStr &    SetFileExtension( const char *extension );  // set the given file extension
  idStr &    StripFileExtension( void );      // remove any file extension
  idStr &    StripAbsoluteFileExtension( void );    // remove any file extension looking from front (useful if there are multiple .'s)
  idStr &    DefaultFileExtension( const char *extension ); // if there's no file extension use the default
  idStr &    DefaultPath( const char *basepath );   // if there's no path use the default
  void    AppendPath( const char *text );     // append a partial path
  idStr &    StripFilename( void );       // remove the filename from a path
  idStr &    StripPath( void );        // remove the path from the filename
  void    ExtractFilePath( idStr &dest ) const;   // copy the file path to another string
  void    ExtractFileName( idStr &dest ) const;   // copy the filename to another string
  void    ExtractFileBase( idStr &dest ) const;   // copy the filename minus the extension to another string
  void    ExtractFileExtension( idStr &dest ) const;  // copy the file extension to another string
  bool    CheckExtension( const char *ext );

  // char * methods to replace library functions
  static int   Length( const char *s );
  static char *  ToLower( char *s );
  static char *  ToUpper( char *s );
  static bool   IsNumeric( const char *s );
  static bool   HasLower( const char *s );
  static bool   HasUpper( const char *s );
  static int   Cmp( const char *s1, const char *s2 );
  static int   Cmpn( const char *s1, const char *s2, int n );
  static int   Icmp( const char *s1, const char *s2 );
  static int   Icmpn( const char *s1, const char *s2, int n );
  static int   IcmpPath( const char *s1, const char *s2 );   // compares paths and makes sure folders come first
  static int   IcmpnPath( const char *s1, const char *s2, int n ); // compares paths and makes sure folders come first
  static void   Append( char *dest, int size, const char *src );
  static void   Copynz( char *dest, const char *src, int destsize );
  static int   snPrintf( char *dest, int size, const char *fmt, ... );
  static int   vsnPrintf( char *dest, int size, const char *fmt, va_list argptr );
  static int   FindChar( const char *str, const char c, int start = 0, int end = -1 );
  static int   FindText( const char *str, const char *text, bool casesensitive = true, int start = 0, int end = -1 );
  static bool   Filter( const char *filter, const char *name, bool casesensitive );
  static void   StripMediaName( const char *name, idStr &mediaName );
  static bool   CheckExtension( const char *name, const char *ext );
  static const char * FloatArrayToString( const float *array, const int length, const int precision );

  // hash keys
  static int   Hash( const char *string );
  static int   Hash( const char *string, int length );
  static int   IHash( const char *string );     // case insensitive
  static int   IHash( const char *string, int length );  // case insensitive

  // character methods
  static char   ToLower( char c );
  static char   ToUpper( char c );
  static bool   CharIsPrintable( int c );
  static bool   CharIsLower( int c );
  static bool   CharIsUpper( int c );
  static bool   CharIsAlpha( int c );
  static bool   CharIsNumeric( int c );
  static bool   CharIsNewLine( char c );
  static bool   CharIsTab( char c );

  friend int   sprintf( idStr &dest, const char *fmt, ... );
  friend int   vsprintf( idStr &dest, const char *fmt, va_list ap );

  void    ReAllocate( int amount, bool keepold );    // reallocate string data buffer
  void    FreeData( void );         // free allocated string memory

  static void   InitMemory( void );
  static void   ShutdownMemory( void );
  static void   PurgeMemory( void );

  int     DynamicMemoryUsed() const;
  static idStr  FormatNumber( int number );

  int     len;
  char *    data;
  int     alloced;
  char    baseBuffer[ STR_ALLOC_BASE ];

  void    Init( void );          // initialize string using base buffer
  void    EnsureAlloced( int amount, bool keepold = true ); // ensure string data buffer is large anough
};

char *     va( const char *fmt, ... );

