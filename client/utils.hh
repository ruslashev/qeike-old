#pragma once

#include <cstdarg>
#include <string>
#include <fstream>

#ifdef __PRETTY_FUNCTION__
#define info() \
  printf("%s:%d in %s", __FILE__, __LINE__, __PRETTY_FUNCTION__);
#else
#define info() \
  printf("%s:%d in %s", __FILE__, __LINE__, __func__);
#endif

#define assertf(X, ...) \
  do { \
    if (!(X)) { \
      printf("assertion `%s' failed at ", #X); \
      info() \
      printf(": "); \
      printf(__VA_ARGS__); \
      puts(""); \
      exit(1); \
    } \
  } while (0)

#define die(...) do { printf(__VA_ARGS__); puts(""); exit(1); } while (0)

#define normal_exit(...) do { printf(__VA_ARGS__); puts(""); exit(0); } while (0)

#define warning(...) do { \
  printf("warning: "); \
  printf(__VA_ARGS__); \
  puts(""); \
} while (0)

#define warning_ln(...) do { \
  printf("warning at "); \
  info() \
  printf(": "); \
  printf(__VA_ARGS__); \
  puts(""); \
} while (0)

#define _glsl(X) "#version 120\n" #X

inline void print_packet(uint8_t *packet, size_t len
    , const char *msg = "packet") {
  printf("%s: ", msg);
  for (size_t i = 0; i < len; i++) {
    int numbits = 8;
    while (--numbits >= 0)
      printf("%c", (packet[i] & ((uint8_t)1 << numbits)) ? '1' : '0');
    printf(" ");
  }
  printf("\n");
}

template <typename T>
inline T clamp(T value, T min, T max) {
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

struct file_reader {
  size_t length;
  char *buffer;

  file_reader()
    : length(0)
    , buffer(nullptr) {
  }

  file_reader(const std::string &filename) {
    open(filename);
  }

  void open(const std::string &filename) {
    // printf("filename=%s\n", filename.c_str());
    std::ifstream ifs(filename, std::ios::binary);
    assertf(ifs, "failed to open file \"%s\"", filename.c_str());

    ifs.seekg(0, ifs.end);
    length = ifs.tellg();
    ifs.seekg(0, ifs.beg);
    buffer = new char [length + 1];
    ifs.read(buffer, length);
    buffer[length] = 0;
    ifs.close();
  }

  ~file_reader() {
    delete[] buffer;
  }
};

