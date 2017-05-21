#pragma once

#include "utils.hh"
#include <GL/glew.h>
#include <vector>
#include <string>
#include <cstring>

void gl_error_description(GLenum err) {
  switch (err) {
    case GL_INVALID_ENUM:
      puts("GL_INVALID_ENUM: An unacceptable value is specified for an\n"
          "enumerated argument. The offending command is ignored and has\n"
          "no other side effect than to set the error flag.");
      break;
    case GL_INVALID_VALUE:
      puts("GL_INVALID_VALUE: A numeric argument is out of range. The\n"
          "offending command is ignored and has no other side effect than\n"
          "to set the error flag.");
      break;
    case GL_INVALID_OPERATION:
      puts("GL_INVALID_OPERATION: The specified operation is not allowed in\n"
          "the current state. The offending command is ignored and has no\n"
          "other side effect than to set the error flag.");
      break;
    case GL_INVALID_FRAMEBUFFER_OPERATION:
      puts("GL_INVALID_FRAMEBUFFER_OPERATION: The framebuffer object is not\n"
          "complete. The offending command is ignored and has no other side\n"
          "effect than to set the error flag.");
      break;
    case GL_OUT_OF_MEMORY:
      puts("GL_OUT_OF_MEMORY: There is not enough memory left to execute the\n"
          "command. The state of the GL is undefined, except for the state\n"
          "of the error flags, after this error is recorded.");
      break;
    case GL_STACK_UNDERFLOW:
      puts("GL_STACK_UNDERFLOW: An attempt has been made to perform an\n"
          "operation that would cause an internal stack to underflow.");
      break;
    case GL_STACK_OVERFLOW:
      puts("GL_STACK_OVERFLOW: An attempt has been made to perform an\n"
          "operation that would cause an internal stack to overflow.");
      break;
    default:
      break;
  }
}

#define gl_check_errors() do { \
  GLenum err = glGetError(); \
  while (err != GL_NO_ERROR) { \
    fprintf(stderr, "glError %d at %s:%d\n", err, __FILE__, __LINE__); \
    gl_error_description(err); \
    err = glGetError(); \
  } \
} while (0)

class ogl_buffer {
protected:
  GLuint _id;
  GLenum _type;
public:
  ogl_buffer(GLenum n_type) : _type(n_type) {
    glGenBuffers(1, &_id);
  }
  ~ogl_buffer() {
    glDeleteBuffers(1, &_id);
  }
  void bind() const {
    glBindBuffer(_type, _id);
  }
  void unbind() const {
    glBindBuffer(_type, 0);
  }
};

class array_buffer : public ogl_buffer {
public:
  array_buffer() : ogl_buffer(GL_ARRAY_BUFFER) {}
  void upload(const std::vector<GLfloat> &data) {
    glBufferData(GL_ARRAY_BUFFER
        , static_cast<GLsizeiptr>(data.size() * sizeof(data[0])), data.data()
        , GL_STATIC_DRAW);
  }
};

class element_array_buffer : public ogl_buffer {
public:
  element_array_buffer() : ogl_buffer(GL_ELEMENT_ARRAY_BUFFER) {}
  void upload(const std::vector<GLushort> &data) {
    glBufferData(GL_ELEMENT_ARRAY_BUFFER
        , static_cast<GLsizeiptr>(data.size() * sizeof(data[0])), data.data()
        , GL_STATIC_DRAW);
  }
};

static std::string get_ogl_shader_err(GLint loglen
    , void (*ogl_errmsg_func)(GLuint, GLsizei, GLsizei*, GLchar*)
    , GLuint id) {
  char msg[loglen + 1];
  ogl_errmsg_func(id, loglen, nullptr, msg);
  msg[loglen] = 0;
  std::string msgstr(msg);
  msgstr.pop_back(); // strip trailing newline
  /*
  // indent every line
  int indent = 3;
  msgstr.insert(msgstr.begin(), indent, ' ');
  for (size_t i = 0; i < msgstr.size(); i++) {
    if (msgstr[i] != '\n')
      continue;
    msgstr.insert(i, indent, ' ');
  }
  */
  return msgstr;
}

struct shader {
  GLuint type;
  GLuint id;
  shader(std::string source, GLuint ntype) : type(ntype) {
    id = glCreateShader(type);
    const char *csrc = source.c_str();
    glShaderSource(id, 1, &csrc, NULL);
    glCompileShader(id);
    GLint loglen;
    glGetShaderiv(id, GL_INFO_LOG_LENGTH, &loglen);
    if (loglen != 0) {
      std::string msg = get_ogl_shader_err(loglen, glGetShaderInfoLog, id);
      GLint compilesucc;
      glGetShaderiv(id, GL_COMPILE_STATUS, &compilesucc);
      if (compilesucc != GL_TRUE)
        die("failed to compile %s shader:\n###\n%s###"
            , type == GL_VERTEX_SHADER ? "vertex" : "fragment"
            , msg.c_str());
      else
        printf("### %s shader diagnostic message:\n%s\n### diagnostic message end\n"
            , type == GL_VERTEX_SHADER ? "vertex" : "fragment"
            , msg.c_str());
    }
  }
  ~shader() {
    glDeleteShader(id);
  }
};

struct shaderprogram {
  GLuint id;
  shaderprogram(const shader &vert, const shader &frag) {
    assertf(vert.type == GL_VERTEX_SHADER && frag.type == GL_FRAGMENT_SHADER
          , "order of shaders in shaderprogram's constructor is reversed");
    id = glCreateProgram();
    glAttachShader(id, vert.id);
    glAttachShader(id, frag.id);
    glLinkProgram(id);
    GLint loglen;
    glGetProgramiv(id, GL_INFO_LOG_LENGTH, &loglen);
    if (loglen != 0) {
      std::string msg = get_ogl_shader_err(loglen, glGetProgramInfoLog, id);
      GLint linksucc;
      glGetProgramiv(id, GL_LINK_STATUS, &linksucc);
      if (linksucc != GL_TRUE) {
        glDetachShader(id, vert.id);
        glDetachShader(id, frag.id);
        die("failed to link a program:\n%s"
            , msg.c_str());
      } else
        printf("### shader program diagnostic message:\n%s\n### diagnostic message end\n"
            , msg.c_str());
    }
  }
  ~shaderprogram() {
    glDeleteProgram(id);
  }
  void vertexattribptr(const array_buffer &buffer, const char *name,
      GLint size, GLenum type, GLboolean normalized, GLsizei stride,
      const GLvoid *ptr) {
    buffer.bind();
    GLint attr = glGetAttribLocation(id, name);
    glEnableVertexAttribArray(static_cast<GLuint>(attr));
    glVertexAttribPointer(static_cast<GLuint>(attr), size, type, normalized
        , stride, ptr);
    buffer.unbind();
  }
  GLuint bind_attrib(const char *name) {
    GLint prev_active_prog;
    glGetIntegerv(GL_CURRENT_PROGRAM, &prev_active_prog);
    use_this_prog();
    GLint attr = glGetAttribLocation(id, name);
    assertf(glGetError() == GL_NO_ERROR, "failed to bind attribute %s", name);
    if (attr == -1)
      printf("warning: failed to bind attribute %s\n", name);
    glUseProgram(static_cast<GLuint>(prev_active_prog));
    return static_cast<GLuint>(attr);
  }
  GLint bind_uniform(const char *name) {
    GLint prev_active_prog;
    glGetIntegerv(GL_CURRENT_PROGRAM, &prev_active_prog);
    use_this_prog();
    GLint unif = glGetUniformLocation(id, name);
    assertf(glGetError() == GL_NO_ERROR, "failed to bind uniform %s", name);
    if (0) // (unif == -1)
      printf("warning: failed to bind uniform %s\n", name);
    glUseProgram(static_cast<GLuint>(prev_active_prog));
    return unif;
  }
  void use_this_prog() {
    glUseProgram(id);
  }
  void dont_use_this_prog() {
    glUseProgram(0);
  }
};

struct vertexarray {
  GLuint id;
  vertexarray() {
    glGenVertexArrays(1, &id);
  }
  ~vertexarray() {
    glDeleteVertexArrays(1, &id);
  }
  void bind() const {
    glBindVertexArray(id);
  }
  void unbind() const {
    glBindVertexArray(0);
  }
};

