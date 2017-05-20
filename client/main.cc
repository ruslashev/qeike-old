#include "camera.hh"
#include "ogl.hh"
#include "screen.hh"
#include "utils.hh"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

static shaderprogram *sp;
static GLuint vattr;
static array_buffer *cube_vbuf;
static element_array_buffer *cube_ebuf;
static GLint resolution_unif, time_unif, modelmat_unif, viewmat_unif, color_unif;
static shader *vs, *fs;
static vertexarray *vao;
static camera *cam;
static float fov = 106.26f;
static int move, strafe;

static void graphics_load(screen *s) {
  s->lock_mouse();

  glEnable(GL_DEPTH_TEST);

  glClearColor(0.055f, 0.055f, 0.055f, 1);

  const char *vsrc = _glsl(
    attribute vec3 vertex_pos;
    uniform mat4 model;
    uniform mat4 view;
    uniform mat4 projection;
    void main() {
      gl_Position = projection * view * model * vec4(vertex_pos, 1.0);
    }
  );
  const char *fsrc = _glsl(
    uniform vec2 iResolution;
    uniform float iGlobalTime;
    uniform vec3 color;
    void main() {
      gl_FragColor = vec4(color, 1.0);
    }
  );

  vs = new shader(vsrc, GL_VERTEX_SHADER);
  fs = new shader(fsrc, GL_FRAGMENT_SHADER);
  sp = new shaderprogram(*vs, *fs);

  vattr = sp->bind_attrib("vertex_pos");
  resolution_unif = sp->bind_uniform("iResolution");
  time_unif = sp->bind_uniform("iGlobalTime");
  modelmat_unif = sp->bind_uniform("model");
  viewmat_unif = sp->bind_uniform("view");
  color_unif = sp->bind_uniform("color");

  vao = new vertexarray;
  cube_vbuf = new array_buffer;
  cube_ebuf = new element_array_buffer;
  vao->bind();
  cube_vbuf->bind();
  const std::vector<float> cube_verts = {
    -1, -1,  1, 1, -1,  1, 1,  1,  1, -1,  1,  1,
    -1, -1, -1, 1, -1, -1, 1,  1, -1, -1,  1, -1,
  };
  cube_vbuf->upload(cube_verts);
  glVertexAttribPointer(vattr, 3, GL_FLOAT, GL_FALSE, 0, 0);
  cube_ebuf->bind();
  const std::vector<GLushort> cube_elements = {
    0, 1, 2, 2, 3, 0, 1, 5, 6, 6, 2, 1, 7, 6, 5, 5, 4, 7,
    4, 0, 3, 3, 7, 4, 4, 5, 1, 1, 0, 4, 3, 2, 6, 6, 7, 3,
  };
  cube_ebuf->upload(cube_elements);
  glEnableVertexAttribArray(vattr);
  vao->unbind();
  glDisableVertexAttribArray(vattr);
  cube_vbuf->unbind();
  cube_ebuf->unbind();

  sp->use_this_prog();
  glUniform2f(resolution_unif, static_cast<GLfloat>(s->window_width)
      , static_cast<GLfloat>(s->window_height));

  glm::mat4 projection_mat = glm::perspective(fov, (float)s->window_width
      / (float)s->window_height, 0.1f, 100.0f);
  glUniformMatrix4fv(sp->bind_uniform("projection"), 1, GL_FALSE
      , glm::value_ptr(projection_mat));
  sp->dont_use_this_prog();

  cam = new camera(2, 1, -10);
}

static void load(screen *s) {
  graphics_load(s);
}

static void key_event(char key, bool down) {
  static bool forward = false, backward = false, left = false, right = false;
  switch (key) {
    case 'w':
      forward = down;
      break;
    case 's':
      backward = down;
      break;
    case 'd':
      right = down;
      break;
    case 'a':
      left = down;
      break;
    default:
      break;
  }
  if (forward == backward)
    move = 0;
  else if (forward)
    move = 1;
  else if (backward)
    move = -1;
  if (right == left)
    strafe = 0;
  else if (right)
    strafe = 1;
  else if (left)
    strafe = -1;
}

static void mouse_motion_event(float xrel, float yrel, int x, int y) {
  cam->update_view_angles(xrel, yrel);
}

static void mouse_button_event(int button, bool down) {
}

static void update(double dt, double t, screen *s) {
  sp->use_this_prog();
  glUniform1f(time_unif, static_cast<GLfloat>(t));
  sp->dont_use_this_prog();
  cam->update_position(dt, move, strafe);
}

static void draw_cube(glm::vec3 pos, glm::vec2 size, float rotation
    , glm::vec3 color) {
  glm::mat4 model;
  glm::vec3 rotation_pivot(0, 0, 0);
  model = glm::translate(model, pos);
  model = glm::rotate(model, to_radians(rotation), glm::vec3(0, 0, 1));
  model = glm::scale(model, glm::vec3(size, 1));
  model = glm::translate(model, -rotation_pivot);

  glUniformMatrix4fv(modelmat_unif, 1, GL_FALSE, glm::value_ptr(model));
  glUniform3f(color_unif, color.x, color.y, color.z);
  vao->bind();
  int ebuf_size;
  glGetBufferParameteriv(GL_ELEMENT_ARRAY_BUFFER, GL_BUFFER_SIZE, &ebuf_size);
  glDrawElements(GL_TRIANGLES, ebuf_size
      / static_cast<GLsizei>(sizeof(GLushort)), GL_UNSIGNED_SHORT, 0);
  vao->unbind();
}

static void draw(double alpha) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  sp->use_this_prog();
  glUniformMatrix4fv(viewmat_unif, 1, GL_FALSE
      , glm::value_ptr(cam->compute_view_mat()));

  draw_cube(glm::vec3(0, 0, 0), glm::vec2(1, 1), 90, glm::vec3(1, 0, 1));
  sp->dont_use_this_prog();
}

static void cleanup() {
  delete vs;
  delete fs;
  delete sp;
  delete cube_vbuf;
  delete vao;
}

int main() {
  try {
    screen s("woof", 700, 525);

    s.mainloop(load, key_event, mouse_motion_event, mouse_button_event, update
        , draw, cleanup);
  } catch (const std::exception &e) {
    die("exception exit: %s", e.what());
  } catch (...) {
    die("unknown exception exit");
  }

  return 0;
}

