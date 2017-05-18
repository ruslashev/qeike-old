#include "ogl.hh"
#include "screen.hh"
#include "utils.hh"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

static shaderprogram *sp;
static GLint vattr;
static array_buffer *model_vbuf;
static GLint resolution_unif, time_unif, modelmat_unif, color_unif;
static shader *vs, *fs;
static vertexarray *vao;

static void graphics_load(screen *s) {
  // s->lock_mouse();

  glClearColor(0.05f, 0.05f, 0.05f, 1);

  const char *vsrc = _glsl(
    attribute vec2 vertex_pos;
    uniform mat4 model;
    uniform mat4 projection;
    void main() {
      gl_Position = projection * model * vec4(vertex_pos, 0.0, 1.0);
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
  modelmat_unif = sp->bind_uniform("model");
  color_unif = sp->bind_uniform("color");

  const std::vector<float> cube_verts = {
    0, 1, 1, 0, 0, 0,
    0, 1, 1, 1, 1, 0
  };
  vao = new vertexarray;
  model_vbuf = new array_buffer;
  vao->bind();
  model_vbuf->bind();
  model_vbuf->upload(cube_verts);
  glVertexAttribPointer(vattr, 2, GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(vattr);
  vao->unbind();
  glDisableVertexAttribArray(vattr);
  model_vbuf->unbind();

  time_unif = sp->bind_uniform("iGlobalTime");

  sp->use_this_prog();
  glUniform2f(resolution_unif, s->window_width, s->window_height);

  glm::mat4 projection_mat = glm::ortho(0.f, (float)s->window_width
      , (float)s->window_height, 0.f, -1.f, 1.f);
  glUniformMatrix4fv(sp->bind_uniform("projection"), 1, GL_FALSE
      , glm::value_ptr(projection_mat));
  sp->dont_use_this_prog();
}

static void load(screen *s) {
  graphics_load(s);
}

static void key_event(char key, bool down) {
}

static void mouse_motion_event(float xrel, float yrel, int x, int y) {
}

static void mouse_button_event(int button, bool down) {
}

static void update(double dt, double t, screen *s) {
  sp->use_this_prog();
  glUniform1f(time_unif, t);
  sp->dont_use_this_prog();
}

static void draw_square(glm::vec2 pos, glm::vec2 size, float rotation
    , glm::vec3 color) {
  glm::mat4 model;
  glm::vec2 rotation_pivot(0.5, 0.5);
  model = glm::translate(model, glm::vec3(pos, 0.f));
  model = glm::rotate(model, to_radians(rotation), glm::vec3(0.f, 0.f, 1.f));
  model = glm::scale(model, glm::vec3(size, 1.f));
  model = glm::translate(model, glm::vec3(-rotation_pivot, 0.f));

  sp->use_this_prog();
  glUniformMatrix4fv(modelmat_unif, 1, GL_FALSE, glm::value_ptr(model));
  glUniform3f(color_unif, color.x, color.y, color.z);
  vao->bind();
  glDrawArrays(GL_TRIANGLES, 0, 6);
  vao->unbind();
  sp->dont_use_this_prog();
}

static void draw(double alpha) {
  glClear(GL_COLOR_BUFFER_BIT);

  draw_square(glm::vec2(100, 200), glm::vec2(10, 10), 90, glm::vec3(1, 0, 1));
}

static void cleanup() {
  delete vs;
  delete fs;
  delete sp;
  delete model_vbuf;
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

