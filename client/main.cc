#include "bsp.hh"
#include "camera.hh"
#include "ogl.hh"
#include "screen.hh"
#include "shaders.hh"
#include "utils.hh"

static float fov = 45, screen_aspect_ratio;
// static GLint vertex_normal_attr;
// static array_buffer *cube_vbuf;
// static GLint resolution_unif, time_unif, object_color_unif, view_pos_unif
//     , light_pos_unif;
// static vertex_array *vao;
static int move, strafe;
static camera *cam;
static bsp *b;
static bool wireframe = false;

static void graphics_load(screen *s) {
  s->lock_mouse();

  screen_aspect_ratio = static_cast<float>(s->window_width)
      / static_cast<float>(s->window_height);

  cam = new camera(0, 5, 0);

  glEnable(GL_DEPTH_TEST);

  glClearColor(0.051f, 0.051f, 0.051f, 1);

  b = new bsp("mapz/ztn3tourney1.bsp", 32.f, 10);
}

static void load(screen *s) {
  graphics_load(s);
}

static void key_event(char key, bool down) {
  static bool forward = false, backward = false, left = false, right = false;
  switch (key) {
    case 'w': forward  = down; break;
    case 's': backward = down; break;
    case 'd': right    = down; break;
    case 'a': left     = down; break;
    case 'f': if (down) wireframe = !wireframe; break;
    default: break;
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
#if 0
  sp->use_this_prog();
  glUniform1f(time_unif, static_cast<GLfloat>(t));
  sp->dont_use_this_prog();
#endif
  cam->update_position(dt, move, strafe);
  // light_pos.x = cos(t) * 2;
  // light_pos.z = sin(t) * 2;
}

#if 0
static void draw_cube(glm::vec3 pos, glm::vec3 size, glm::vec3 color) {
  glm::vec3 pivot(0, 0, 0);
  glm::mat4 model = glm::translate(glm::mat4(), pos);
  model = glm::scale(model, size);
  model = glm::translate(model, -pivot);

  glUniformMatrix4fv(model_mat_unif, 1, GL_FALSE, glm::value_ptr(model));
  // glUniform3f(object_color_unif, color.x, color.y, color.z);
  vao->bind();
  glDrawArrays(GL_TRIANGLES, 0, 36);
  vao->unbind();
}
#endif

static void draw(double alpha) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (wireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glm::mat4 projection = glm::perspective(glm::radians(fov), screen_aspect_ratio
          , 0.1f, 1000.f)
    , view = cam->compute_view_mat(), model = glm::mat4()
    , mvp = projection * view * model;

  b->render(cam->pos(), mvp);
}

static void cleanup() {
  // delete cube_vbuf;
  // delete vao;
  delete cam;
  delete b;
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

