#include "camera.hh"
#include "ogl.hh"
#include "screen.hh"
#include "utils.hh"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

static shaderprogram *sp;
static GLuint vertex_pos_attr, vertex_normal_attr;
static array_buffer *cube_vbuf;
static GLint resolution_unif, time_unif, model_mat_unif, view_mat_unif
    , projection_mat_unif, object_color_unif, view_pos_unif, light_pos_unif;
static shader *vs, *fs;
static vertexarray *vao;
static camera *cam;
static float fov = 106.26f, screen_aspect_ratio;
static int move, strafe;

static void graphics_load(screen *s) {
  s->lock_mouse();

  glEnable(GL_DEPTH_TEST);

  glClearColor(0.055f, 0.055f, 0.055f, 1);

  const char *vsrc = _glsl(
    attribute vec3 vertex_pos;
    attribute vec3 vertex_normal;
    uniform mat4 model;
    uniform mat4 view;
    uniform mat4 projection;
    varying vec3 frag_normal;
    varying vec3 frag_pos;
    void main() {
      frag_pos = vec3(model * vec4(vertex_pos, 1.0f));
      frag_normal = normalize(mat3(transpose(inverse(model))) * vertex_normal);
      gl_Position = projection * view * model * vec4(vertex_pos, 1.0);
    }
  );
  const char *fsrc = _glsl(
    uniform vec2 iResolution;
    uniform float iGlobalTime;
    uniform vec3 light_pos;
    uniform vec3 view_pos;
    uniform vec3 object_color;
    varying vec3 frag_normal;
    varying vec3 frag_pos;
    void main() {
      vec3 light_color = vec3(1, 1, 1);
      float ambient_strength = 0.1f;
      vec3 ambient = ambient_strength * light_color;

      vec3 light_dir = normalize(light_pos - frag_pos);
      float diff = max(dot(frag_normal, light_dir), 0.0);
      vec3 diffuse = diff * light_color;

      float specular_strength = 0.5f;
      vec3 view_dir = normalize(view_pos - frag_pos);
      vec3 reflect_dir = reflect(-light_dir, frag_normal);
      float spec = pow(max(dot(view_dir, reflect_dir), 0.0), 32);
      vec3 specular = specular_strength * spec * light_color;

      vec3 result = (ambient + diffuse + specular) * object_color;
      gl_FragColor = vec4(result, 1.0);
    }
  );

  vs = new shader(vsrc, GL_VERTEX_SHADER);
  fs = new shader(fsrc, GL_FRAGMENT_SHADER);
  sp = new shaderprogram(*vs, *fs);

  vertex_pos_attr = sp->bind_attrib("vertex_pos");
  vertex_normal_attr = sp->bind_attrib("vertex_normal");
  resolution_unif = sp->bind_uniform("iResolution");
  time_unif = sp->bind_uniform("iGlobalTime");
  model_mat_unif = sp->bind_uniform("model");
  view_mat_unif = sp->bind_uniform("view");
  projection_mat_unif = sp->bind_uniform("projection");
  object_color_unif = sp->bind_uniform("object_color");
  view_pos_unif = sp->bind_uniform("view_pos");
  light_pos_unif = sp->bind_uniform("light_pos");

  vao = new vertexarray;
  cube_vbuf = new array_buffer;
  vao->bind();
  cube_vbuf->bind();
  const std::vector<float> cube_verts = {
    -.5f, -.5f, -.5f,  0,  0, -1,  .5f, -.5f, -.5f,  0,  0, -1,
     .5f,  .5f, -.5f,  0,  0, -1,  .5f,  .5f, -.5f,  0,  0, -1,
    -.5f,  .5f, -.5f,  0,  0, -1, -.5f, -.5f, -.5f,  0,  0, -1,
    -.5f, -.5f,  .5f,  0,  0,  1,  .5f, -.5f,  .5f,  0,  0,  1,
     .5f,  .5f,  .5f,  0,  0,  1,  .5f,  .5f,  .5f,  0,  0,  1,
    -.5f,  .5f,  .5f,  0,  0,  1, -.5f, -.5f,  .5f,  0,  0,  1,
    -.5f,  .5f,  .5f, -1,  0,  0, -.5f,  .5f, -.5f, -1,  0,  0,
    -.5f, -.5f, -.5f, -1,  0,  0, -.5f, -.5f, -.5f, -1,  0,  0,
    -.5f, -.5f,  .5f, -1,  0,  0, -.5f,  .5f,  .5f, -1,  0,  0,
     .5f,  .5f,  .5f,  1,  0,  0,  .5f,  .5f, -.5f,  1,  0,  0,
     .5f, -.5f, -.5f,  1,  0,  0,  .5f, -.5f, -.5f,  1,  0,  0,
     .5f, -.5f,  .5f,  1,  0,  0,  .5f,  .5f,  .5f,  1,  0,  0,
    -.5f, -.5f, -.5f,  0, -1,  0,  .5f, -.5f, -.5f,  0, -1,  0,
     .5f, -.5f,  .5f,  0, -1,  0,  .5f, -.5f,  .5f,  0, -1,  0,
    -.5f, -.5f,  .5f,  0, -1,  0, -.5f, -.5f, -.5f,  0, -1,  0,
    -.5f,  .5f, -.5f,  0,  1,  0,  .5f,  .5f, -.5f,  0,  1,  0,
     .5f,  .5f,  .5f,  0,  1,  0,  .5f,  .5f,  .5f,  0,  1,  0,
    -.5f,  .5f,  .5f,  0,  1,  0, -.5f,  .5f, -.5f,  0,  1,  0
  };
  cube_vbuf->upload(cube_verts);
  glVertexAttribPointer(vertex_pos_attr, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), 0);
  glEnableVertexAttribArray(vertex_pos_attr);
  glVertexAttribPointer(vertex_normal_attr, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (GLvoid*)(3 * sizeof(GLfloat)));
  glEnableVertexAttribArray(vertex_normal_attr);
  vao->unbind();
  glDisableVertexAttribArray(vertex_pos_attr);
  glDisableVertexAttribArray(vertex_normal_attr);
  cube_vbuf->unbind();

  sp->use_this_prog();
  glUniform2f(resolution_unif, static_cast<GLfloat>(s->window_width)
      , static_cast<GLfloat>(s->window_height));
  sp->dont_use_this_prog();

  screen_aspect_ratio = static_cast<float>(s->window_width)
      / static_cast<float>(s->window_height);

  cam = new camera(-10, 0, 0);
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
  sp->use_this_prog();
  glUniform1f(time_unif, static_cast<GLfloat>(t));
  sp->dont_use_this_prog();
  cam->update_position(dt, move, strafe);
}

static void draw_cube(glm::vec3 pos, glm::vec3 size, glm::vec3 color) {
  glm::vec3 pivot(0, 0, 0);
  glm::mat4 model = glm::translate(glm::mat4(), -pos);
  model = glm::scale(model, size);
  model = glm::translate(model, -pivot);

  glUniformMatrix4fv(model_mat_unif, 1, GL_FALSE, glm::value_ptr(model));
  glUniform3f(object_color_unif, color.x, color.y, color.z);
  vao->bind();
  glDrawArrays(GL_TRIANGLES, 0, 36);
  vao->unbind();
}

static void draw(double alpha) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  sp->use_this_prog();
  glUniformMatrix4fv(view_mat_unif, 1, GL_FALSE
      , glm::value_ptr(cam->compute_view_mat()));

  glUniformMatrix4fv(projection_mat_unif, 1, GL_FALSE
      , glm::value_ptr(glm::perspective(fov, screen_aspect_ratio, 0.1f, 100.0f)));

  glUniform3f(view_pos_unif, cam->pos_x(), cam->pos_y(), cam->pos_z());
  glUniform3f(light_pos_unif, 0, 1.5f, 0);

  draw_cube(glm::vec3(0, 0, 0), glm::vec3(1, 1, 1), glm::vec3(1, 0, 1));

  draw_cube(glm::vec3(0, 1.5, 0), glm::vec3(0.03), glm::vec3(1));

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

