#include "bsp.hh"
#include "camera.hh"
#include "ogl.hh"
#include "screen.hh"
#include "shaders.hh"
#include "utils.hh"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

static float fov = 45, screen_aspect_ratio;
static shaderprogram *sp;
static shader *vs, *fs;
static GLint vertex_pos_attr, texture_coord_attr, lightmap_coord_attr /*, vertex_normal_attr */;
// static array_buffer *cube_vbuf;
static GLint /* resolution_unif, time_unif, */ model_mat_unif, view_mat_unif
    , projection_mat_unif /* , object_color_unif, view_pos_unif, light_pos_unif */
    , texture_sampler_unif, lightmap_sampler_unif;
static vertexarray *vao;
static int move, strafe;
static camera *cam;
static bsp *b;

static void graphics_load(screen *s) {
  s->lock_mouse();

  screen_aspect_ratio = static_cast<float>(s->window_width)
      / static_cast<float>(s->window_height);

  cam = new camera(0, 5, 0);

  glEnable(GL_DEPTH_TEST);

  glClearColor(0.051f, 0.051f, 0.051f, 1);

  vs = new shader(shaders::map_vert, GL_VERTEX_SHADER);
  fs = new shader(shaders::map_frag, GL_FRAGMENT_SHADER);
  sp = new shaderprogram(*vs, *fs);
  sp->use_this_prog();

  vertex_pos_attr = sp->bind_attrib("vertex_pos");
  texture_coord_attr = sp->bind_attrib("texture_coord");
  lightmap_coord_attr = sp->bind_attrib("lightmap_coord");
  model_mat_unif = sp->bind_uniform("model");
  view_mat_unif = sp->bind_uniform("view");
  projection_mat_unif = sp->bind_uniform("projection");

  glUniform1i(glGetUniformLocation(sp->id, "texture_sampler"), 0);
  glUniform1i(glGetUniformLocation(sp->id, "lightmap_sampler"), 1);

  b = new bsp("mapz/test1.bsp");
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

  sp->use_this_prog();

  glUniformMatrix4fv(view_mat_unif, 1, GL_FALSE
      , glm::value_ptr(cam->compute_view_mat()));

  glUniformMatrix4fv(projection_mat_unif, 1, GL_FALSE
      , glm::value_ptr(glm::perspective(glm::radians(fov), screen_aspect_ratio
          , 0.1f, 1000.f)));

  glUniformMatrix4fv(model_mat_unif, 1, GL_FALSE, glm::value_ptr(glm::mat4()));

  glEnableVertexAttribArray(vertex_pos_attr);
  glEnableVertexAttribArray(texture_coord_attr);
  glEnableVertexAttribArray(lightmap_coord_attr);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  b->set_visible_faces(cam->pos());

  for (int i = 0; i < b->faces.size(); i++) {
    if (!b->visible_faces[i] || (b->faces[i].type != 1 && b->faces[i].type != 3))
      continue;
    glVertexAttribPointer(vertex_pos_attr, 3, GL_FLOAT, GL_FALSE
        , sizeof(bsp_vertex), &b->vertices[b->faces[i].vertex].position);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, b->texture_ids[b->faces[i].texture]);
    glVertexAttribPointer(texture_coord_attr, 2, GL_FLOAT, GL_FALSE
        , sizeof(bsp_vertex), &b->vertices[b->faces[i].vertex].decal);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, b->lightmap_texture_ids[b->faces[i].lm_index]);
    glVertexAttribPointer(lightmap_coord_attr, 2, GL_FLOAT, GL_FALSE
        , sizeof(bsp_vertex), &b->vertices[b->faces[i].vertex].lightmap);
    glDrawElements(GL_TRIANGLES, b->faces[i].n_meshverts, GL_UNSIGNED_INT
        , &b->meshverts[b->faces[i].meshvert].offset);
  }

  glDisableVertexAttribArray(vertex_pos_attr);
  glDisableVertexAttribArray(texture_coord_attr);
  glDisableVertexAttribArray(lightmap_coord_attr);

  sp->dont_use_this_prog();
}

static void cleanup() {
  delete vs;
  delete fs;
  delete sp;
  // delete cube_vbuf;
  delete vao;
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

