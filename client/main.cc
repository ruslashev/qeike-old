#include "camera.hh"
#include "../engine/frustum.hh"
#include "../engine/ogl.hh"
#include "../engine/screen.hh"
#include "../engine/shaders.hh"
#include "../engine/utils.hh"

static float fov = 60, screen_aspect_ratio;
static int move, strafe;
// static axis_drawer *ad;
static entity *cube_ent;
static camera *cam;
static bool wireframe = false, noclip = true, update_frustum_culling = true;
static int movement_switch = 0;
static frustum f;

static void graphics_load(screen *s) {
  screen_aspect_ratio = static_cast<float>(s->window_width)
      / static_cast<float>(s->window_height);

  // ad = new axis_drawer;

  cube_ent = new entity(glm::vec3(1, 2, 5));
  cam = new camera(cube_ent);

  glEnable(GL_DEPTH_TEST);
  // glEnable(GL_CULL_FACE);

  glClearColor(0.051f, 0.051f, 0.051f, 1);

  // b = new bsp("mapz/q3/q3dm6.bsp", 15);

  s->lock_mouse();
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
    case 'x': if (down) noclip = !noclip; break;
    case 'z': if (down) update_frustum_culling = !update_frustum_culling; break;
    case 'c': if (down) movement_switch = (movement_switch + 1) % 4; break;
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

static void mouse_button_event(int button, bool down, int x, int y) {
}

static void update(double dt, double t, screen *s) {
  cube_ent->update(dt, t);
  if (noclip) {
    cam->update_position(dt, move, strafe);
    cam->vel = glm::vec3(0);
  } else {
    auto accelerate = [](glm::vec3 prev_velocity, glm::vec3 wish_dir, float dtf) {
      float wish_speed = 200.f, accel = 1.f;
#if 0
      if (movement_switch == 0) {
        // puts("Q3 movement");
        float proj_vel = glm::dot(prev_velocity, wish_dir)
          , accel_speed = accel * dtf * wish_speed;
        if (accel_speed > wish_speed - proj_vel)
          accel_speed = wish_speed - proj_vel;
        return prev_velocity + accel_speed * wish_dir;
      } else if (movement_switch == 1) {
        // puts("Article movement");
        float proj_vel = glm::dot(prev_velocity, wish_dir)
          , accel_speed = accel * dtf;
        if (accel_speed > wish_speed - proj_vel)
          accel_speed = wish_speed - proj_vel;
        return prev_velocity + accel_speed * wish_dir;
      } else if (movement_switch == 2) {
#endif
        // puts("Q3 fixed movement");
        glm::vec3 wish_velocity = wish_dir * wish_speed
          , push_vec = wish_velocity - prev_velocity;
        if (glm::length(push_vec) < 0.1f)
          return prev_velocity;
        glm::vec3 push_dir = glm::normalize(push_vec);
        float push_len = glm::length(push_dir)
          , can_push = accel * dtf * wish_speed;
        if (can_push > push_len)
          can_push = push_len;
        return prev_velocity + push_dir * can_push;
#if 0
      } else
        return prev_velocity;
#endif
    };

    /* vec3 old_vel = vel;
       vel += accel * dt;
       pos += (old_vel + vel) * 0.5 * dt;
     */
    glm::vec3 old_vel = cam->vel;
    cam->vel = accelerate(cam->vel, cam->compute_view_dir(), dt);
    cam->pos += (old_vel + cam->vel) * 0.5f * (float)dt;
  }
}

static void draw(double alpha) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (wireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glm::mat4 projection = glm::perspective(glm::radians(fov), screen_aspect_ratio
          , 0.1f, 10000.f), view = cam->compute_view_mat();

  if (update_frustum_culling) {
    f.extract_planes(projection * view);
    f.position = cam->pos;
  }

  // ad->draw(projection * view);
}

static void cleanup() {
  delete cam;
  delete cube_ent;
}

int main() {
  try {
    screen s("qeike", 700, 525);

    s.mainloop(load, key_event, mouse_motion_event, mouse_button_event, update
        , draw, cleanup);
  } catch (const std::exception &e) {
    die("exception exit: %s", e.what());
  } catch (...) {
    die("unknown exception exit");
  }

  return 0;
}

