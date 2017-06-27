#include "bsp.hh"
#include "camera.hh"
#include "frustum.hh"
#include "ogl.hh"
#include "screen.hh"
#include "shaders.hh"
#include "utils.hh"

static float fov = 60, screen_aspect_ratio;
static int move, strafe;
static cube_drawer *cd;
static sphere_drawer *sd;
static axis_drawer *ad;
static entity *cube_ent;
static camera *cam;
static bsp *b;
static bool wireframe = false, noclip = true, update_frustum_culling = true;
static int movement_switch = 0;
static frustum f;

static void graphics_load(screen *s) {
  screen_aspect_ratio = static_cast<float>(s->window_width)
      / static_cast<float>(s->window_height);

  cd = new cube_drawer;
  sd = new sphere_drawer;
  ad = new axis_drawer;

  cube_ent = new entity(glm::vec3(1, 2, 5));
  cam = new camera(cube_ent);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

  glClearColor(0.051f, 0.051f, 0.051f, 1);

  b = new bsp("mapz/ztn3tourney1.bsp", 32.f, 10);

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

static void mouse_button_event(int button, bool down) {
}

static void update(double dt, double t, screen *s) {
  cube_ent->update(dt, t);
  if (noclip) {
    cam->update_position(dt, move, strafe);
    cam->vel = glm::vec3(0);
  } else {
    auto accelerate = [](glm::vec3 prev_velocity, glm::vec3 wish_dir, float dtf) {
      float wish_speed = 4.f, accel = 10.f;
      if (movement_switch == 0) {
        // Q3 movement
        puts("Q3 movement");
        float proj_vel = glm::dot(prev_velocity, wish_dir)
          , accel_speed = accel * dtf * wish_speed;
        if (accel_speed > wish_speed - proj_vel)
          accel_speed = wish_speed - proj_vel;
        return prev_velocity + accel_speed * wish_dir;
      } else if (movement_switch == 1) {
        // Article movement
        puts("Article movement");
        float proj_vel = glm::dot(prev_velocity, wish_dir)
          , accel_speed = accel * dtf;
        if (accel_speed > wish_speed - proj_vel)
          accel_speed = wish_speed - proj_vel;
        return prev_velocity + accel_speed * wish_dir;
      } else if (movement_switch == 2) {
        // Q3 fixed movement
        // TODO: bugful
        puts("Q3 fixed movement");
        glm::vec3 wish_velocity = wish_dir * wish_speed
          , push_dir = glm::normalize(wish_velocity - prev_velocity);
        float push_len = glm::length(push_dir)
          , can_push = accel * dtf * wish_speed;
        if (can_push > push_len)
          can_push = push_len;
        return prev_velocity + push_dir * can_push;
      } else {
        // Q3 fixed movement article like
        puts("Q3 fixed movement article like");
        glm::vec3 wish_velocity = wish_dir * wish_speed
          , push_dir = glm::normalize(wish_velocity - prev_velocity);
        float push_len = glm::length(push_dir)
          , can_push = accel * dtf;
        if (can_push > push_len)
          can_push = push_len;
        return prev_velocity + push_dir * can_push;
      }
    };

    auto clip_velocity = [](glm::vec3 velocity, glm::vec3 normal) {
      const float overclip = 1.01f;
      float back_off = glm::dot(velocity, normal);
      if (back_off < 0)
        back_off *= overclip;
      else
        back_off /= overclip;
      return velocity - normal * back_off;
    };

    auto slide_move = [clip_velocity](glm::vec3 &position, glm::vec3 old_vel, glm::vec3 &vel, float dtf) {
      const int max_clip_planes = 5, num_bumps = 5;
      glm::vec3 planes[max_clip_planes];
      int num_planes = 0, i;
      planes[num_planes] = glm::normalize(vel);
      ++num_planes;
      float time_left = dtf;
      glm::vec3 end_vel = vel;
      for (int bump_count = 0; bump_count < num_bumps; ++bump_count) {
        // integration step
        glm::vec3 wish_pos = position + (old_vel + vel) * 0.5f * (float)time_left;
        trace_result tr;
        b->trace_sphere(&tr, position, wish_pos, 0.25f);
        if (tr.fraction > 0)
          position = tr.end;
        if (tr.fraction == 1)
          break;
        time_left -= time_left * tr.fraction;
        if (num_planes >= max_clip_planes) {
          vel = glm::vec3(0);
          return;
        }
        for (i = 0; i < num_planes; ++i)
          if (glm::dot(tr.clip_plane_normal, planes[i]) > 0.99f) {
            vel += tr.clip_plane_normal;
            break;
          }
        if (i < num_planes)
          continue;
        planes[num_planes] = tr.clip_plane_normal;
        ++num_planes;
        for (i = 0 ; i < num_planes; ++i) {
          float into = glm::dot(vel, planes[i]);
          if (into >= 0.1f)
            continue;

          // if (-into > cam->impact_speed)
          //   cam->impact_speed = -into;

          glm::vec3 clip_vel = clip_velocity(vel, planes[i])
            , end_clip_vel = clip_velocity(end_vel, planes[i]);

          // see if there is a second plane that the new move enters
          for (int j = 0; j < num_planes; ++j) {
            if (j == i)
              continue;
            if (glm::dot(clip_vel, planes[j]) >= 0.1f)
              continue; // move doesn't interact with the plane

            // try clipping the move to the plane
            clip_vel = clip_velocity(clip_vel, planes[j]);
            end_clip_vel = clip_velocity(end_clip_vel, planes[j]);

            // see if it goes back into the first clip plane
            if (glm::dot(clip_vel, planes[i]) >= 0)
              continue;

            // slide the original velocity along the crease
            glm::vec3 dir = glm::normalize(glm::cross(planes[i], planes[j]));
            float d = glm::dot(dir, vel);
            clip_vel = dir * d;
            d = glm::dot(dir, end_vel);
            end_clip_vel = dir * d;

            // see if there is a third plane the the new move enters
            for (int k = 0; k < num_planes; ++k) {
              if (k == i || k == j)
                continue;
              if (glm::dot(clip_vel, planes[k]) >= 0.1f)
                continue; // move doesn't interact with the plane
              // stop dead at a triple plane interaction
              vel = glm::vec3(0);
            }
          }

          // if we have fixed all interactions, try another move
          vel = clip_vel;
          end_vel = end_clip_vel;
          break;
        }
      }
    };
    /* vec3 old_vel = vel;
       vel += accel * dt;
       pos += (old_vel + vel) * 0.5 * dt;
     */
    glm::vec3 old_vel = cam->vel;
    cam->vel = accelerate(cam->vel, cam->compute_view_dir(), dt);
    slide_move(cam->pos, old_vel, cam->vel, (float)dt);
  }
}

static void draw(double alpha) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (wireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glm::mat4 projection = glm::perspective(glm::radians(fov), screen_aspect_ratio
          , 0.1f, 100.f), view = cam->compute_view_mat();

  if (update_frustum_culling)
    f.extract_planes(projection * view);

  b->draw(cam->pos, projection * view, f);
  cd->draw(projection * view * cube_ent->compute_model_mat());
  ad->draw(projection * view);
}

static void cleanup() {
  delete cam;
  delete b;
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

