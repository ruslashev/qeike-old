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
static entity *player;
static camera *cam;
static bsp *b;
static bool wireframe = false, noclip = true, update_frustum_culling = true;
static frustum f;

static void graphics_load(screen *s) {
  screen_aspect_ratio = static_cast<float>(s->window_width)
      / static_cast<float>(s->window_height);

  cd = new cube_drawer;

  player = new entity();
  cam = new camera(player);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_FRONT);

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
  player->update(dt, t);
  // search for CM_BoxTrace and CM_TransformedBoxTrace
  if (noclip) {
    cam->update_position(dt, move, strafe);
    return;
  } else {
    // glm::vec3 old_pos = player->pos;
    // trace_result tr;
    // b->trace_sphere(&tr, old_pos, player->pos, 0.25f);
    // glm::vec3 dot = glm::dot()
    /*
	dot = DotProduct( velocity, trace->plane.normal );
	VectorMA( velocity, -2*dot, trace->plane.normal, le->pos.trDelta );

	VectorScale( le->pos.trDelta, le->bounceFactor, le->pos.trDelta );

	VectorCopy( trace->endpos, le->pos.trBase );
	le->pos.trTime = cg.time;


	// check for stop, making sure that even on low FPS systems it doesn't bobble
	if ( trace->allsolid ||
		( trace->plane.normal[2] > 0 &&
		( le->pos.trDelta[2] < 40 || le->pos.trDelta[2] < -cg.frametime * le->pos.trDelta[2] ) ) ) {
		le->pos.trType = TR_STATIONARY;
     */
  }
}

static void draw(double alpha) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (wireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glm::mat4 projection = glm::perspective(glm::radians(fov), screen_aspect_ratio
          , 0.1f, 100.f)
    , view = cam->compute_view_mat(), model = player->compute_model_mat(alpha)
    , cube_mvp = projection * view * model;

  if (update_frustum_culling)
    f.extract_planes(projection * view);

  b->draw(cam->pos, projection * view, f);
  cd->draw(cube_mvp);
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

