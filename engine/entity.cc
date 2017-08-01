#include "entity.hh"
#include "math.hh"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/quaternion.hpp>

physics_state::physics_state()
  : position(glm::vec3(0, 0, 0))
  , momentum(glm::vec3(0, 0, 0))
  , orientation(glm::quat())
  , angular_momentum(glm::vec3(0, 0, 0))
  , scale(1)
  , mass(1)
  , inertia_tensor(mass * scale * scale * 1.0f / 6.0f) {
  recalculate();
}

void physics_state::recalculate() {
  // recalculate secondary state values from primary
  velocity = momentum / mass;
  angular_velocity = angular_momentum / inertia_tensor;
  orientation = glm::normalize(orientation);
  spin = glm::quat(0, angular_velocity) * orientation * 0.5f;
}

static physics_state interpolate(const physics_state &a, const physics_state &b
    , float alpha) {
  physics_state state = b;
  state.position = a.position * (1 - alpha) + b.position * alpha;
  state.momentum = a.momentum * (1 - alpha) + b.momentum * alpha;
  state.orientation = slerp(a.orientation, b.orientation, alpha);
  state.angular_momentum = a.angular_momentum * (1 - alpha)
    + b.angular_momentum * alpha;
  state.recalculate();
  return state;
}

entity::entity()
 : entity(glm::vec3(0, 0, 0)) {
}

entity::entity(const glm::vec3 &n_position) {
  _current.position = n_position;
  _current.recalculate();
  _previous = _current;
}

void entity::update(double dt, double t) {
  _previous = _current;
  _integrate(_current, dt, t);
}

void entity::_integrate(physics_state &state, double dt, double t) {
  physics_state_deriv a = _evaluate(state, t, 0, nullptr)
    , b = _evaluate(state, t, dt * 0.5, &a)
    , c = _evaluate(state, t, dt * 0.5, &b)
    , d = _evaluate(state, t, dt, &c);

  float dtf = (float)dt; // temporary
  state.position += 1.0f/6.0f * dtf * (a.velocity + 2.0f * (b.velocity
      + c.velocity) + d.velocity);
  state.momentum += 1.0f/6.0f * dtf * (a.force + 2.0f * (b.force + c.force)
      + d.force);
  state.orientation += 1.0f/6.0f * dtf * (a.spin + 2.0f * (b.spin + c.spin)
      + d.spin);
  state.angular_momentum += 1.0f/6.0f * dtf * (a.torque + 2.0f * (b.torque
      + c.torque) + d.torque);

  state.recalculate();
}

// Evaluate physics_state_deriv at future time t + dt using the specified set
// of derivatives to advance dt seconds from the specified physics state.
physics_state_deriv entity::_evaluate(physics_state state, double t, double dt
    , const physics_state_deriv *derivative) {
  if (derivative) {
    float dtf = (float)dt;
    state.position += derivative->velocity * dtf;
    state.momentum += derivative->force * dtf;
    state.orientation += derivative->spin * dtf;
    state.angular_momentum += derivative->torque * dtf;
    state.recalculate();
  }

  physics_state_deriv output;
  output.velocity = state.velocity;
  output.spin = state.spin;
  _forces(state, t + dt, output.force, output.torque);
  return output;
}

void entity::_forces(const physics_state &state, double t, glm::vec3 &force
    , glm::vec3 &torque) {
  torque.x = 1.0f * (float)glm::sin((float)t * 0.9f + 0.5f);
  torque.y = 1.1f * (float)glm::sin((float)t * 0.5f + 0.4f);
  torque.z = 1.2f * (float)glm::sin((float)t * 0.7f + 0.9f);
  torque -= 0.2f * state.angular_velocity;
}

glm::mat4 entity::compute_model_mat(float alpha) const {
  physics_state state = interpolate(_previous, _current, alpha);
  glm::mat4 id = glm::mat4(), translation = glm::translate(id, state.position)
    , rotation = toMat4(state.orientation);
  return translation * rotation;
}

