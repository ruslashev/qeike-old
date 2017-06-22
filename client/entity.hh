#pragma once

#include <glm/gtc/quaternion.hpp>

struct physics_state {
  // primary state
  glm::vec3 position;
  glm::vec3 momentum;
  glm::quat orientation;
  glm::vec3 angular_momentum;
  // secondary state
  glm::vec3 velocity;
  glm::quat spin;
  glm::vec3 angular_velocity;
  /// constant state
  float scale;
  float mass;
  float inertia_tensor;

  physics_state();
  void recalculate();
};

struct physics_state_deriv {
  glm::vec3 velocity;
  glm::vec3 force;
  glm::quat spin;
  glm::vec3 torque;
};

class entity {
  physics_state _previous, _current;
  void _integrate(physics_state &state, double dt, double t);
  physics_state_deriv _evaluate(const physics_state &state, double t);
  physics_state_deriv _evaluate(physics_state state, double t, double dt
      , const physics_state_deriv &derivative);
  void _forces(const physics_state &state, double t, glm::vec3 &force
      , glm::vec3 &torque);
public:
  entity();
  entity(const glm::vec3 &n_position);
  void update(double dt, double t);
  glm::mat4 compute_model_mat(float alpha = 1.0f) const;
};

