#include <math.h>

#include "madgwick.h"
#include "attitude.h"

madgwick::madgwick()
{
  this->beta = 0.00;
  this->q[0] = 1.0;
  this->q[1] = 0.0;
  this->q[2] = 0.0;
  this->q[3] = 0.0;
}

madgwick::madgwick(double beta)
{
  this->beta = beta;
}

void madgwick::set_beta(double beta)
{
  this->beta = beta;
}

// Set initial quaternion
void madgwick::set_q0(const double q0[4])
{
  this->q[0] = q0[0];
  this->q[1] = q0[1];
  this->q[2] = q0[2];
  this->q[3] = q0[3];
}

void madgwick::update(const double g[3], const double a[3], const double dt)
{
  // q_dot using angular rate
  double dq[4] = {0.0, 0.0, 0.0};
  const double w[3] = {g[0], g[1], g[2]};
  quat_rate(q, w, dq);

  // Auxiliary variables for gradient descent
  const double _4q0 = 4.0f * q[0]; const double _2q0 = 2.0f * q[0]; const double q0q0 = q[0] * q[0];
  const double _4q1 = 4.0f * q[1]; const double _2q1 = 2.0f * q[1]; const double q1q1 = q[1] * q[1];
  const double _4q2 = 4.0f * q[2]; const double _2q2 = 2.0f * q[2]; const double q2q2 = q[2] * q[2];
  const double _8q1 = 8.0f * q[1]; const double _2q3 = 2.0f * q[3]; const double q3q3 = q[3] * q[3];
  const double _8q2 = 8.0f * q[2];

  // Gradient decent algorithm corrective step
  double dq_hat[4];
  dq_hat[0] = _4q0 * q2q2 + _2q2 * a[0] + _4q0 * q1q1 - _2q1 * a[1];
  dq_hat[1] = _4q1 * q3q3 - _2q3 * a[0] + 4.0f * q0q0 * q[1] - _2q0 * a[1] - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * a[2];
  dq_hat[2] = 4.0f * q0q0 * q[2] + _2q0 * a[0] + _4q2 * q3q3 - _2q3 * a[1] - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * a[2];
  dq_hat[3] = 4.0f * q1q1 * q[3] - _2q1 * a[0] + 4.0f * q2q2 * q[3] - _2q2 * a[1];

  // Estimated rate of quaternion
  // f_grad / norm(f_grad)
  quat_normalize(dq_hat);

  for (int i = 0; i < 4; i++)
  {
    dq[i] -= beta * dq_hat[i]; // Feedback correction
    q[i] += dq[i] * dt;        // Numerical integration
  }

  quat_normalize(q);
}
