#include <array>
#include <numeric>
#include <algorithm>

#define R2D 57.2957795131
#define D2R 1.0 / 57.2957795131

#include "attitude.h"
#include "madgwick.h"
#include "madmadgwick.h"

// Normalize a 3D vector
neo_vec3 normalize(const neo_vec3& v)
{
  double norm = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  return {v.x / norm, v.y / norm, v.z / norm};
}

std::vector<madgwick_t> madmadgwick(
                        const std::vector<neo_vec3>& a,
                        const std::vector<neo_vec3>& g,
                        const std::vector<double>& t,
                        double beta)
{
  size_t N = a.size();
  std::vector<neo_quat> q_fwd(N);
  std::vector<neo_quat> q_bwd(N);
  std::vector<madgwick_t> output(N);

  /**************************
   * ---- Forward pass ---- *
   **************************/

  madgwick md_fwd;

  double q0_fwd[4] = {1, 0, 0, 0};
  md_fwd.set_q0(q0_fwd);

  for (size_t i = 0; i < N; ++i)
  {
    neo_vec3 a_norm = normalize(a[i]);
    double dt = (i > 0) ? (t[i] - t[i-1]) : 0;

    double g_in[] = {g[i].x, g[i].y, g[i].z};
    double a_in[] = {a_norm.x, a_norm.y, a_norm.z};

    md_fwd.set_beta(beta);
    md_fwd.update(g_in, a_in, dt);

    std::copy(std::begin(md_fwd.q), std::end(md_fwd.q), reinterpret_cast<double*>(&q_fwd[i]));
  }

  /**************************
   * ---- Reverse pass ---- *
   **************************/

  madgwick md_bwd;

  double q0_bwd[4] = {q_fwd[N-1].q0, q_fwd[N-1].q1, q_fwd[N-1].q2, q_fwd[N-1].q3};
  md_bwd.set_q0(q0_bwd);

  for (size_t i = N; i-- > 0;)
  {
    neo_vec3 a_norm = normalize(a[i]);
    double dt = (i < N-1) ? (t[i+1] - t[i]) : 0;

    double g_in[] = {-g[i].x, -g[i].y, -g[i].z};
    double a_in[] = {a_norm.x, a_norm.y, a_norm.z};

    md_bwd.set_beta(beta);
    md_bwd.update(g_in, a_in, dt);

    std::copy(std::begin(md_bwd.q), std::end(md_bwd.q), reinterpret_cast<double*>(&q_bwd[i]));
  }

  // Reverse
  for (size_t i = N; i-- > 0;)
  {
    // Quaternion to Euler angle
    double ypr[3] = {0.0};
    double q[4] = {q_bwd[i].q0, q_bwd[i].q1, q_bwd[i].q2, q_bwd[i].q3};
    quat_to_euler(q, ypr, EULER_ZYX);

    // Flip the sign of yaw and prepare output
    output[i].ypr = {.y = -ypr[0] * R2D, .p = ypr[1] * R2D, .r = ypr[2] * R2D};
    output[i].q = {.q0 = q[0], .q1 = q[1], .q2 = q[2], .q3 = q[3]};
  }

  // LPF the quaternion output
  const double alpha = 0.01;
  for (size_t i = 1; i < N; i++)
  {
    const double q1[4] = {output[i-1].q.q0, output[i-1].q.q1, output[i-1].q.q2, output[i-1].q.q3};
    const double q2[4] = {output[i].q.q0, output[i].q.q1, output[i].q.q2, output[i].q.q3};
    double q[4];

    quat_slerp(q1, q2, alpha, q);
    output[i].q = {.q0 = q[0] , .q1 = q[1], .q2 = q[2], .q3 = q[3]};

    // Flip the sign of yaw and prepare output
    double ypr[3] = {0.0};
    quat_to_euler(q, ypr, EULER_ZYX);
    output[i].ypr = {.y = -ypr[0] * R2D, .p = ypr[1] * R2D, .r = ypr[2] * R2D};
  }

  return output;
}

std::vector<neo_euler> get_euler(const std::vector<madgwick_t> &m)
{
  std::vector<neo_euler> ypr;
  ypr.reserve(m.size());

  for (const auto& mi : m)
  {
    ypr.push_back(mi.ypr);
  }

  return ypr;
}

std::vector<neo_quat> get_quat(const std::vector<madgwick_t> &m)
{
  std::vector<neo_quat> q;
  q.reserve(m.size());

  for (const auto& mi : m)
  {
    q.push_back(mi.q);
  }

  return q;
}
