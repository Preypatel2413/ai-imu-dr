/**
 * @brief Madgwick filter adopted from the original implementation by Madgwick himself.
 *
 * @cite Madgwick - An Efficient Orientation Filter for Inertial and Inertial-Magnetic sensor Arrays
 *       (2010)
 *
 * @author rms
 * @date 2021/08/03
 */

#ifndef _MADGWICK_H_
#define _MADGWICK_H_

#include "math.h"

class madgwick
{
private:
  double beta; // Filter gain

public:
  double q[4]; // Quaternion estimate

  madgwick();
  madgwick(double beta);
  void set_beta(double beta);
  void set_q0(const double q0[4]);
  void update(const double g[3], const double a[3], const double dt);
};

#endif // madgwick.h
