// Least squares based magnetometer calibration.
// rms (2025-05-16 10:31:45 AM)

#ifndef _MAGCAL_H_
#define _MAGCAL_H_

#include <vector>

#include "Eigen"
#include "data_types.h"

using namespace Eigen;

typedef struct
{
  Matrix3d Ainv; // Soft-iron correction matrix
  Vector3d b;    // Hard-iron bias vector
  double r;      // Magnetic field magnitude
  double err;    // Fit error
  bool is_pd;    // Is positive definite
} magcal_t;

namespace magcal
{
  magcal_t compute(const VectorXd &x, const VectorXd &y, const VectorXd &z);
  void calibrate(const magcal_t &magcal,
                 const VectorXd &x, const VectorXd &y, const VectorXd &z,
                 VectorXd &x_cal, VectorXd &y_cal, VectorXd &z_cal);

  magcal_t compute(const std::vector<neo_vec3> &m);
  std::vector<neo_vec3> magcal(const std::vector<neo_vec3> &m);
  std::vector<neo_vec3> calibrate(const magcal_t params, const std::vector<neo_vec3> &m);
};

#endif // magcal.h
