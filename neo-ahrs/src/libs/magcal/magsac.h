/**
 * magcal() with RANSAC
 *
 * This is a wrapper of RANSAC which uses magcal() from magcal.cpp/h as is for calibration.
**/

#ifndef _MAGSAC_H_
#define _MAGSAC_H_

#include <vector>
#include <random>
#include <memory>

#include "Eigen"
#include "magcal.h"
#include "data_types.h"

namespace magsac
{
  struct magsac_opts
  {
    magsac_opts()
        : min_num_iterations(100),
          max_num_iterations(10000),
          success_probability(0.99),
          inlier_threshold(0.1),
          min_inlier_percent(0.5),
          random_seed(0) {}

    uint32_t min_num_iterations; // Maximum number of RANSAC iterations
    uint32_t max_num_iterations; // Minimum number of RANSAC iterations
    double success_probability;  // Success probability of RANSAC
    double inlier_threshold;     // Error tolerance for inliers
    double min_inlier_percent;   // Minimum fraction of inliers to accept model
    unsigned int random_seed;    // Controls the random sampling
  };

  struct magsac_stats
  {
    bool success;
    double inlier_ratio;
    int best_num_inliers;
    uint32_t num_iterations;
    double best_model_score;
    std::vector<int> inlier_indices;
  };

  class magsac
  {
  protected:
    const Eigen::MatrixXd &data_;
    const int min_sample_size_;

  public:
    magsac(const Eigen::MatrixXd &data, int min_sample_size = 10)
        : data_(data), min_sample_size_(min_sample_size) {}

    int min_sample_size() const { return min_sample_size_; }

    // Ransac methods
    bool solve_magcal_params(const std::vector<int> &sample_indices, magcal_t *model) const;
    double evaluate_magcal_on_point(const magcal_t &model, int point_index) const;

    // Try to refit all inliers with standard magcal
    bool refit_inliers(const std::vector<int> &inlier_indices, magcal_t *model) const
    {
      if (inlier_indices.size() < min_sample_size_)
      {
        return false;
      }

      Eigen::VectorXd x(inlier_indices.size()), y(inlier_indices.size()), z(inlier_indices.size());
      for (size_t i = 0; i < inlier_indices.size(); ++i)
      {
        int idx = inlier_indices[i];
        x(i) = data_(idx, 0);
        y(i) = data_(idx, 1);
        z(i) = data_(idx, 2);
      }

      try
      {
        *model = magcal::compute(x, y, z);
        return true;
      }
      // todo: Throw error from magcal::compute()
      catch (...)
      {
        return false;
      }
    }
  };

  // RANSAC calibration function
  magcal_t compute_ransac(const std::vector<neo_vec3> &measurements,
                          const magsac_opts &options = magsac_opts(),
                          magsac_stats *stats = nullptr);

  // RANSAC calibration function overload for Eigen vectors
  magcal_t compute_ransac(const Eigen::VectorXd &x, const Eigen::VectorXd &y, const Eigen::VectorXd &z,
                          const magsac_opts &options = magsac_opts(),
                          magsac_stats *stats = nullptr);

} // namespace magcal

#endif // magsac.h
