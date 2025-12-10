#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm>

#include "magsac.h"

namespace magsac
{
  // Calculate required number of RANSAC iterations
  uint32_t compute_iteration_num(double inlier_ratio, double failure_probability,
                                 int min_sample_size, uint32_t min_iterations,
                                 uint32_t max_iterations)
  {
    if (inlier_ratio <= 0.0)
      return max_iterations;

    double sample_failure_prob = 1.0 - std::pow(inlier_ratio, min_sample_size);
    if (sample_failure_prob <= std::numeric_limits<double>::epsilon())
    {
      return min_iterations;
    }

    double num_iter = std::log(failure_probability) / std::log(sample_failure_prob);
    num_iter = std::max(static_cast<double>(min_iterations), num_iter);
    num_iter = std::min(static_cast<double>(max_iterations), num_iter);

    return static_cast<uint32_t>(std::ceil(num_iter));
  }

  // Random sampling without replacement
  void sample_randomly(int total_points, int sample_size, std::mt19937 *rng,
                       std::vector<int> *sample)
  {
    sample->resize(sample_size);
    if (sample_size >= total_points)
    {
      for (int i = 0; i < total_points; ++i)
      {
        (*sample)[i] = i;
      }
      return;
    }

    std::vector<int> indices(total_points);
    for (int i = 0; i < total_points; ++i)
    {
      indices[i] = i;
    }

    for (int i = 0; i < sample_size; ++i)
    {
      std::uniform_int_distribution<int> dist(i, total_points - 1);
      int j = dist(*rng);
      std::swap(indices[i], indices[j]);
      (*sample)[i] = indices[i];
    }
  }

  // Perform magcal for the sampled indices and return the valid (positive definite) model.
  bool magsac::solve_magcal_params(const std::vector<int> &sample_indices, magcal_t *model) const
  {
    if (sample_indices.size() < static_cast<size_t>(min_sample_size_))
    {
      return false;
    }

    try
    {
      Eigen::VectorXd x(sample_indices.size()), y(sample_indices.size()), z(sample_indices.size());
      for (size_t i = 0; i < sample_indices.size(); ++i)
      {
        int idx = sample_indices[i];
        x(i) = data_(idx, 0);
        y(i) = data_(idx, 1);
        z(i) = data_(idx, 2);
      }

      *model = magcal::compute(x, y, z);
      return model->is_pd;
    }
    catch (...)
    {
      return false;
    }
  }

  // Perform magcal for a sample and return the residue.
  double magsac::evaluate_magcal_on_point(const magcal_t &model, int point_index) const
  {
    // Compute corrected magnetic field
    Eigen::Vector3d m_raw(data_(point_index, 0), data_(point_index, 1), data_(point_index, 2));
    Eigen::Vector3d m_corr = model.Ainv * (m_raw - model.b);

    // Error is the deviation from ideal sphere radius
    double mag = m_corr.norm();
    return std::abs(mag - model.r);
  }

  magcal_t compute_ransac(const std::vector<neo_vec3> &measurements,
                          const magsac_opts &options,
                          magsac_stats *statistics)
  {
    // Convert to Eigen matrix
    Eigen::MatrixXd data(measurements.size(), 3);
    for (size_t i = 0; i < measurements.size(); ++i)
    {
      data(i, 0) = measurements[i].x;
      data(i, 1) = measurements[i].y;
      data(i, 2) = measurements[i].z;
    }

    return compute_ransac(data.col(0), data.col(1), data.col(2), options, statistics);
  }

  magcal_t compute_ransac(const Eigen::VectorXd &x, const Eigen::VectorXd &y, const Eigen::VectorXd &z,
                          const magsac_opts &options,
                          magsac_stats *statistics)
  {
    // Initialize statistics
    magsac_stats stats;
    stats.num_iterations = 0;
    stats.best_num_inliers = 0;
    stats.best_model_score = std::numeric_limits<double>::max();
    stats.inlier_ratio = 0.0;
    stats.success = false;

    // Prepare data matrix
    const int N = x.size();
    Eigen::MatrixXd data(N, 3);
    data.col(0) = x;
    data.col(1) = y;
    data.col(2) = z;

    // Sanity checks
    const int min_sample_size = 10;

    if (N < min_sample_size)
    {
      std::cerr << "RANSAC: Not enough data points (" << N << " < " << min_sample_size << ")" << std::endl;

      if (statistics)
      {
        *statistics = stats;
      }

      // Fallback to standard calibration
      return magcal::compute(x, y, z);
    }

    // Initialize RANSAC solver and RNG
    magsac solver(data, min_sample_size);
    std::mt19937 rng(options.random_seed);

    // Best model tracking
    magcal_t best_model;
    double best_score = std::numeric_limits<double>::max();
    std::vector<int> best_inliers;
    int max_inliers = 0;

    // RANSAC main loop
    uint32_t max_iterations = options.max_num_iterations;

    for (stats.num_iterations = 0; stats.num_iterations < max_iterations; ++stats.num_iterations)
    {
      // Random sampling
      std::vector<int> sample;
      sample_randomly(N, min_sample_size, &rng, &sample);

      // Fit model to sample
      magcal_t model;
      if (!solver.solve_magcal_params(sample, &model))
      {
        continue;
      }

      // Evaluate model on all points
      int num_inliers = 0;
      double total_error = 0.0;
      std::vector<int> inlier_indices;

      for (int i = 0; i < N; ++i)
      {
        double error = solver.evaluate_magcal_on_point(model, i);

        if (error < options.inlier_threshold)
        {
          num_inliers++;
          total_error += error;
          inlier_indices.push_back(i);
        }
      }

      // Check if this is a good model
      double inlier_ratio = static_cast<double>(num_inliers) / N;
      double model_score = total_error / std::max(1, num_inliers); // Average error

      bool good_model = (num_inliers >= static_cast<int>((options.min_inlier_percent / 100.0) * N)) &&
                        (num_inliers > max_inliers ||
                         (num_inliers == max_inliers && model_score < best_score));

      if (good_model)
      {
        // Refit using all inliers for better accuracy
        magcal_t refined_model;
        if (solver.refit_inliers(inlier_indices, &refined_model))
        {
          // Re-evaluate refined model
          int refined_inliers = 0;
          double refined_total_error = 0.0;
          std::vector<int> refined_inlier_indices;

          for (int i = 0; i < N; ++i)
          {
            double error = solver.evaluate_magcal_on_point(refined_model, i);
            if (error < options.inlier_threshold)
            {
              refined_inliers++;
              refined_total_error += error;
              refined_inlier_indices.push_back(i);
            }
          }

          double refined_score = refined_total_error / std::max(1, refined_inliers);

          if (refined_inliers >= num_inliers && refined_score <= model_score)
          {
            model = refined_model;
            num_inliers = refined_inliers;
            model_score = refined_score;
            inlier_indices = refined_inlier_indices;
          }
        }

        // Update best model
        if (num_inliers > max_inliers || (num_inliers == max_inliers && model_score < best_score))
        {
          best_model = model;
          best_score = model_score;
          best_inliers = inlier_indices;
          max_inliers = num_inliers;

          // Update iteration count based on inlier ratio
          stats.inlier_ratio = static_cast<double>(num_inliers) / N;
          max_iterations = compute_iteration_num(
              stats.inlier_ratio, 1.0 - options.success_probability,
              min_sample_size, options.min_num_iterations, options.max_num_iterations);
        }
      }
    }

    // Final refinement with all inliers
    if (max_inliers > 0)
    {
      solver.refit_inliers(best_inliers, &best_model);

      // Final evaluation
      stats.best_num_inliers = 0;
      double final_total_error = 0.0;
      stats.inlier_indices.clear();

      for (int i = 0; i < N; ++i)
      {
        double error = solver.evaluate_magcal_on_point(best_model, i);
        if (error < options.inlier_threshold)
        {
          stats.best_num_inliers++;
          final_total_error += error;
          stats.inlier_indices.push_back(i);
        }
      }

      stats.best_model_score = final_total_error / std::max(1, stats.best_num_inliers);
      stats.inlier_ratio = static_cast<double>(stats.best_num_inliers) / N;
      stats.success = true;

      std::cout << "RANSAC: " << stats.best_num_inliers << " inliers ("
                << (100.0 * stats.inlier_ratio) << "%), RMSE = " << stats.best_model_score << std::endl;
    }
    else
    {
      // Fallback to standard calibration
      std::cerr << "RANSAC: No good model found. Using standard calibration." << std::endl;
      best_model = magcal::compute(x, y, z);
      stats.success = false;
    }

    if (statistics)
      *statistics = stats;
    return best_model;
  }

} // namespace magsac
