#include "magseg.h"
#include <iostream>
#include <numeric>
#include <deque>
#include <algorithm>

namespace magseg
{

/**
 * Compute windowed rolling standard deviation for each axes of magnetometer measurements and merge
 * them to a common magnetic anomaly metric by computing alpha exponent norm. I am not sure if
 * exponent norm is a thing, but a look at how alpha is used explains the term.
 */
std::vector<double> compute_rolling_std(const std::vector<neo_vec3> &measurements, int window_size, double alpha)
{
  int n = measurements.size();
  std::vector<double> rolling_std(n, 0.0);

  if (n == 0 || window_size <= 0)
    return rolling_std;

  for (int i = 0; i < n; ++i)
  {
    int half_left = window_size / 2;
    int half_right = window_size / 2;

    int start = std::max(0, i - half_left);
    int end = std::min(n - 1, i + half_right - 1);

    int len = end - start + 1;

    double mean_x = 0.0, mean_y = 0.0, mean_z = 0.0;
    for (int j = start; j <= end; ++j)
    {
      mean_x += measurements[j].x;
      mean_y += measurements[j].y;
      mean_z += measurements[j].z;
    }
    mean_x /= len;
    mean_y /= len;
    mean_z /= len;

    double var_x = 0.0, var_y = 0.0, var_z = 0.0;
    for (int j = start; j <= end; ++j)
    {
      var_x += (measurements[j].x - mean_x) * (measurements[j].x - mean_x);
      var_y += (measurements[j].y - mean_y) * (measurements[j].y - mean_y);
      var_z += (measurements[j].z - mean_z) * (measurements[j].z - mean_z);
    }
    if (len > 1)
    {
      var_x /= (len - 1);
      var_y /= (len - 1);
      var_z /= (len - 1);
    }
    else
    {
      var_x = var_y = var_z = 0.0;
    }

    double std_x = std::sqrt(var_x);
    double std_y = std::sqrt(var_y);
    double std_z = std::sqrt(var_z);

    rolling_std[i] = std::sqrt(
        std::pow(std_x, 2 * alpha) +
        std::pow(std_y, 2 * alpha) +
        std::pow(std_z, 2 * alpha));
  }

  return rolling_std;
}

// Computes the relative peak based on the percentage threshold of the mean of the signal.
std::vector<Peak> find_relative_peaks(const std::vector<double> &signal,
                                      const std::vector<double> &timestamps,
                                      int peak_threshold_percent)
{
  int n = signal.size();
  std::vector<Peak> peaks;

  if (n < 3)
  {
    return peaks;
  }

  double sum = 0.0;
  for (double val : signal)
  {
    sum += val;
  }

  double mean_val = sum / n;
  double threshold = mean_val * (1 + peak_threshold_percent / 100.0);

  std::cout << "Peak detection: mean=" << mean_val
            << ", threshold=" << threshold
            << " (" << peak_threshold_percent << "% above mean)" << std::endl;

  std::vector<Peak> all_peaks;
  for (int i = 1; i < n - 1; ++i)
  {
    if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1])
    {
      all_peaks.push_back({i, signal[i], timestamps[i]});
    }
  }

  for (const auto &peak : all_peaks)
  {
    if (peak.value > threshold)
    {
      peaks.push_back(peak);
    }
  }

  std::cout << "Found " << all_peaks.size() << " local peaks, "
            << peaks.size() << " above threshold" << std::endl;

  return peaks;
}

// If the peaks are very close, consider them as one.
std::vector<int> filter_peaks_by_time(const std::vector<Peak> &peaks,
                                      const std::vector<double> &timestamps,
                                      double min_time_between_peaks)
{
  if (peaks.empty())
  {
    return {};
  }

  std::vector<int> filtered_locs;
  filtered_locs.push_back(peaks[0].index);

  for (size_t i = 1; i < peaks.size(); ++i)
  {
    if ((peaks[i].timestamp - timestamps[filtered_locs.back()]) >= min_time_between_peaks)
    {
      filtered_locs.push_back(peaks[i].index);
    }
  }

  return filtered_locs;
}

// Perform segmented magnetometer calibration.
std::vector<neo_vec3> calibrate_segmented_ransac(
    const std::vector<double> &timestamps,
    const std::vector<neo_vec3> &measurements,
    const magseg_opts &options)
{
  int n = measurements.size();
  std::vector<neo_vec3> calibrated_measurements(n);

  std::cout << "Computing rolling standard deviation..." << std::endl;
  auto rolling_std = compute_rolling_std(measurements, options.window_size, options.alpha);

  std::cout << "Finding peaks..." << std::endl;
  auto peaks = find_relative_peaks(rolling_std, timestamps, options.peak_threshold_percent);
  std::cout << "Found " << peaks.size() << " potential peaks" << std::endl;

  auto valid_locs = filter_peaks_by_time(peaks, timestamps, options.min_time_between_peaks);
  std::cout << "Filtered to " << valid_locs.size() << " peaks after time filtering" << std::endl;

  std::vector<int> boundaries;
  boundaries.push_back(0);
  for (int loc : valid_locs)
  {
    boundaries.push_back(loc);
  }
  boundaries.push_back(n);

  int num_segments = boundaries.size() - 1;
  std::cout << "Divided data into " << num_segments << " segments" << std::endl;

  // Step 5: Calibrate each segment with RANSAC
  for (int seg = 0; seg < num_segments; ++seg)
  {
    int idx_start = boundaries[seg];
    int idx_end = boundaries[seg + 1] - 1;
    int segment_size = idx_end - idx_start + 1;

    std::cout << std::endl;
    std::cout << "Calibrating segment " << (seg + 1) << "/" << num_segments
              << " (indices " << idx_start << " to " << idx_end
              << ", size: " << segment_size << ")" << std::endl;

    if (segment_size < 10)
    {
      std::cout << "Segment too small for RANSAC (min " << 10 << " required), using original data" << std::endl;
      for (int i = idx_start; i <= idx_end; ++i)
      {
        calibrated_measurements[i] = measurements[i];
      }
      continue;
    }

    // Extract segment data
    Eigen::VectorXd x(segment_size), y(segment_size), z(segment_size);
    for (int i = 0; i < segment_size; ++i)
    {
      int global_idx = idx_start + i;
      x(i) = measurements[global_idx].x;
      y(i) = measurements[global_idx].y;
      z(i) = measurements[global_idx].z;
    }

    // Run RANSAC calibration on this segment
    magsac::magsac_stats stats;
    magcal_t calibration = compute_ransac(x, y, z, options.ransac_opts, &stats);
    std::cout << calibration.Ainv << std::endl;
    std::cout << calibration.b << std::endl;

    if (stats.success)
    {
      std::cout << "Segment " << (seg + 1) << " RANSAC success: "
                << stats.best_num_inliers << " inliers ("
                << (100.0 * stats.inlier_ratio) << "%)" << std::endl;

      // Apply calibration to current segment
      for (int i = 0; i < segment_size; ++i)
      {
        int global_idx = idx_start + i;
        Eigen::Vector3d m_raw(measurements[global_idx].x,
                              measurements[global_idx].y,
                              measurements[global_idx].z);
        Eigen::Vector3d m_cal = calibration.Ainv * (m_raw - calibration.b);
        calibrated_measurements[global_idx] = {m_cal(0), m_cal(1), m_cal(2)};
      }
    }
    else
    {
      std::cout << "Segment " << (seg + 1) << " RANSAC failed, using original data" << std::endl;
      for (int i = idx_start; i <= idx_end; ++i)
      {
        calibrated_measurements[i] = measurements[i];
      }
    }
  }

  std::cout << "Segmented calibration complete!" << std::endl;
  return calibrated_measurements;
}

// Overload for Eigen matrices
Eigen::MatrixXd calibrate_segmented_ransac(
    const Eigen::VectorXd &timestamps,
    const Eigen::MatrixXd &measurements,
    const magseg_opts &options)
{
  // Convert Eigen matrix to vector of neo_vec3
  std::vector<double> ts(timestamps.data(), timestamps.data() + timestamps.size());
  std::vector<neo_vec3> meas(measurements.rows());

  for (int i = 0; i < measurements.rows(); ++i)
  {
    meas[i] = {measurements(i, 0), measurements(i, 1), measurements(i, 2)};
  }

  // Call the main function
  auto calibrated_vec = calibrate_segmented_ransac(ts, meas, options);

  // Convert back to Eigen matrix
  Eigen::MatrixXd result(calibrated_vec.size(), 3);
  for (size_t i = 0; i < calibrated_vec.size(); ++i)
  {
    result(i, 0) = calibrated_vec[i].x;
    result(i, 1) = calibrated_vec[i].y;
    result(i, 2) = calibrated_vec[i].z;
  }

  return result;
}

} // namespace magseg
