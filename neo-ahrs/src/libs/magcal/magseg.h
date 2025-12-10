
/**
 * magcal() with RANSAC and drive SEGmentation
 *
 * We segment a drive into sections based on magnetic anomaly. The anomaly is determined by the
 * rolling windowed standard deviation followed by with exponent biasing alpha.
 **/
#ifndef _MAGSEG_H_
#define _MAGSEG_H_

#include <cmath>
#include <vector>
#include <algorithm>

#include "Eigen"
#include "magsac.h"
#include "data_types.h"

namespace magseg
{

struct magseg_opts
{
  magseg_opts()
      : window_size(100),
        alpha(2.0),
        peak_threshold_percent(200), // percentage above mean (like Octave)
        min_time_between_peaks(60.0) // seconds
  {
  }

  magsac::magsac_opts ransac_opts;
  int window_size;               // Window size to compute rolling standard deviation
  double alpha;                  // Exponent weighting of signal for anomaly detection
  int peak_threshold_percent;    // Percentage above mean for peak detection
  double min_time_between_peaks; // If the peaks are too close, consider them as one
};

// Peak detection structure
struct Peak
{
  int index;
  double value;
  double timestamp;
};

// Main segmentation and calibration function
std::vector<neo_vec3> calibrate_segmented_ransac(
    const std::vector<double> &timestamps,
    const std::vector<neo_vec3> &measurements,
    const magseg_opts &options = magseg_opts());

Eigen::MatrixXd calibrate_segmented_ransac(
    const Eigen::VectorXd &timestamps,
    const Eigen::MatrixXd &measurements, // N x 3
    const magseg_opts &options = magseg_opts());

std::vector<Peak> find_relative_peaks(const std::vector<double> &signal,
                                      const std::vector<double> &timestamps,
                                      int peak_threshold_percent);

std::vector<double> compute_rolling_std(const std::vector<neo_vec3> &measurements,
                                        int window_size, double alpha);

std::vector<int> filter_peaks_by_time(const std::vector<Peak> &peaks,
                                      const std::vector<double> &timestamps,
                                      double min_time_between_peaks);

} // namespace magseg

#endif // magseg.h
