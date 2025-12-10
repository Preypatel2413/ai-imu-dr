/**
 * Implementation of second order Butterworth LPF
 *
 * References:
 *   https://github.com/mofr/butterworth
 *   https://www.meme.net.au/butterworth.html
 *
 *   rms (2023-01-30)
 */
#ifndef _BUTTERWORTH_H_
#define _BUTTERWORTH_H_

#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>

#include "data_types.h"

class butterworth
{
public:
  butterworth(double cutoff_freq, double sample_rate);

  double filter(double input);
  std::vector<double> filter(const std::vector<double>& input);
  std::vector<double> filtfilt(const std::vector<double>& input);
  std::vector<neo_vec3> filtfilt(const std::vector<neo_vec3>& input);

private:
  void compute_coefficients(double fc, double fs);
  void reset_state();
  std::vector<double> filter_impl(const std::vector<double>& input,
                                  const double si[], double scale);

  // Filter coefficients
  double b[3]; // numerator coefficients
  double a[3]; // denominator coefficients

  // Internal state (Transposed Direct Form II)
  double z[2];

  // Filter order
  static constexpr int n = 2;
};

#endif // _BUTTERWORTH_H_
