#include "butterworth.h"
#include <iostream>

butterworth::butterworth(double cutoff_freq, double sample_rate)
{
  reset_state();
  compute_coefficients(cutoff_freq, sample_rate);
}

void butterworth::reset_state()
{
  z[0] = z[1] = 0.0;
}

void butterworth::compute_coefficients(double fc, double fs)
{
  // Normalized frequency
  double omega = 2.0 * M_PI * fc / fs;
  double tan_omega = tan(omega / 2);
  double sqrt2 = 1.4142135623730951;

  // Bilinear transform with pre-warping
  double a0 = 1.0 + sqrt2 * tan_omega + tan_omega * tan_omega;

  // Normalized coefficients (matches MATLAB/Octave)
  b[0] = tan_omega * tan_omega / a0;
  b[1] = 2.0 * b[0];
  b[2] = b[0];

  a[0] = 1.0;
  a[1] = 2.0 * (tan_omega * tan_omega - 1.0) / a0;
  a[2] = (1.0 - sqrt2 * tan_omega + tan_omega * tan_omega) / a0;
}

double butterworth::filter(double input)
{
  // Transposed Direct Form II implementation
  double output = b[0] * input + z[0];
  z[0] = b[1] * input - a[1] * output + z[1];
  z[1] = b[2] * input - a[2] * output;

  return output;
}

std::vector<double> butterworth::filter(const std::vector<double> &input)
{
  std::vector<double> output;
  output.reserve(input.size());

  for (double x : input)
  {
    output.push_back(filter(x));
  }

  return output;
}

std::vector<double> butterworth::filtfilt(const std::vector<double> &input)
{
  if (input.size() <= 6)
  {
    throw std::runtime_error("Input signal too short for filtfilt");
  }

  // Reflection length
  const size_t lrefl = 3 * n;

  // Compute initial conditions using DC gain
  double sum_b = b[0] + b[1] + b[2];
  double sum_a = a[0] + a[1] + a[2];
  double kdc = (sum_a != 0) ? sum_b / sum_a : 0;

  double si[n] = {(b[1] - a[1] * kdc) + (b[2] - a[2] * kdc), (b[2] - a[2] * kdc)};

  // Create reflected input
  std::vector<double> padded;
  padded.reserve(input.size() + 2 * lrefl);

  // Beginning reflection
  for (size_t i = 0; i < lrefl; i++)
  {
    padded.push_back(2 * input[0] - input[lrefl - i]);
  }

  // Original signal
  padded.insert(padded.end(), input.begin(), input.end());

  // End reflection
  for (size_t i = 0; i < lrefl; i++)
  {
    padded.push_back(2 * input.back() - input[input.size() - 2 - i]);
  }

  // Forward filter
  std::vector<double> forward = filter_impl(padded, si, padded[0]);

  // Reverse filter
  std::reverse(forward.begin(), forward.end());
  std::vector<double> backward = filter_impl(forward, si, forward[0]);
  std::reverse(backward.begin(), backward.end());

  if (backward.size() < 2 * lrefl + input.size())
    throw std::runtime_error("Unexpected filtfilt output size");

  // Extract central portion
  return std::vector<double>(backward.begin() + lrefl, backward.begin() + lrefl + input.size());
}

std::vector<double> butterworth::filter_impl(const std::vector<double> &input,
                                             const double si[], double scale)
{
  // Set initial state
  z[0] = si[0] * scale;
  z[1] = si[1] * scale;

  return filter(input);
}

std::vector<neo_vec3> butterworth::filtfilt(const std::vector<neo_vec3> &input)
{
  size_t N = input.size();
  std::vector<double> x_in(N), y_in(N), z_in(N);

  // Separate axes
  for (size_t i = 0; i < N; ++i)
  {
    x_in[i] = input[i].x;
    y_in[i] = input[i].y;
    z_in[i] = input[i].z;
  }

  // Filter each axis
  std::vector<double> x_out = filtfilt(x_in);
  std::vector<double> y_out = filtfilt(y_in);
  std::vector<double> z_out = filtfilt(z_in);

  // Recombine
  std::vector<neo_vec3> result(N);
  for (size_t i = 0; i < N; ++i)
  {
    result[i].x = x_out[i];
    result[i].y = y_out[i];
    result[i].z = z_out[i];
  }

  return result;
}
