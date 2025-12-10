/**
 * @brief Implementation of k-vector range searching algorithm
 * @cite Mortari, Neta - k-Vector range searching techniques (2014)
 * @author rms
 * @date 2023-04-10
**/

#ifndef _K_VECTOR_H_
#define _K_VECTOR_H_

#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>
#include <algorithm>

class k_vector
{
private:
  double m, q;           // z(x) = m * x + q
  std::vector<size_t> k; // k-vector
  std::vector<double> s; // Sorted database
  size_t n;              // Length of s

  void generate_kvector();
  size_t get_higher_integer(double x) const;
  size_t get_lower_integer(double x) const;
  void check_boundaries(double ya, double yb, size_t* k_start, size_t* k_end) const;
  bool safe_double_compare(double a, double b) const;

public:
  explicit k_vector(const std::vector<double>& sorted_db);
  std::vector<double> search(double ya, double yb) const;

  // Safety parameters
  static constexpr double K_VECTOR_TOLERANCE = 0.02;
  static constexpr double EPSILON = std::numeric_limits<double>::epsilon() * 100;
  static constexpr size_t MAX_DATABASE_SIZE = 1000000; // Safety limit
};

#endif // k_vector.h
