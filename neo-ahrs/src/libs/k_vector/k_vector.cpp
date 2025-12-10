#include "k_vector.h"
#include <cassert>

k_vector::k_vector(const std::vector<double>& sorted_db) : s(sorted_db)
{
  if (s.empty())
  {
    throw std::invalid_argument("Input database cannot be empty");
  }

  // Safety check for maximum size
  if (s.size() > MAX_DATABASE_SIZE)
  {
    throw std::invalid_argument("Input database exceeds maximum allowed size");
  }

  // Verify the input is sorted
  for (size_t i = 1; i < s.size(); ++i)
  {
    if (s[i] < s[i-1])
    {
      throw std::invalid_argument("Input database must be sorted");
    }
  }

  n = s.size();
  generate_kvector();
}

std::vector<double> k_vector::search(double ya, double yb) const
{
  if (std::isnan(ya) || std::isnan(yb))
  {
    throw std::invalid_argument("Search range cannot be NaN");
  }

  if (ya > yb) std::swap(ya, yb); // Ensure ya <= yb

  // Limit the search range to database bounds
  ya = std::max(ya, s.front());
  yb = std::min(yb, s.back());

  if (ya > s.back() || yb < s.front())
  {
    return {}; // No possible matches
  }

  // Calculate indices with safety checks
  double jb_d = (ya - q) / m;
  double jt_d = (yb - q) / m;

  if (!std::isfinite(jb_d) || !std::isfinite(jt_d))
  {
    return {}; // Handle numerical overflow
  }

  size_t jb = get_lower_integer(jb_d);
  size_t jt = get_higher_integer(jt_d);

  // Clamp indices to valid range
  jb = std::min(jb, n - 1);
  jt = std::min(jt, n - 1);

  size_t k_start = (jb == 0) ? 0 : k[jb - 1];
  size_t k_end = (jt == 0) ? 0 : k[jt - 1];

  // Adjust boundaries for more accurate results
  check_boundaries(ya, yb, &k_start, &k_end);

  // Return the range with bounds checking
  if (k_end >= s.size() || k_start > k_end)
  {
    return {};
  }

  return std::vector<double>(s.begin() + k_start, s.begin() + k_end + 1);
}

void k_vector::generate_kvector()
{
  const double y_min = s.front();
  const double y_max = s.back();
  const double xi = K_VECTOR_TOLERANCE;

  // Handle case where all values are equal
  if (safe_double_compare(y_min, y_max))
  {
    m = 0.0;
    q = y_min;
    k.assign(n, n-1);
    return;
  }

  m = (y_max - y_min + 2 * xi) / static_cast<double>(n - 1);
  q = y_min - m - xi;

  k.resize(n);
  k[0] = 0;

  size_t cum_sum = 0;
  for (size_t i = 1; i < n; ++i)
  {
    const double z = m * (i + 1) + q;
    while (cum_sum < n && s[cum_sum] < z)
    {
      ++cum_sum;
    }
    k[i] = std::min(cum_sum, n - 1); // Ensure we don't exceed bounds
  }
}

size_t k_vector::get_higher_integer(double x) const
{
  if (!std::isfinite(x))
  {
    return 0;
  }

  double int_part;
  if (std::modf(x, &int_part) > EPSILON)
  {
    if (int_part >= static_cast<double>(std::numeric_limits<size_t>::max()) - 1)
    {
      return std::numeric_limits<size_t>::max();
    }
    return static_cast<size_t>(int_part + 1);
  }
  return static_cast<size_t>(int_part);
}

size_t k_vector::get_lower_integer(double x) const
{
  if (!std::isfinite(x))
  {
    return 0;
  }
  if (x < 0) return 0;
  return static_cast<size_t>(std::floor(x));
}

void k_vector::check_boundaries(double ya, double yb, size_t* k_start, size_t* k_end) const
{
  if (k_start == nullptr || k_end == nullptr) return;

  // Adjust lower bound
  while (*k_start < n && s[*k_start] < ya)
  {
    ++(*k_start);
  }

  // Adjust upper bound
  while (*k_end > 0 && s[*k_end] > yb)
  {
    --(*k_end);
  }

  // Ensure valid range
  if (*k_start >= n || *k_end >= n || *k_start > *k_end)
  {
    *k_start = 0;
    *k_end = 0;
  }
}

bool k_vector::safe_double_compare(double a, double b) const
{
  return std::fabs(a - b) <= EPSILON * std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}
