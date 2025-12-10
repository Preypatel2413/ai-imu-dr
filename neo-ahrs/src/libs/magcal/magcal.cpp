#include "magcal.h"

#include <vector>
#include <tuple>
#include <cmath>
#include <stdexcept>

using namespace Eigen;

VectorXd magcal_fit_ellipsoid10(const VectorXd &x, const VectorXd &y, const VectorXd &z);
VectorXd magcal_fit_ellipsoid7(const VectorXd &x, const VectorXd &y, const VectorXd &z);
VectorXd magcal_fit_ellipsoid4(const VectorXd &x, const VectorXd &y, const VectorXd &z);

std::tuple<Matrix3d, Matrix3d, Vector3d> get_calibration_params(const VectorXd &v);
double get_field_intensity(const VectorXd &v, const Matrix3d &A, const Vector3d &b);
std::pair<double, bool> get_fit_quality(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                                        const Matrix3d &A, const Matrix3d &Ainv, const Vector3d &b, double r);
std::tuple<Matrix3d, Matrix3d, Vector3d, double, double, bool> get_all_params(const VectorXd &x, const VectorXd &y, const VectorXd &z, int idx);

VectorXd magcal_fit_ellipsoid10(const VectorXd &x, const VectorXd &y, const VectorXd &z)
{
  // Check input sizes
  assert(x.size() == y.size() && y.size() == z.size());
  int n = x.size();

  // Design matrix D (10 x n)
  MatrixXd D(10, n);
  D.row(0) = x.array().square().transpose();
  D.row(1) = y.array().square().transpose();
  D.row(2) = z.array().square().transpose();
  D.row(3) = (2.0 * y.array() * z.array()).transpose();
  D.row(4) = (2.0 * x.array() * z.array()).transpose();
  D.row(5) = (2.0 * x.array() * y.array()).transpose();
  D.row(6) = (2.0 * x).transpose();
  D.row(7) = (2.0 * y).transpose();
  D.row(8) = (2.0 * z).transpose();
  D.row(9) = VectorXd::Ones(n).transpose();

  MatrixXd invC(6, 6);
  invC << 0, 0.5, 0.5, 0, 0, 0,
      0.5, 0, 0.5, 0, 0, 0,
      0.5, 0.5, 0, 0, 0, 0,
      0, 0, 0, -0.25, 0, 0,
      0, 0, 0, 0, -0.25, 0,
      0, 0, 0, 0, 0, -0.25;

  MatrixXd S = D * D.transpose();
  MatrixXd S11 = S.block(0, 0, 6, 6); // 6x6
  MatrixXd S12 = S.block(0, 6, 6, 4); // 6x4
  MatrixXd S22 = S.block(6, 6, 4, 4); // 4x4

  JacobiSVD<MatrixXd> svd(S22, ComputeThinU | ComputeThinV);
  double tolerance = std::numeric_limits<double>::epsilon() * std::max(S22.cols(), S22.rows()) * svd.singularValues()(0);
  MatrixXd Sx = svd.matrixV() *
                (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() *
                svd.matrixU().adjoint() * S12.transpose();

  // Eqn(14) and Eqn(15)
  MatrixXd M = invC * (S11 - S12 * Sx);
  EigenSolver<MatrixXd> es(M);

  // Index of the 'only' positive eigenvalue
  int max_cidx;
  es.eigenvalues().real().maxCoeff(&max_cidx);

  // Eigenvector corresponding to max_cidx
  VectorXd u1 = es.eigenvectors().real().col(max_cidx);
  VectorXd u2 = -Sx * u1;
  VectorXd u(10);
  u << u1, u2;

  Matrix3d A;
  A << u(0), u(5), u(4),
      u(5), u(1), u(3),
      u(4), u(3), u(2);

  // Ensure positive definiteness
  if (A.determinant() < 0)
  {
    u = -u;
  }

  return u;
}

VectorXd magcal_fit_ellipsoid7(const VectorXd &x, const VectorXd &y, const VectorXd &z)
{
  // Check input sizes
  assert(x.size() == y.size() && y.size() == z.size());
  int n = x.size();

  // Design matrix D (n x 7)
  MatrixXd D(n, 7);
  D.col(0) = x.array().square(); // x^2
  D.col(1) = y.array().square(); // y^2
  D.col(2) = z.array().square(); // z^2
  D.col(3) = 2 * x;              // 2x
  D.col(4) = 2 * y;              // 2y
  D.col(5) = 2 * z;              // 2z
  D.col(6) = VectorXd::Ones(n);  // 1

  // Covariance matrix S (7x7)
  MatrixXd S = D.transpose() * D;

  // Eigen decomposition to find smallest eigenvalue/vector
  SelfAdjointEigenSolver<MatrixXd> eigensolver(S);

  if (eigensolver.info() != Success)
  {
    throw std::runtime_error("Eigen decomposition failed");
  }

  // Find index of smallest eigenvalue
  int min_idx;
  eigensolver.eigenvalues().minCoeff(&min_idx);
  VectorXd v = eigensolver.eigenvectors().col(min_idx);

  // Extract parameters to 10-element vector (with zeros for cross-terms)
  VectorXd u = VectorXd::Zero(10);
  u(0) = v(0); // a (x^2)
  u(1) = v(1); // b (y^2)
  u(2) = v(2); // c (z^2)
  // u(3)-u(5) remain 0 (f, g, h - cross terms)
  u(6) = v(3); // p (x term coefficient)
  u(7) = v(4); // q (y term coefficient)
  u(8) = v(5); // r (z term coefficient)
  u(9) = v(6); // d (constant term)

  // Check positive definiteness of the ellipsoid matrix
  Matrix3d A;
  A << u(0), u(5), u(4), // Note: u(5)=h=0, u(4)=g=0
      u(5), u(1), u(3),  // u(3)=f=0
      u(4), u(3), u(2);

  if (A.determinant() < 0)
  {
    u = -u;
  }

  return u;
}

VectorXd magcal_fit_ellipsoid4(const VectorXd &x, const VectorXd &y, const VectorXd &z)
{
  // Check input sizes
  assert(x.size() == y.size() && y.size() == z.size());
  const int n = x.size();

  // Construct the X matrix (n x 4)
  MatrixXd X(n, 4);
  X.col(0) = 2.0 * x;           // 2x
  X.col(1) = 2.0 * y;           // 2y
  X.col(2) = 2.0 * z;           // 2z
  X.col(3) = VectorXd::Ones(n); // 1

  // Construct Y vector (x² + y² + z²)
  VectorXd Y = x.array().square() + y.array().square() + z.array().square();

  // Solve the linear system XᵀX v = -XᵀY
  MatrixXd XTX = X.transpose() * X;
  VectorXd XTY = X.transpose() * Y;
  VectorXd v = -XTX.ldlt().solve(XTY); // Using LDLT decomposition

  // Extract parameters to 10-element vector (with most terms zero)
  VectorXd u = VectorXd::Zero(10);
  u(0) = 1.0; // a = 1 (x² term)
  u(1) = 1.0; // b = 1 (y² term)
  u(2) = 1.0; // c = 1 (z² term)
  // f,g,h remain 0 (cross terms)
  u(6) = v(0); // p (x term coefficient)
  u(7) = v(1); // q (y term coefficient)
  u(8) = v(2); // r (z term coefficient)
  u(9) = v(3); // d (constant term)

  return u;
}

// Magnitude of magnetic field intensity (normalized radius)
double get_field_intensity(const VectorXd &v, const Matrix3d &A, const Vector3d &b)
{
  double terms = A(0, 0) * b(0) * b(0) + 2 * A(0, 1) * b(0) * b(1) + A(1, 1) * b(1) * b(1) + 2 * A(0, 2) * b(0) * b(2) + A(2, 2) * b(2) * b(2) + 2 * A(1, 2) * b(1) * b(2) - v(9);
  return sqrt(abs(terms)) / sqrt(cbrt(A.determinant()));
}

std::tuple<Matrix3d, Matrix3d, Vector3d> get_calibration_params(const VectorXd &v)
{
  Matrix3d A;
  A << v(0), v(5), v(4),
      v(5), v(1), v(3),
      v(4), v(3), v(2);

  Vector3d k(v(6), v(7), v(8));

  double detA = A.determinant();
  double detA_third = std::cbrt(detA);

  SelfAdjointEigenSolver<Matrix3d> es(A / detA_third);
  Matrix3d sqrtA = es.operatorSqrt();
  Matrix3d Ainv = sqrtA;

  Vector3d b = -A.lu().solve(k);

  return std::make_tuple(A, Ainv, b);
}

// Quality of ellipsoid fit
std::pair<double, bool> get_fit_quality(const VectorXd &x, const VectorXd &y, const VectorXd &z,
                                        const Matrix3d &A, const Matrix3d &Ainv, const Vector3d &b, double r)
{
  // Scaled RMSE: Deviation from a perfect sphere
  MatrixXd m(x.size(), 3);
  m.col(0) = x;
  m.col(1) = y;
  m.col(2) = z;

  MatrixXd m_hat = (Ainv * (m.transpose().colwise() - b)).transpose();
  VectorXd r_sq = m_hat.array().square().rowwise().sum();
  VectorXd err = r_sq.array() - r * r;
  double rmse = sqrt(err.array().square().mean());
  double srmse = rmse / (2 * r * r);

  // Check if positive definite
  LLT<Matrix3d> llt(A);
  bool is_pd = (llt.info() == Success);

  return std::make_pair(srmse, is_pd);
}

std::tuple<Matrix3d, Matrix3d, Vector3d, double, double, bool> get_all_params(const VectorXd &x, const VectorXd &y, const VectorXd &z, int idx)
{
  VectorXd v;

  if (idx == 4)
  {
    v = magcal_fit_ellipsoid4(x, y, z);
  }
  else if (idx == 7)
  {
    v = magcal_fit_ellipsoid7(x, y, z);
  }
  else
  {
    v = magcal_fit_ellipsoid10(x, y, z);
  }

  Matrix3d A, Ainv;
  Vector3d b;
  std::tie(A, Ainv, b) = get_calibration_params(v);

  double r = get_field_intensity(v, A, b);

  double err;
  bool is_pd;
  std::tie(err, is_pd) = get_fit_quality(x, y, z, A, Ainv, b, r);

  // Check for complex components by examining norm of imaginary parts
  if (A.norm() != A.real().norm() ||
      Ainv.norm() != Ainv.real().norm() ||
      b.norm() != b.real().norm())
  {
    is_pd = false;
  }

  return std::make_tuple(A.real(), Ainv.real(), b.real(), r, err, is_pd);
}

magcal_t magcal::compute(const VectorXd &x, const VectorXd &y, const VectorXd &z)
{
  Matrix3d A, Ainv;
  Vector3d b;
  double r, err;
  bool is_pd;

  std::tie(A, Ainv, b, r, err, is_pd) = get_all_params(x, y, z, 4);

  Matrix3d A7, Ainv7;
  Vector3d b7;
  double r7, err7;
  bool is_pd7;
  std::tie(A7, Ainv7, b7, r7, err7, is_pd7) = get_all_params(x, y, z, 7);

  if (is_pd7 && (err7 < err))
  {
    Ainv = Ainv7;
    b = b7;
    r = r7;
    err = err7;
    is_pd = is_pd7;
  }

  Matrix3d A10, Ainv10;
  Vector3d b10;
  double r10, err10;
  bool is_pd10;
  std::tie(A10, Ainv10, b10, r10, err10, is_pd10) = get_all_params(x, y, z, 10);

  if (is_pd10 && (err10 < err))
  {
    Ainv = Ainv10;
    b = b10;
    r = r10;
    err = err10;
    is_pd = is_pd10;
  }

  return {Ainv, b, r, err, is_pd};
}

void magcal::calibrate(const magcal_t &magcal,
                       const VectorXd &x, const VectorXd &y, const VectorXd &z,
                       VectorXd &x_cal, VectorXd &y_cal, VectorXd &z_cal)
{
  const int n = x.size();

  // Create 3xN matrix of raw measurements
  MatrixXd raw(3, n);
  raw.row(0) = x.transpose();
  raw.row(1) = y.transpose();
  raw.row(2) = z.transpose();

  // Subtract bias (broadcast operation)
  MatrixXd centered = raw.colwise() - magcal.b;

  // Apply soft-iron correction
  MatrixXd calibrated = magcal.Ainv * centered;

  // Unpack results
  x_cal = calibrated.row(0).transpose();
  y_cal = calibrated.row(1).transpose();
  z_cal = calibrated.row(2).transpose();
}

#include <iostream>

magcal_t magcal::compute(const std::vector<neo_vec3> &m)
{
  const size_t n = m.size();
  VectorXd x(n), y(n), z(n);

  for (size_t i = 0; i < n; i++)
  {
    x[i] = m[i].x;
    y[i] = m[i].y;
    z[i] = m[i].z;
  }

  return magcal::compute(x, y, z);
}

std::vector<neo_vec3> magcal::magcal(const std::vector<neo_vec3> &m)
{
  const size_t n = m.size();
  VectorXd x(n), y(n), z(n);

  for (size_t i = 0; i < n; i++)
  {
    x[i] = m[i].x;
    y[i] = m[i].y;
    z[i] = m[i].z;
  }

  magcal_t cal_params = magcal::compute(x, y, z);
  magcal::calibrate(cal_params, x, y, z, x, y, z);

  // std::cout << cal_params.Ainv << std::endl;
  // std::cout << cal_params.b << std::endl;

  std::vector<neo_vec3> m_hat(n);

  for (size_t i = 0; i < n; i++)
  {
    m_hat[i] = {x[i], y[i], z[i]};
  }

  return m_hat;
}

std::vector<neo_vec3> magcal::calibrate(const magcal_t params, const std::vector<neo_vec3> &m)
{
    const size_t n = m.size();

  VectorXd x(n), y(n), z(n);

  for (size_t i = 0; i < n; i++)
  {
    x[i] = m[i].x;
    y[i] = m[i].y;
    z[i] = m[i].z;
  }

  magcal::calibrate(params, x, y, z, x, y, z);

  // std::cout << cal_params.Ainv << std::endl;
  // std::cout << cal_params.b << std::endl;

  std::vector<neo_vec3> m_hat(n);

  for (size_t i = 0; i < n; i++)
  {
    m_hat[i] = {x[i], y[i], z[i]};
  }

  return m_hat;
}
