// Custom data types for NeoAHRS.
// rms (2025-05-24 5:13:48 PM)

#ifndef _NEOAHRS_DATA_TYPES_H_
#define _NEOAHRS_DATA_TYPES_H_

#include <array>
#include <vector>

enum neo_status
{
  neo_status_FILE_READ,
  neo_status_CSV_NO_COLUMN,
  neo_status_FILE_WRITE,
  NEO_SUCCESS
};

enum neo_output_flags
{
  NEOAHRS_OUTPUT_NONE            = 0,
  NEOAHRS_OUTPUT_YAW_COMPASS     = 1 << 0, // ym
  NEOAHRS_OUTPUT_YAW_SENSOR_GYRO = 1 << 1, // ys
  NEOAHRS_OUTPUT_TILT_SENSOR     = 1 << 1, // ps, rs
  NEOAHRS_OUTPUT_ACCEL_LINEAR    = 1 << 2  // alx, aly, alz
};

struct neo_output
{
  std::vector<double> t;
  std::vector<double> ym;
  std::vector<double> ys;
  std::vector<double> ps;
  std::vector<double> rs;
  std::vector<double> alx;
  std::vector<double> aly;
  std::vector<double> alz;
  std::vector<double> afx;
  std::vector<double> afy;
  std::vector<double> afz;
  std::vector<double> gfx;
  std::vector<double> gfy;
  std::vector<double> gfz;

  neo_status status;
};

struct neo_vec3
{
  double x; // X-axis component
  double y; // Y-axis component
  double z; // Z-axis component
};

struct neo_log
{
  // Timestamps
  double t;  // Android        [ms]
  double ta; // Accelerometer  [ms]
  double tg; // Gyroscope      [ms]
  double tm; // Magnetometer   [ms]

  // Raw sensor measurements
  neo_vec3 a; // Accelerometer [m/s^2]
  neo_vec3 g; // Gyroscope     [deg/s]
  neo_vec3 m; // Magnetometer  [uT]
};

struct neo_sensor
{
  double t; // Base time       [s]

  // Sensor measurements
  neo_vec3 a; // Accelerometer [m/s^2]
  neo_vec3 g; // Gyroscope     [deg/s]
  neo_vec3 m; // Magnetometer  [uT]
};

struct neo_euler
{
  double y; // Yaw             [deg]
  double p; // Pitch           [deg]
  double r; // Roll            [deg]
};

struct neo_quat
{
  // q = q0 + q1 i + q2 j + q3 k
  double q0;
  double q1;
  double q2;
  double q3;
};

inline std::vector<double> get_yaw(const std::vector<neo_euler>& yprs)
{
    std::vector<double> result;
    result.reserve(yprs.size());
    for (const auto& ypr : yprs)
        result.push_back(ypr.y);
    return result;
}

inline std::vector<double> get_pitch(const std::vector<neo_euler>& yprs)
{
    std::vector<double> result;
    result.reserve(yprs.size());
    for (const auto& ypr : yprs)
        result.push_back(ypr.p);
    return result;
}

inline std::vector<double> get_roll(const std::vector<neo_euler>& yprs)
{
    std::vector<double> result;
    result.reserve(yprs.size());
    for (const auto& ypr : yprs)
        result.push_back(ypr.r);
    return result;
}

inline std::vector<double> get_x(const std::vector<neo_vec3>& vecs)
{
  std::vector<double> result;
  result.reserve(vecs.size());
  for (const auto& v : vecs)
    result.push_back(v.x);
  return result;
}

inline std::vector<double> get_y(const std::vector<neo_vec3>& vecs)
{
  std::vector<double> result;
  result.reserve(vecs.size());
  for (const auto& v : vecs)
    result.push_back(v.y);
  return result;
}

inline std::vector<double> get_z(const std::vector<neo_vec3>& vecs)
{
  std::vector<double> result;
  result.reserve(vecs.size());
  for (const auto& v : vecs)
    result.push_back(v.z);
  return result;
}

inline std::vector<neo_vec3> get_accel(const std::vector<neo_sensor>& sensor)
{
  std::vector<neo_vec3> result;
  result.reserve(sensor.size());
  for (const auto& s : sensor)
    result.push_back(s.a);
  return result;
}

inline std::vector<neo_vec3> get_gyro(const std::vector<neo_sensor>& sensor)
{
  std::vector<neo_vec3> result;
  result.reserve(sensor.size());
  for (const auto& s : sensor)
    result.push_back(s.g);
  return result;
}

inline std::vector<neo_vec3> get_magn(const std::vector<neo_sensor>& sensor)
{
  std::vector<neo_vec3> result;
  result.reserve(sensor.size());
  for (const auto& s : sensor)
    result.push_back(s.m);
  return result;
}

#endif // data_types.h
