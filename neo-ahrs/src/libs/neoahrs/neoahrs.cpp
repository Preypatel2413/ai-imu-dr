#include <cmath>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <iostream>

#include "config.h"
#include "magcal.h"
#include "magsac.h"
#include "magseg.h"
#include "neoahrs.h"
#include "k_vector.h"
#include "rapidcsv.h"
#include "attitude.h"
#include "butterworth.h"
#include "madmadgwick.h"

#define R2D 57.2957795131
#define D2R 1.0 / 57.2957795131

// Corresponding variable names
std::string neo_output_name[] =
    {
        "t",
        "ym",
        "ys",
        "ps",
        "rs",
        "alx",
        "aly",
        "alz",
        "afx",
        "afy",
        "afz",
        "gfx",
        "gfy",
        "gfz",
     };

/**
 * @brief Reads and parses a Neolog CSV file.
 *
 * @param[in] fname Input CSV filename.
 * @param[out] data Vector of structs containing all necessary information from the Neolog file.
 *
 * @return Error code indicating the execution status of the function.
 */
neo_status read_file(const std::string fname, std::vector<neo_log> &data)
{
  try
  {
    rapidcsv::Document csv(fname);

    // Read all columns from CSV
    std::vector<double> t = csv.GetColumn<double>("Time(ms)");
    std::vector<double> tg = csv.GetColumn<double>("Gyro Sensor Time(ms)");
    std::vector<double> ta = csv.GetColumn<double>("Accel Sensor Time(ms)");
    std::vector<double> tm = csv.GetColumn<double>("Magn Sensor Time(ms)");

    std::vector<double> ax = csv.GetColumn<double>("Ax");
    std::vector<double> ay = csv.GetColumn<double>("Ay");
    std::vector<double> az = csv.GetColumn<double>("Az");

    std::vector<double> gx = csv.GetColumn<double>("Gx");
    std::vector<double> gy = csv.GetColumn<double>("Gy");
    std::vector<double> gz = csv.GetColumn<double>("Gz");

    std::vector<double> mx = csv.GetColumn<double>("Mx");
    std::vector<double> my = csv.GetColumn<double>("My");
    std::vector<double> mz = csv.GetColumn<double>("Mz");

    // Check if all columns have the same size
    size_t num_rows = t.size();
    if (ta.size() != num_rows || tg.size() != num_rows || tm.size() != num_rows ||
        ax.size() != num_rows || ay.size() != num_rows || az.size() != num_rows ||
        gx.size() != num_rows || gy.size() != num_rows || gz.size() != num_rows ||
        mx.size() != num_rows || my.size() != num_rows || mz.size() != num_rows)
    {
      return neo_status_CSV_NO_COLUMN;
    }

    data.clear();
    data.reserve(num_rows);

    // Populate data vector
    for (size_t i = 0; i < num_rows; ++i)
    {
      neo_log entry;

      // Set timestamps
      entry.t = t[i];
      entry.ta = ta[i];
      entry.tg = tg[i];
      entry.tm = tm[i];

      // Set accelerometer data
      entry.a.x = ax[i];
      entry.a.y = ay[i];
      entry.a.z = az[i];

      // Set gyroscope data
      entry.g.x = gx[i];
      entry.g.y = gy[i];
      entry.g.z = gz[i];

      // Set magnetometer data
      entry.m.x = mx[i];
      entry.m.y = my[i];
      entry.m.z = mz[i];

      data.push_back(entry);
    }

    return NEO_SUCCESS;
  }
  catch (const std::exception &)
  {
    return neo_status_FILE_READ;
  }
}

/**
 * @brief Performs time synchronization of 3-axis sensors to a common time base
 * @param nl Input sensor data logs
 * @param[out] nd Output synchronized sensor data
 * @param dt Time step for output synchronization
 */
void precondition(const std::vector<neo_log> &nl, std::vector<neo_sensor> &nd, double dt)
{
  if (nl.empty())
  {
    nd.clear();
    return;
  }

  // Phase 1: Extract data from neo_log
  std::vector<double> ta, tg, tm;
  std::vector<std::array<double, 3>> ar, gr, mr;

  // Reserve space upfront
  const size_t initial_size = nl.size();
  ta.reserve(initial_size);
  tg.reserve(initial_size);
  tm.reserve(initial_size);
  ar.reserve(initial_size);
  gr.reserve(initial_size);
  mr.reserve(initial_size);

  // Extract raw data
  for (size_t i = 0; i < nl.size(); ++i)
  {
    const auto &data = nl[i];
    {
      ta.push_back(data.ta);
      tg.push_back(data.tg);
      tm.push_back(data.tm);
      ar.push_back({data.a.x, data.a.y, data.a.z});
      gr.push_back({data.g.x, data.g.y, data.g.z});
      mr.push_back({data.m.x, data.m.y, data.m.z});
    }
  }

  // Phase 2: Find unique and valid measurements
  auto filter_valid = [](const auto &times, const auto &data)
  {
    std::vector<double> valid_times;
    std::vector<std::array<double, 3>> valid_data;

    if (times.empty())
    {
      return std::make_pair(valid_times, valid_data);
    }

    valid_times.reserve(times.size());
    valid_data.reserve(data.size());

    // First element (if valid)
    if (times[0] > 0)
    {
      valid_times.push_back(times[0] * 0.001); // Convert to seconds
      valid_data.push_back(data[0]);
    }

    // Subsequent elements
    for (size_t i = 1; i < times.size(); ++i)
    {
      if (times[i] > times[i - 1] && times[i] > 0)
      {
        valid_times.push_back(times[i] * 0.001);
        valid_data.push_back(data[i]);
      }
    }

    return std::make_pair(valid_times, valid_data);
  };

  auto [ta_valid, a_valid] = filter_valid(ta, ar);
  auto [tg_valid, g_valid] = filter_valid(tg, gr);
  auto [tm_valid, m_valid] = filter_valid(tm, mr);

  // Early exit if any sensor has no valid data
  if (ta_valid.empty() || tg_valid.empty() || tm_valid.empty())
  {
    nd.clear();
    return;
  }

  // Phase 3: Build k_vectors for fast searching
  k_vector kva(ta_valid);
  k_vector kvg(tg_valid);
  k_vector kvm(tm_valid);

  // Phase 4: Determine common time range
  const double time_start = std::max({ta_valid.front(), tg_valid.front(), tm_valid.front()});
  const double time_end = std::min({ta_valid.back(), tg_valid.back(), tm_valid.back()});

  if (time_start > time_end)
  {
    nd.clear();
    return;
  }

  // Phase 5: Generate common time vector
  std::vector<double> t;
  const size_t expected_points = static_cast<size_t>((time_end - time_start) / dt) + 2;
  t.reserve(expected_points);

  for (double current = time_start; current <= time_end; current += dt)
  {
    t.push_back(current);
  }

  // Ensure we include the end point
  if (t.empty() || std::abs(t.back() - time_end) > 1e-9)
  {
    t.push_back(time_end);
  }

  // Phase 6: Interpolation
  nd.resize(t.size());
  const double t0 = t[0];

  for (size_t i = 0; i < t.size(); ++i)
  {
    const double current_t = t[i];

    // Interpolation helper
    auto interpolate = [](const k_vector &kv, const auto &t_vec, const auto &data, double t)
    {
      auto bounds = kv.search(t, t);

      if (bounds.empty())
      {
        return (t <= t_vec.front()) ? data.front() : data.back();
      }

      size_t j = std::lower_bound(t_vec.begin(), t_vec.end(), t) - t_vec.begin();
      j = std::max<size_t>(1, j) - 1;

      const double alpha = (t - t_vec[j]) / (t_vec[j + 1] - t_vec[j]);
      std::array<double, 3> result;

      for (int axis = 0; axis < 3; ++axis)
      {
        result[axis] = data[j][axis] * (1 - alpha) + data[j + 1][axis] * alpha;
      }

      return result;
    };

    // Interpolate all sensors
    auto a_interp = interpolate(kva, ta_valid, a_valid, current_t);
    auto g_interp = interpolate(kvg, tg_valid, g_valid, current_t);
    auto m_interp = interpolate(kvm, tm_valid, m_valid, current_t);

    // Store results
    nd[i] = {
        current_t - t0,
        {a_interp[0], a_interp[1], a_interp[2]},
        {g_interp[0], g_interp[1], g_interp[2]},
        {m_interp[0], m_interp[1], m_interp[2]}};
  }
}

/**
 * @brief Performs tilt correction for 3-axis sensors (accelerometer, gyroscope, and compass).
 *
 * @note Only tilt (not yaw) is corrected; sensor measurements are derotated accordingly.
 *
 * @param s 3-axis sensor measurements in an arbitrary orientation.
 * @param ypr Sensor orientation relative to the NED frame, expressed as 321 Euler angles.
 */
std::vector<neo_sensor> correct_tilt(const std::vector<neo_sensor> &s, const std::vector<neo_euler> &ypr)
{
  std::vector<neo_sensor> output;
  output.reserve(s.size());

  for (size_t i = 0; i < s.size(); ++i)
  {
    const neo_euler &ei = ypr[i];

    // Only use pitch and roll for tilt compensation
    double e[3] = {0.0, ei.p * D2R, ei.r * D2R}; // [yaw, pitch, roll] in radians

    double dcm[3][3], dcm_t[3][3];
    euler_to_dcm(e, EULER_ZYX, dcm); // Returns DCM body->NED
    dcm_trans(dcm, dcm_t);           // Get DCM NED->body (transpose)

    auto apply_dcm = [&](const neo_vec3 &v) -> neo_vec3
    {
      double vin[3] = {v.x, v.y, v.z};
      double vout[3];
      dcm_rotate(dcm_t, vin, vout);
      return {vout[0], vout[1], vout[2]};
    };

    neo_sensor so;
    so.t = s[i].t;
    so.a = apply_dcm(s[i].a);
    so.g = apply_dcm(s[i].g);
    so.m = apply_dcm(s[i].m);

    output.push_back(so);
  }

  return output;
}

/**
 * @brief Calibrates sensor coordinate frames from the Android frame of reference to align them with
 *        the NeoAHRS reference frame.
 *
 * Converting the Android frame to the desired phone frame requires rotation of 90 degrees about the
 * Z-axis.
 *
 *  @note The goal is to define the coordinate frames of the vehicle and the phone so that
 *        discussions and visualizations of orientation become intuitive. The way we define the
 *        frames is based on four key motivations:
 *
 *        1. The vehicle should have zero roll, pitch, and yaw when resting on a flat road. The road
 *           is considered flat relative to gravity.
 *
 *        2. A phone placed on a flat surface (with respect to gravity), screen facing up, should
 *           exhibit zero roll and pitch.
 *
 *        3. The Madgwick filter expects the accelerometer to measure [0, 0, +g] when the sensor
 *           frame is aligned with NED (gravity pointing downward). However, an accelerometer at
 *           rest (e.g., placed on a table) measures the reaction force from the surface rather than
 *           gravity itself, resulting in [0, 0, -g]. To align the measurement with the filter's
 *           expected reference frame, we flip the z-axis of the IMU.
 *
 *        4. The phone's X-axis (vector from rear to front) serves as the reference for North
 *           tracking i.e., heading is 0 degree when the X-axis points toward geographic North. The
 *           same convention applies to the vehicle: North is defined relative to the X-axis.
 *
 * @verbatim
 *           ____________________________________________________________________________
 *          |      (Y)               |        (X)               |              (X)       |
 *          |       ^                |         ^                |               ^        |
 *          |   ____!____            |     ____!____            |           ____!____    |
 *          |  || FRONT ||           |   (|| FRONT ||)          |          || FRONT ||   |
 *          |  ||       ||           |   (||       ||)          |          ||       ||   |
 *          |  ||       ||           |    ||       ||           |          ||       ||   |
 *          |  ||  (.)  ||----> (X)  |    ||  (x)  ||----> (Y)  | (Y) <----||  (.)  ||   |
 *          |  ||   Z   ||           |    ||   Z   ||           |          ||   Z   ||   |
 *          |  ||       ||           |   (||       ||)          |          ||       ||   |
 *          |  ||==(O)==||           |   (||_______||)          |          ||==(O)==||   |
 *          |                        |        REAR              |                        |
 *          |                        |                          |                        |
 *          |  Fig 1. Android        | Fig 2. Vehicle (NeoAHRS) | Fig 3. Phone (NeoAHRS) |
 *          |________________________|__________________________|________________________|
 * @endverbatim
 *
 * @param s Time-stamp corrected sensor measurement.
 *
 * @return Sensor measurements which aligns with NeoAHRS phone coordinate system.
 */
std::vector<neo_sensor> calibrate_frame(const std::vector<neo_sensor> &s)
{
  std::vector<neo_sensor> output;
  output.reserve(s.size());

  // Construct rotation matrix
  double R[3][3];
  double e[3] = {0.0, 0.0, 90 * D2R};
  euler_to_dcm(e, EULER_ZYX, R);

  auto apply_dcm = [&](const neo_vec3 &v) -> neo_vec3
  {
    double vin[3] = {v.x, v.y, v.z};
    double vout[3];
    dcm_rotate(R, vin, vout);
    return {vout[0], vout[1], vout[2]};
  };

  for (size_t i = 0; i < s.size(); ++i)
  {
    // Rotate sensor
    neo_sensor so;
    so.t = s[i].t;
    so.a = apply_dcm(s[i].a);
    so.g = apply_dcm(s[i].g);
    so.m = apply_dcm(s[i].m);

    output.push_back(so);
  }

  // Perform second order Butterworth LPF
  const double fc = NEOAHRS_BUTTERWORTH_FC; // Cut-off frequency
  const double fs = NEOAHRS_BUTTERWORTH_FS; // Sampling frequency
  butterworth filta(fc, fs);
  butterworth filtg(fc, fs);
  butterworth filtm(fc, fs);

  std::vector<neo_vec3> a = filta.filtfilt(get_accel(output));
  std::vector<neo_vec3> g = filtg.filtfilt(get_gyro(output));
  std::vector<neo_vec3> m = filtm.filtfilt(get_magn(output));

  for (size_t i = 0; i < s.size(); ++i)
  {
    output[i].a = {a[i].x, a[i].y, a[i].z};
    output[i].g = {g[i].x, g[i].y, g[i].z};
    output[i].m = {m[i].x, m[i].y, m[i].z};
  }

  return output;
}

/**
 * @brief Computes the magnetometer heading relative to magnetic North (in degrees) using tilt-corrected measurements in the NED frame.
 *
 * The heading is the angle between magnetic North and the x-axis of the tilt-corrected magnetometer.
 * Ensure the sensor is tilt-corrected before calling this function.
 */
std::vector<double> get_mag_heading(const std::vector<neo_vec3> &m_vec)
{
  std::vector<double> heading;
  heading.reserve(m_vec.size());

  for (const auto &m : m_vec)
  {
    double angle_rad = std::atan2(m.y, m.x);
    double angle_deg = angle_rad * R2D;

    // Normalize to [0, 360)
    if (angle_deg < 0)
    {
      angle_deg += 360.0;
    }

    heading.push_back(angle_deg);
  }

  return heading;
}

/**
 * @brief Computes roll and pitch (in degrees) from 3-axis accelerometer measurements.
 *
 * If the accelerometer frame is defined as the NeoAHRS phone frame, the resulting tilt is relative
 * to the NED frame.
 *
 * @param a_vec 3-axis accelerometer measurements.
 * @return 321 Euler angles with roll and pitch populated; yaw is set to zero.
 */
std::vector<neo_euler> get_accel_tilt(const std::vector<neo_vec3> &a_vec)
{
  std::vector<neo_euler> result;
  result.reserve(a_vec.size());

  for (const auto &a : a_vec)
  {
    double ax = a.x;
    double ay = a.y;
    double az = a.z;

    double roll_rad = std::atan2(ay, az);
    double pitch_rad = std::atan2(-ax, std::sqrt(ay * ay + az * az));

    double roll_deg = roll_rad * R2D;
    double pitch_deg = pitch_rad * R2D;

    result.push_back({0.0, pitch_deg, roll_deg});
  }

  return result;
}

/**
 * @brief Computes linear acceleration (without gravity) in NED using accelerometer measurements and
 *        the quaternion orientation of the accelerometer in NED.
 *
 * @param a  3-axis accelerometer measurements.
 * @param q Orientation of the accelerometer in NED.
 * @param g Reference acceleration due to gravity [m/s^2].
 *
 * @return 3-axis linear acceleration of the accelerometer.
 */
std::vector<neo_vec3> get_linaccel(const std::vector<neo_vec3> &a, const std::vector<neo_euler> &ypr, double g)
{
  std::vector<neo_vec3> al_ned;
  size_t n = a.size();

  al_ned.reserve(n);

  for (size_t i = 0; i < n; ++i)
  {
    const neo_euler &ei = ypr[i];
    double e[3] = {0.0, ei.p * D2R, ei.r * D2R};

    double dcm[3][3], dcm_t[3][3];
    euler_to_dcm(e, EULER_ZYX, dcm);
    dcm_trans(dcm, dcm_t);

    // Rotate acceleration into NED frame
    double acc_ned[3];
    double acc_sensor[3] = {a[i].x, a[i].y, a[i].z};
    dcm_rotate(dcm_t, acc_sensor, acc_ned);

    // Subtract gravity
    neo_vec3 lin_accel;
    lin_accel.x = 0.0 - acc_ned[0];
    lin_accel.y = 0.0 - acc_ned[1];
    lin_accel.z = g - acc_ned[2];

    al_ned.push_back(lin_accel);
  }

  return al_ned;
}

neo_output execute(std::vector<neo_sensor> &sensor, int output_flag)
{
  neo_output output;
  double g_accel = NEOAHRS_PHY_G; // Reference acceleration due to gravity

  // Unpack sensor measurements for Madgwick filter
  std::vector<double> t;
  std::vector<neo_vec3> a, g, m;

  for (const auto &s : sensor)
  {
    a.push_back(s.a);
    g.push_back(s.g);
    m.push_back(s.m);
    t.push_back(s.t);
  }

  output.t = t;

  // Perform Madgwick filter
  double beta = NEOAHRS_MADMADGWICK_BETA;
  std::vector<madgwick_t> mad = madmadgwick(a, g, t, beta);

  // Configure segmentation
  magseg::magseg_opts options;
  options.alpha = NEOAHRS_MAGSEG_ALPHA;
  options.window_size = NEOAHRS_MAGSEG_WINDOW_SIZE;
  options.min_time_between_peaks = NEOAHRS_MIN_TIME_BETWEEN_PEAKS_SEC;
  options.peak_threshold_percent = NEOAHRS_THRESHOLD_PERCENT;

  // Configure RANSAC for each segment
  options.ransac_opts.min_num_iterations = NEOAHRS_RANSAC_MIN_NUM_ITERATIONS;
  options.ransac_opts.max_num_iterations = NEOAHRS_RANSAC_MAX_NUM_ITERATIONS;
  options.ransac_opts.inlier_threshold = NEOAHRS_RANSAC_INLIER_THRESHOLD;
  options.ransac_opts.min_inlier_percent = NEOAHRS_RANSAC_INLIER_PERCENT;

  auto mc = magseg::calibrate_segmented_ransac(t, m, options);

  for (size_t i = 0; i < sensor.size(); ++i)
  {
    sensor[i].m.x = mc[i].x;
    sensor[i].m.y = mc[i].y;
    sensor[i].m.z = mc[i].z;
  }

  // Tilt correction
  std::vector<neo_euler> ypr_s = get_euler(mad);
  output.ys = get_yaw(ypr_s);
  output.ps = get_pitch(ypr_s);
  output.rs = get_roll(ypr_s);

  std::vector<neo_sensor> flat_sensor = correct_tilt(sensor, ypr_s);

  std::vector<neo_vec3> af, gf, mf;

  for (size_t i = 0; i < flat_sensor.size(); ++i)
  {
    af.push_back(flat_sensor[i].a);
    gf.push_back(flat_sensor[i].g);
    mf.push_back(flat_sensor[i].m);
  }

   output.afx = get_x(af);
   output.afy = get_y(af);
   output.afz = get_z(af);
   output.gfx = get_x(gf);
   output.gfy = get_y(gf);
   output.gfz = get_z(gf);

   // Compute Magnetometer bearing
  output.ym = get_mag_heading(mf);

  // if (output_flag & NEOAHRS_OUTPUT_ACCEL_LINEAR)
  // {
    // Compute linear acceleration
    std::vector<neo_quat> qs = get_quat(mad);
    std::vector<neo_vec3> al = get_linaccel(a, ypr_s, g_accel);

    output.alx = get_x(al);
    output.aly = get_y(al);
    output.alz = get_z(al);
  // }

  return output;
}

void neoahrs::export_csv(const std::string &filename, const neo_output &output)
{
  std::ofstream file(filename);
  if (!file.is_open())
  {
    std::cerr << "Failed to open output file.\n";
    return;
  }

  // Create a mapping between struct members and their names
  const auto &output_ref = output;

  // List of member pointers paired with their names
  const std::pair<const std::vector<double> *, std::string> field_map[] =
  {
    {&output_ref.t, neo_output_name[0]},
    {&output_ref.ym, neo_output_name[1]},
    {&output_ref.ys, neo_output_name[2]},
    {&output_ref.ps, neo_output_name[3]},
    {&output_ref.rs, neo_output_name[4]},
    {&output_ref.alx, neo_output_name[5]},
    {&output_ref.aly, neo_output_name[6]},
    {&output_ref.alz, neo_output_name[7]},
    {&output_ref.afx, neo_output_name[8]},
    {&output_ref.afy, neo_output_name[9]},
    {&output_ref.afz, neo_output_name[10]},
    {&output_ref.gfx, neo_output_name[11]},
    {&output_ref.gfy, neo_output_name[12]},
    {&output_ref.gfz, neo_output_name[13]}
  };

  // Collect non-empty columns
  std::vector<std::pair<const std::vector<double> *, std::string>> columns;
  for (const auto &[vec_ptr, name] : field_map)
  {
    if (!vec_ptr->empty())
    {
      columns.emplace_back(vec_ptr, name);
    }
  }

  if (columns.empty())
  {
    std::cerr << "No data to write.\n";
    return;
  }

  // Write header
  for (size_t i = 0; i < columns.size(); ++i)
  {
    if (i > 0)
      file << ",";
    file << columns[i].second;
  }
  file << "\n";

  // Get row count and validate consistent sizes
  const size_t N = columns[0].first->size();
  for (const auto &[vec_ptr, _] : columns)
  {
    if (vec_ptr->size() != N)
    {
      std::cerr << "Error: Columns have different sizes.\n";
      return;
    }
  }

  std::cout << std::endl;
  std::cout << "variable: size" << std::endl;
  for (const auto &[vec_ptr, name] : field_map)
  {
    std::cout << name << ": " << vec_ptr->size() << "\n";
  }
  std::cout << std::endl;

  // Write data
  for (size_t row = 0; row < N; ++row)
  {
    for (size_t col = 0; col < columns.size(); ++col)
    {
      if (col > 0)
      {
        file << ",";
      }
      file << (*columns[col].first)[row];
    }
    file << "\n";
  }

  file.close();
}

// Perform NeoAHRS computation for input Neolog file and return the output specified in output_flag.
neo_output neoahrs::neoahrs(std::string &filename, int output_flag)
{
  // Parse the raw Neolog measurements
  std::vector<neo_log> neolog_data;

  if (read_file(filename, neolog_data) != NEO_SUCCESS)
  {
    std::cerr << "Error reading input CSV file." << std::endl;
  }

  return neoahrs(neolog_data, output_flag);
}

// Perform NeoAHRS computation for input data and return the output specified in output_flag.
neo_output neoahrs::neoahrs(std::vector<neo_log> &data, int output_flag)
{
  // Preprocess sensor measurements
  std::vector<neo_sensor> sensor_raw;
  double dt = NEOAHRS_SENSOR_SAMPLE_TIME_SEC; // Sensor sample time [s]

  precondition(data, sensor_raw, dt);
  std::vector<neo_sensor> sensor = calibrate_frame(sensor_raw);

  neo_output output = execute(sensor, output_flag);

  return output;
}
