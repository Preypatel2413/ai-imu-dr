import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, detrend

"""
stage3.py
=========

Stage 3: orientation calibration (roll, pitch, yaw) from synchronized IMU+GPS.

This stage:
  - Loads the time-aligned CSV from stage2 (expects columns: time, Ax/Ay/Az, Gx/Gy/Gz).
  - Estimates device roll & pitch by low-pass filtering accelerometer to extract
    the gravity vector and fitting it to [0,0,+g].
  - Rotates accelerometer and gyroscope measurements into a leveled horizontal frame.
  - Searches for the yaw offset that minimizes lateral motion (non-holonomic cost)
    by rotating the horizontal accelerations and minimizing lateral velocity variance.
  - Optionally plots velocity/acceleration traces and the cost curve.
  - Returns (roll_deg, pitch_deg, best_yaw_deg) in degrees.
"""


def stage3(input_path, segment_start = 0, segment_end = None, prnt = False):
    df = pd.read_csv(input_path)
    
    if(segment_start!=0 or segment_end!=None):
        segment_end = len(df) if segment_end==None else segment_end
        df = df[segment_start:segment_end]

    time = df['time'].values
    ax = df['Ax'].values
    ay = df['Ay'].values
    az = df['Az'].values
    gx = df['Gx'].values
    gy = df['Gy'].values
    gz = df['Gz'].values

    fs = 1 / np.mean(np.diff(time))
    print(f"Sampling frequency: {fs:.1f} Hz")

    # ============================= STEP 1: Estimate Roll & Pitch from Gravity =============================
    # Low-pass filter to get gravity vector in device frame (during all periods, including motion)
    cutoff = 0.2  # Hz
    nyq = 0.5 * fs
    b, a = butter(4, cutoff / nyq, btype='low')

    gx_f = filtfilt(b, a, ax)
    gy_f = filtfilt(b, a, ay)
    gz_f = filtfilt(b, a, az)

    # Average gravity vector over the whole sequence (or only static parts if you know them)
    g_vec = np.array([np.mean(gx_f), np.mean(gy_f), np.mean(gz_f)])
    g_norm = np.linalg.norm(g_vec)
    if(prnt):
      print(f"Average gravity vector (device frame): [{g_vec[0]:.3f}, {g_vec[1]:.3f}, {g_vec[2]:.3f}] m/s²")
      print(f"Gravity magnitude: {g_norm:.3f} m/s²")

    # Compute roll and pitch that align this gravity vector to [0, 0, +9.81]
    roll = np.arctan2(g_vec[1], np.sqrt(g_vec[0]**2 + g_vec[2]**2))
    pitch = np.arctan2(-g_vec[0], np.sqrt(g_vec[1]**2 + g_vec[2]**2))

    roll_deg = np.degrees(roll)
    pitch_deg = np.degrees(pitch)

    if(prnt):
      print(f"\n=== ESTIMATED DEVICE MOUNTING ===")
      print(f"Roll  (rotation around X):  {roll_deg:+6.3f}°")
      print(f"Pitch (rotation around Y):  {pitch_deg:+6.3f}°")

    # Build rotation matrix: R_device_to_horizontal = R_pitch @ R_roll
    cr, sr = np.cos(roll), np.sin(roll)
    cp, sp = np.cos(pitch), np.sin(pitch)

    R_roll = np.array([[1, 0, 0],
                    [0, cr, -sr],
                    [0, sr, cr]])

    R_pitch = np.array([[cp, 0, sp],
                        [0, 1, 0],
                        [-sp, 0, cp]])

    R_dev_to_horiz = R_pitch @ R_roll  # Apply roll first, then pitch

    # ============================= STEP 2: Rotate Accel & Gyro into Horizontal Car Frame =============================
    # Rotate linear acceleration (first remove gravity in device frame, then rotate)
    lin_acc_dev = np.vstack((ax, ay, az)).T - np.vstack((gx_f, gy_f, gz_f)).T
    lin_acc_horiz = (R_dev_to_horiz @ lin_acc_dev.T).T   # Shape: (N, 3)

    # Rotate gyroscope (angular rates transform with the same rotation)
    gyro_dev = np.vstack((gx, gy, gz)).T
    gyro_horiz = (R_dev_to_horiz @ gyro_dev.T).T

    a_hx, a_hy, a_hz = lin_acc_horiz[:,0], lin_acc_horiz[:,1], lin_acc_horiz[:,2]
    g_hx, g_hy, g_hz = gyro_horiz[:,0], gyro_horiz[:,1], gyro_horiz[:,2]

    # ============================= STEP 3: Yaw Calibration in Horizontal Frame =============================
    dt = np.mean(np.diff(time))
    def yaw_cost_in_horizontal_frame(yaw_deg):
        theta = np.deg2rad(yaw_deg)
        c, s = np.cos(theta), np.sin(theta)

        # Rotate horizontal accelerations into candidate car frame
        a_forward =  c * a_hx + s * a_hy
        a_lateral = -s * a_hx + c * a_hy

        v_forward = np.cumsum(a_forward) * dt
        v_lateral = np.cumsum(a_lateral) * dt

        cost = np.mean(v_lateral**2)   # Main nonholonomic term
        return cost

    # Coarse search
    # print("\nRunning yaw grid search in leveled horizontal frame...")
    yaws = np.arange(-180, 180, 0.8)
    costs = [yaw_cost_in_horizontal_frame(y) for y in yaws]

    best_coarse_idx = np.argmin(costs)
    best_yaw_coarse = yaws[best_coarse_idx]

    # Fine refinement
    fine_yaws = np.arange(best_yaw_coarse - 4, best_yaw_coarse + 4.1, 0.1)
    fine_costs = [yaw_cost_in_horizontal_frame(y) for y in fine_yaws]
    best_idx = np.argmin(fine_costs)
    best_yaw_deg = fine_yaws[best_idx]

    if(prnt):
      print(f"\n=== FINAL RESULT (FULL 3D CALIBRATION) ===")
      print(f"Device mounting  → Roll: {roll_deg:+6.3f}° | Pitch: {pitch_deg:+6.3f}°")
      print(f"Remaining yaw offset (horizontal → car forward): {best_yaw_deg:+6.3f}°")

    # ============================= OPTIONAL: Plot best yaw =============================
    if(prnt):
        theta = np.deg2rad(best_yaw_deg)
        c, s = np.cos(theta), np.sin(theta)
        a_fwd = c * a_hx + s * a_hy
        a_lat = -s * a_hx + c * a_hy
        v_fwd = detrend(np.cumsum(a_fwd) * dt)
        v_lat = detrend(np.cumsum(a_lat) * dt)

        plt.figure(figsize=(12,8))
        plt.suptitle(f"Best Alignment → Roll {roll_deg:+.2f}° | Pitch {pitch_deg:+.2f}° | Yaw {best_yaw_deg:+.2f}°")

        plt.subplot(3,1,1)
        plt.plot(time - time[0], v_fwd, label='Forward velocity')
        plt.plot(time - time[0], v_lat, label='Lateral velocity')
        plt.legend(); plt.grid(); plt.ylabel('Velocity (m/s)')

        plt.subplot(3,1,2)
        plt.plot(time - time[0], a_fwd, label='a_forward')
        plt.plot(time - time[0], a_lat, label='a_lateral')
        plt.legend(); plt.grid(); plt.ylabel('Accel (m/s²)')

        plt.subplot(3,1,3)
        plt.plot(yaws, costs, 'b.-', markersize=3)
        plt.axvline(best_yaw_deg, color='r', linestyle='--')
        plt.xlabel('Yaw offset (deg)')
        plt.ylabel('Cost')
        plt.grid()
        plt.tight_layout()
        plt.show()
    
    return roll_deg, pitch_deg, best_yaw_deg



