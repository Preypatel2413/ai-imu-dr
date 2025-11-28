import numpy as np
import pandas as pd

"""
stage4.py
=========

Stage 4: rotate IMU measurements from device frame into vehicle (car) frame.

This stage:
  - Loads the time-aligned CSV (from stage2).
  - Builds a rotation from estimated mounting angles (roll, pitch, yaw).
  - Applies the rotation to accelerometer and gyroscope columns (Ax/Ay/Az, Gx/Gy/Gz).
  - Writes a new CSV with rotated IMU columns (other columns preserved).

Notes:
  - Inputs: yaw_deg, roll_deg, pitch_deg are expected in degrees.
  - No filtering or other changes are performed; only a rigid-body rotation is applied.
  - The function writes IMU columns in the output CSV.
"""



def stage4(input_path, output_path, yaw_deg, roll_deg = 0, pitch_deg = 0):
    INPUT = input_path
    OUTPUT = output_path

    roll_dg = -1* roll_deg
    pitch_dg = -1 * pitch_deg
    yaw_dg = yaw_deg


    roll_rad = np.deg2rad(roll_dg)
    pitch_rad = np.deg2rad(pitch_dg)
    yaw_rad = np.deg2rad(yaw_dg)


    imu_acc_cols = ["Ax", "Ay", "Az"]
    imu_gyr_cols = ["Gx", "Gy", "Gz"]
    # ----------------------------

    def Rx(phi):
        c = np.cos(phi); s = np.sin(phi)
        return np.array([[1.0, 0.0, 0.0],
                        [0.0,   c,  -s],
                        [0.0,   s,   c]], dtype=float)

    def Ry(theta):
        c = np.cos(theta); s = np.sin(theta)
        return np.array([[  c, 0.0,   s],
                        [0.0, 1.0,  0.0],
                        [ -s, 0.0,   c]], dtype=float)

    def Rz(psi):
        c = np.cos(psi); s = np.sin(psi)
        return np.array([[  c, -s, 0.0],
                        [  s,  c, 0.0],
                        [0.0, 0.0, 1.0]], dtype=float)

    # Build full device->car rotation matrix.
    # We remove device roll,pitch,yaw by applying inverse rotations:
    def rotation_device_to_car(roll, pitch, yaw):
        return Rz(-yaw) @ Ry(-pitch) @ Rx(-roll) 

    # load CSV
    df = pd.read_csv(INPUT)

    # compute R_total from your estimated angles
    R_total = rotation_device_to_car(roll_rad, pitch_rad, yaw_rad)

    def rotate_columns(df, cols, R):
        """
        Rotate 3-vector columns (x,y,z) using rotation R: v_new = R @ v_old.
        Returns rotated Nx3 numpy array. Does not modify df in place.
        """
        arr = df[cols].to_numpy(dtype=float)   # shape (N,3)
        # apply rotation: for each row v_new = R @ v_old -> arr_rot = (R @ arr.T).T = arr.dot(R.T)
        arr_rot = arr.dot(R.T)
        return arr_rot

    # rotate IMU accel and gyro
    acc_rot = rotate_columns(df, imu_acc_cols, R_total)
    gyr_rot = rotate_columns(df, imu_gyr_cols, R_total)

    # write back
    df[imu_acc_cols] = acc_rot
    df[imu_gyr_cols] = gyr_rot

    df.to_csv(OUTPUT, index=False)
    print("Wrote rotated IMU to:", OUTPUT)
