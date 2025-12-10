import pandas as pd
import numpy as np

def stage1_neo(input_stage_1, input_neo, output_file):
    df = pd.read_csv(input_stage_1)

    cols = ["AST", "Lat", "Lon", "GPS Altitude", "GPS Speed",
            "Map Bearing", "Snapped Lat", "Snapped Lon"]
    df = df[cols]

    df["GPST_sec"] = df["AST"] / 1000.0
    df["GPST_sec"] -= df["GPST_sec"].min()

    df = df.groupby("GPST_sec", as_index=False).mean()

    t_max = df["GPST_sec"].max()
    uniform_time = np.arange(0, t_max, 0.05)

    gps_cols = ["Lat", "Lon", "GPS Altitude", "GPS Speed",
                "Map Bearing", "Snapped Lat", "Snapped Lon"]

    df_interp = pd.DataFrame({"GPST": uniform_time})

    for col in gps_cols:
        df_interp[col] = np.interp(uniform_time, df["GPST_sec"], df[col])

    df_imu = pd.read_csv(input_neo)
    new_cols = [["AST", "t"], ["Ax", "afx"], ["Ay", "afy"], ["Az", "afz"], ["GST", "t"], ["Gx", "gfx"], ["Gy", "gfy"], ["Gz", "gfz"]]
    imu_time = df_imu["t"].values   # IMU time column
    target_time = uniform_time 

    for new_col, imu_col in new_cols:
        if imu_col == "t":
            df_interp[new_col] = np.interp(target_time, imu_time, df_imu[imu_col])
        else:
            df_interp[new_col] = np.interp(target_time, imu_time, df_imu[imu_col])

    df_interp["GPST"] = df_interp["GPST"]*1000
    df_interp["AST"] = df_interp["AST"]*1000
    df_interp["GST"] = df_interp["GST"]*1000

    # Save final aligned GPS dataset
    df_interp.to_csv(output_file, index=False)
