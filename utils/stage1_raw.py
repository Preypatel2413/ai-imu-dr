import numpy as np
import pandas as pd

"""
stage1.py
=========
Stage 1 preprocessing for IMU/GPS pipeline.

This stage:
  - loads raw CSV exported from device,
  - keeps only the columns we care about,
  - optionally skips initial rows (to align initial timestamps),
  - renames timestamp columns to short names (GPST, AST, GST, MST, PST),
  - converts timestamps to a relative timebase by subtracting a rounded
    `time_offset` (ms),
  - writes the reduced CSV to disk.

The function `stage1` is intentionally small and deterministic — it does
not attempt to interpolate or resample. Later stages handle interpolation,
synchronization and filtering.
"""


def stage1(input_path, output_path, prnt = False):

    INPUT_FILE = input_path
    OUTPUT_FILE = output_path

    df = pd.read_csv(INPUT_FILE, sep=",", engine= "python")

    cols_to_keep = [
        "GPS Time(ms)", "Lat", "Lon", "GPS Altitude", "GPS Speed", "Map Bearing",
        "Snapped Lat", "Snapped Lon",
        "Accel Sensor Time(ms)", "Ax", "Ay", "Az",
        "Gyro Sensor Time(ms)", "Gx", "Gy", "Gz",
        "Magn Sensor Time(ms)", "Mx", "My", "Mz",
        "Pressure Sensor Time(ms)", "Pressure"
    ]

    df = df[cols_to_keep]
    df = df[200:]                   # Skip initial rows to avoid mismatched initial timing bursts.

    df.rename(columns={
        'GPS Time(ms)': 'GPST',
        'Accel Sensor Time(ms)': 'AST',
        'Gyro Sensor Time(ms)': 'GST',
        'Magn Sensor Time(ms)': 'MST',
        'Pressure Sensor Time(ms)': 'PST'
    }, inplace=True)

    # All timestamps are epoch-based, so subtract the minimum (rounded to nearest second)
    mod = 1000
    time_offset = (min(df['GPST'].iloc[0], df['AST'].iloc[0], df['GST'].iloc[0], df['MST'].iloc[0], df['PST'].iloc[0])//mod )* mod

    print("Time offset : ", time_offset)
    df['GPST'] = df['GPST'] - time_offset
    df['AST'] = df['AST'] - time_offset
    df['GST'] = df['GST'] - time_offset
    df['MST'] = df['MST'] - time_offset
    df['PST'] = df['PST'] - time_offset

    df.to_csv(OUTPUT_FILE, index=False)
    if(prnt):
        print(f"Saved necessary columns -> {OUTPUT_FILE}")
        print(df.head(10))
