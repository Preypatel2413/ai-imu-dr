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
"""


def stage1(input_path, output_path, prnt = False):

    INPUT_FILE = input_path
    OUTPUT_FILE = output_path

    df = pd.read_csv(INPUT_FILE)

    cols_to_keep = [
        "GPS Time(ms)", "Lat", "Lon", "GPS Altitude", "GPS Speed", "Map Bearing",
        "Snapped Lat", "Snapped Lon",
        "Accel Sensor Time(ms)", "Ax", "Ay", "Az",
        "Gyro Sensor Time(ms)", "Gx", "Gy", "Gz",
        "Magn Sensor Time(ms)", "Mx", "My", "Mz",
        "Pressure Sensor Time(ms)", "Pressure"
    ]

    df = df[cols_to_keep]
    df.rename(columns={
        'GPS Time(ms)': 'GPST',
        'Accel Sensor Time(ms)': 'AST',
        'Gyro Sensor Time(ms)': 'GST',
        'Magn Sensor Time(ms)': 'MST',
        'Pressure Sensor Time(ms)': 'PST'
    }, inplace=True)

    max_skip_rows = 1000   # scan up to this many rows to find the first aligned timestamps
    sync_time = 100         # max allowed timestamp spread (ms) to treat a row as synced
    
    cols = ["GPST", "AST", "GST", "MST", "PST"]
    limit = min(len(df), max_skip_rows)
    ts = df[cols].iloc[:limit].to_numpy(dtype=float)
    row_max = np.nanmax(ts, axis=1)
    row_min = np.nanmin(ts, axis=1)
    row_range = row_max - row_min
    candidates = np.where(row_range < sync_time)[0]

    if len(candidates) > 0:
        i = int(candidates[0])   
    else:
        i = limit - 1         
        
    df = df[i:].reset_index(drop=True)


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
