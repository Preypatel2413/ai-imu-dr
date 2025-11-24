import numpy as np
import pandas as pd

try:
    from pyproj import Transformer # type: ignore
    HAVE_PYPROJ = True
except Exception:
    HAVE_PYPROJ = False

"""
stage2.py
=========

Stage 2 preprocessing for the IMU/GPS pipeline.

This stage:
  - Loads the cleaned CSV from stage1.
  - Chooses a master timebase (AST → GST → GPST).
  - Converts GPS lat/lon to local ENU coordinates.
  - Interpolates GPS position, altitude, speed, and bearing to IMU timestamps.
  - Computes ground-truth velocity (vx, vy, vz) and yaw from GPS bearing.
  - Synchronizes raw accel/gyro data to the master timeline using a ±window
    mean (handles sensor batching).
  - Builds a unified, time-aligned dataset containing:
        time (s), p_gt_*, v_gt_*, ang_gt_*, Ax/Ay/Az, Gx/Gy/Gz.
  - Detects large time gaps and fills them with interpolated samples.
  - Writes the synchronized and gap-filled CSV (data_stage_2).

Output is fully aligned IMU + GPS ground truth.
"""


def stage2(input_path, output_path, prnt=False):
    INPUT = input_path
    OUTPUT = output_path

    # Which IMU timestamp to use as the master timebase
    PREFERRED_IMU_TIME = ["AST", "GST", "GPST"]  # accel -> gyro -> gps fallback
    USE_ENU = True
    ROLL_PITCH_NOISE_STD = 0.0
    MAX_GAP_S = 0.5

    # Window for accel/gyro synchronization (ms)
    WINDOW_HALF_MS = 8   # +-15 ms is safe and works perfectly on all Android phones

    # requested column groups / names
    gt_cols = ["p_gt_x","p_gt_y","p_gt_z","v_gt_x","v_gt_y","v_gt_z",
            "ang_gt_roll","ang_gt_pitch","ang_gt_yaw"]
    imu_cols = ["Ax","Ay","Az","Gx","Gy","Gz"]
    time_col = "time"

    df = pd.read_csv(INPUT)

    # Choose master timebase (prefers accel -> gyro -> gps)
    imu_time = None
    for cand in PREFERRED_IMU_TIME:
        if cand in df.columns:
            imu_time = cand
            break
    if imu_time is None:
        raise ValueError(f"None of {PREFERRED_IMU_TIME} found. Columns: {df.columns.tolist()}")

    t_target_ms = df[imu_time].values.astype(np.float64)   # master timestamps in ms

    # ============ GPS -> ENU + interpolation ============
    if "GPST" not in df.columns:
        df["GPST"] = df[imu_time]

    gps_df = df[["GPST", "Lat", "Lon", "GPS Altitude", "GPS Speed", "Map Bearing"]].copy()
    gps_df = gps_df.dropna(subset=["GPST"])
    gps_df = gps_df.groupby("GPST", as_index=False).first()     # when multiple rows have same GPST we keep the first (groupby.first)
    gps_df = gps_df.sort_values("GPST")

    t_gps = gps_df["GPST"].values.astype(np.float64)
    lat_gps = gps_df["Lat"].astype(float).values
    lon_gps = gps_df["Lon"].astype(float).values
    alt_gps = gps_df["GPS Altitude"].astype(float).values if "GPS Altitude" in gps_df.columns else np.full_like(t_gps, np.nan)
    speed_gps = gps_df["GPS Speed"].astype(float).values if "GPS Speed" in gps_df.columns else np.full_like(t_gps, np.nan)
    bearing_gps = gps_df["Map Bearing"].astype(float).values if "Map Bearing" in gps_df.columns else np.full_like(t_gps, np.nan)

    def latlon_to_enu(lats, lons, lat0=None, lon0=None):
        if lat0 is None: lat0 = float(np.mean(lats))
        if lon0 is None: lon0 = float(np.mean(lons))
        if HAVE_PYPROJ and USE_ENU:
            proj_str = f"+proj=aeqd +R=6378137 +lat_0={lat0} +lon_0={lon0}"
            transformer = Transformer.from_crs("epsg:4326", proj_str, always_xy=True)
            xs, ys = transformer.transform(lons.tolist(), lats.tolist())
            return np.array(xs), np.array(ys), lat0, lon0
        else:
            R = 6378137.0
            lat0_rad = np.deg2rad(lat0)
            dlat = np.deg2rad(lats - lat0)
            dlon = np.deg2rad(lons - lon0)
            x = R * dlon * np.cos(lat0_rad)
            y = R * dlat
            return x, y, lat0, lon0

    if len(lat_gps) == 0:
        raise ValueError("No GPS samples found!")

    x_gps, y_gps, lat0, lon0 = latlon_to_enu(lat_gps, lon_gps)

    order = np.argsort(t_gps)
    t_gps = t_gps[order]
    x_gps = x_gps[order]; y_gps = y_gps[order]
    alt_gps = alt_gps[order]; speed_gps = speed_gps[order]; bearing_gps = bearing_gps[order]

    def circular_interp_deg(t_src, angles_deg, t_tgt):
        ang_rad = np.deg2rad(np.array(angles_deg, dtype=float))
        sin_i = np.interp(t_tgt, t_src, np.sin(ang_rad), left=np.nan, right=np.nan)
        cos_i = np.interp(t_tgt, t_src, np.cos(ang_rad), left=np.nan, right=np.nan)
        ang_i = np.rad2deg(np.arctan2(sin_i, cos_i))
        return (ang_i + 360.0) % 360.0

    # Interpolate GPS fields
    x_target    = np.interp(t_target_ms, t_gps, x_gps, left=np.nan, right=np.nan)
    y_target    = np.interp(t_target_ms, t_gps, y_gps, left=np.nan, right=np.nan)
    alt_target  = np.interp(t_target_ms, t_gps, alt_gps, left=np.nan, right=np.nan)
    speed_target= np.interp(t_target_ms, t_gps, speed_gps, left=np.nan, right=np.nan)
    bearing_target = circular_interp_deg(t_gps, bearing_gps, t_target_ms)

    # Fill leading/trailing with nearest GPS
    if t_target_ms[0] < t_gps[0]:
        first = np.searchsorted(t_target_ms, t_gps[0])
        for arr, val in zip([x_target,y_target,alt_target,speed_target,bearing_target],
                            [x_gps[0],y_gps[0],alt_gps[0],speed_gps[0],bearing_gps[0]]):
            arr[:first] = val
    if t_target_ms[-1] > t_gps[-1]:
        last = np.searchsorted(t_target_ms, t_gps[-1], side="right")
        for arr, val in zip([x_target,y_target,alt_target,speed_target,bearing_target],
                            [x_gps[-1],y_gps[-1],alt_gps[-1],speed_gps[-1],bearing_gps[-1]]):
            arr[last:] = val

    # ============ VELOCITY FROM SPEED + BEARING ============
    bearing_rad = np.deg2rad(bearing_target)
    vx = speed_target * np.cos(bearing_rad)
    vy = speed_target * np.sin(bearing_rad)

    t_s = t_target_ms / 1000.0
    vz = np.zeros_like(alt_target)
    if len(t_s) >= 2:
        dt = np.diff(t_s)
        dalt = np.diff(alt_target)
        valid = dt > 1e-9
        vz_mid = np.zeros_like(dalt)
        vz_mid[valid] = dalt[valid] / dt[valid]
        vz[0] = vz_mid[0]
        vz[-1] = vz_mid[-1]
        if len(vz_mid) > 1:
            vz[1:-1] = 0.5 * (vz_mid[:-1] + vz_mid[1:])

    # ============ BUILD OUTPUT DATAFRAME ============
    df_out = pd.DataFrame()
    df_out[time_col] = t_s  # seconds

    df_out["p_gt_x"] = x_target
    df_out["p_gt_y"] = y_target
    df_out["p_gt_z"] = alt_target
    df_out["v_gt_x"] = vx
    df_out["v_gt_y"] = vy
    df_out["v_gt_z"] = vz

    yaw_rad = (np.deg2rad(bearing_target) + np.pi) % (2*np.pi) - np.pi       # yaw is derived from GPS bearing
    df_out["ang_gt_roll"]  = 0.0
    df_out["ang_gt_pitch"] = 0.0
    df_out["ang_gt_yaw"]   = yaw_rad

    if ROLL_PITCH_NOISE_STD > 0:
        rng = np.random.default_rng(42)
        df_out["ang_gt_roll"]  += rng.normal(0, ROLL_PITCH_NOISE_STD, len(df_out))
        df_out["ang_gt_pitch"] += rng.normal(0, ROLL_PITCH_NOISE_STD, len(df_out))


    # ============ ROBUST IMU SYNCHRONIZATION ============
    if(prnt):
        print(f"Synchronizing IMU using ±{WINDOW_HALF_MS} ms window averaging...")

    imu_parts = []
    if all(c in df.columns for c in ["AST", "Ax", "Ay", "Az"]):
        accel = df[["AST", "Ax", "Ay", "Az"]].rename(columns={"AST": "t_ms"})
        accel["sensor"] = "accel"
        imu_parts.append(accel)
    if all(c in df.columns for c in ["GST", "Gx", "Gy", "Gz"]):
        gyro = df[["GST", "Gx", "Gy", "Gz"]].rename(columns={"GST": "t_ms"})
        gyro["sensor"] = "gyro"
        imu_parts.append(gyro)

    if not imu_parts:
        raise ValueError("No accel or gyro data found!")

    imu_all = pd.concat(imu_parts, ignore_index=True)
    imu_all = imu_all.sort_values("t_ms").reset_index(drop=True)

    def window_mean_sync(t_native, values, t_target, half_win):
        out = np.full((len(t_target), values.shape[1]), np.nan)
        for i, t in enumerate(t_target):
            mask = (t_native >= t - half_win) & (t_native <= t + half_win)
            if mask.any():
                out[i] = np.mean(values[mask], axis=0)
            else:
                # fallback nearest
                idx = np.argmin(np.abs(t_native - t))
                out[i] = values[idx]
        return out

    # Extract raw data
    accel_mask = imu_all["sensor"] == "accel"
    gyro_mask  = imu_all["sensor"] == "gyro"

    accel_t = imu_all.loc[accel_mask, "t_ms"].values
    gyro_t  = imu_all.loc[gyro_mask,  "t_ms"].values

    accel_val = imu_all.loc[accel_mask, ["Ax","Ay","Az"]].values.astype(float)
    gyro_val  = imu_all.loc[gyro_mask,  ["Gx","Gy","Gz"]].values.astype(float)

    synced_accel = window_mean_sync(accel_t, accel_val, t_target_ms, WINDOW_HALF_MS) if len(accel_val) else np.full((len(t_target_ms),3), np.nan)
    synced_gyro  = window_mean_sync(gyro_t,  gyro_val,  t_target_ms, WINDOW_HALF_MS) if len(gyro_val)  else np.full((len(t_target_ms),3), np.nan)

    df_out["Ax"] = synced_accel[:,0]
    df_out["Ay"] = synced_accel[:,1]
    df_out["Az"] = synced_accel[:,2]
    df_out["Gx"] = synced_gyro[:,0]
    df_out["Gy"] = synced_gyro[:,1]
    df_out["Gz"] = synced_gyro[:,2]

    if(prnt):
        print(f"   Accel samples: {len(accel_t)}, Gyro samples: {len(gyro_t)}, Synced to {len(t_target_ms)} timestamps")

    # ============ FINAL CLEANUP & GAP FILLING  ============
    final_cols = [time_col] + gt_cols + imu_cols
    for c in final_cols:
        if c not in df_out.columns:
            df_out[c] = np.nan
    df_out = df_out[final_cols]

    # Gap filling
    t = df_out[time_col].values
    if len(t) == 0:
        raise ValueError("No samples!")

    _, uniq = np.unique(t, return_index=True)
    if len(uniq) != len(t):
        df_out = df_out.iloc[sorted(uniq)].reset_index(drop=True)
        t = df_out[time_col].values

    new_times = [t[0]]
    for i in range(len(t)-1):
        dt = t[i+1] - t[i]
        if dt > MAX_GAP_S:
            inter = np.arange(t[i] + MAX_GAP_S, t[i+1], MAX_GAP_S)
            new_times.extend(inter.tolist())
        new_times.append(t[i+1])
    new_times = np.array(sorted(set(new_times)))

    if len(new_times) == len(t) and np.allclose(new_times, t):
        df_out.to_csv(OUTPUT, index=False)
        print("Saved (no resampling needed):", OUTPUT)
    else:
        df_interp = df_out.set_index(time_col)
        df_new = pd.DataFrame(index=new_times)

        # circular yaw
        sin_y = np.sin(df_interp["ang_gt_yaw"])
        cos_y = np.cos(df_interp["ang_gt_yaw"])
        sin_i = pd.Series(sin_y, index=df_interp.index).reindex(new_times).interpolate('linear').ffill().bfill()
        cos_i = pd.Series(cos_y, index=df_interp.index).reindex(new_times).interpolate('linear').ffill().bfill()
        df_new["ang_gt_yaw"] = np.arctan2(sin_i.values, cos_i.values)

        for col in df_interp.columns:
            if col == "ang_gt_yaw": continue
            s = df_interp[col].reindex(new_times)
            df_new[col] = s.interpolate('linear').ffill().bfill().values

        df_new = df_new.reset_index().rename(columns={"index": time_col})
        for c in final_cols:
            if c not in df_new.columns:
                df_new[c] = np.nan
        df_new = df_new[final_cols]
        df_new.to_csv(OUTPUT, index=False)
        print(f"Resampled & saved: {OUTPUT} ({len(t)} → {len(new_times)} rows)")

    if(prnt):
        print("Final columns:", df_new.columns.tolist() if 'df_new' in locals() else df_out.columns.tolist())
        print(pd.read_csv(OUTPUT, nrows=8))
