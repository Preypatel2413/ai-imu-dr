import numpy as np
import pandas as pd
import os
import pickle

try:
    import torch
except Exception as e:
    raise RuntimeError("This script requires PyTorch (torch). Install with `pip install torch` and re-run.") from e

"""
stage5.py
=========

Stage 5: CSV -> pickle serialization for downstream model code.

This stage:
  - Loads the interpolated CSV produced earlier (stage2/stage4).
  - Validates required columns exist (time, p_gt_*, v_gt_*, ang_gt_*, Ax..Gz).
  - Builds numpy arrays for positions, velocities, attitudes and IMU controls.
  - Reorders IMU into [gx,gy,gz, ax,ay,az] per-sample and stacks into u (N,6).
  - Converts arrays to torch.float32 tensors.
  - Packs everything into a dict and pickles it to disk using Python pickle.
"""


def stage5(input_path, output_path, prnt = False):

    rndm = 0
    pre_csv = input_path
    out_pickle = output_path

    # Name to store inside pickle (if None, it uses basename of pre_csv)
    name_in_pickle = None                # e.g. "2011_09_30_drive_0028_extract" or None to auto-derive

    df = pd.read_csv(pre_csv)
    # df = df[1500:]

    required = [
        "time",
        "p_gt_x","p_gt_y","p_gt_z",
        "v_gt_x","v_gt_y","v_gt_z",
        "ang_gt_roll","ang_gt_pitch","ang_gt_yaw",
        "Ax","Ay","Az","Gx","Gy","Gz"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"pre_pickle.csv is missing required columns: {missing}")

    N = len(df)
    print(f"Rows in CSV: {N}")

    # times are assumed to be seconds relative to t0 (as in your pipeline)
    t_rel = df["time"].to_numpy(dtype=float)              # shape (N,)
    # t0_sec = float(time_offset) / 1000.0                        # convert ms -> seconds
    t0_sec = 0 
    # p_gt: Nx3
    p_gt = df[["p_gt_x","p_gt_y","p_gt_z"]].to_numpy(dtype=float)

    # v_gt: Nx3
    v_gt = df[["v_gt_x","v_gt_y","v_gt_z"]].to_numpy(dtype=float)

    # ang_gt: Nx3 (roll, pitch, yaw) — ensure yaw is in radians already
    ang_gt = df[["ang_gt_roll","ang_gt_pitch","ang_gt_yaw"]].to_numpy(dtype=float)

    # u: reorder into [gx,gy,gz, ax,ay,az] if CSV is Ax,Ay,Az,Gx,Gy,Gz
    # The CSV may have columns Ax..Gz. We'll read them and place into u_order.
    ax = df["Ax"].to_numpy(dtype=float)
    ay = df["Ay"].to_numpy(dtype=float)
    az = df["Az"].to_numpy(dtype=float)
    gx = df["Gx"].to_numpy(dtype=float)
    gy = df["Gy"].to_numpy(dtype=float)
    gz = df["Gz"].to_numpy(dtype=float)

    # Final u ordering (match sample): gyro then accel
    u = np.stack([gx, gy, gz, ax, ay, az], axis=1)   # shape (N,6)

    # Convert numpy arrays to torch tensors (float32)
    ang_gt_t = torch.from_numpy(ang_gt.astype(np.float32))
    p_gt_t   = torch.from_numpy(p_gt.astype(np.float32))
    t_t      = torch.from_numpy(t_rel.astype(np.float32))
    v_gt_t   = torch.from_numpy(v_gt.astype(np.float32))
    u_t      = torch.from_numpy(u.astype(np.float32))

    # Name for pickle
    if name_in_pickle is None:
        base = os.path.splitext(os.path.basename(pre_csv))[0]
        name_in_pickle = base + "_extract"

    # Compose dict matching sample format
    out_dict = {
        "name": name_in_pickle,
        "ang_gt": ang_gt_t,
        "p_gt": p_gt_t,
        "t0": t0_sec,
        "t": t_t,
        "v_gt": v_gt_t,
        "u": u_t
    }

    # Save pickle
    with open(out_pickle, "wb") as f:
        pickle.dump(out_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    # Print summary
    if(prnt):
        print(f"Saved pickle to: {out_pickle}")
        print("Contents summary:")
        print(f" name: {out_dict['name']}")
        print(f" ang_gt: {out_dict['ang_gt'].shape}, dtype={out_dict['ang_gt'].dtype}")
        print(f" p_gt: {out_dict['p_gt'].shape}, dtype={out_dict['p_gt'].dtype}")
        print(f" t0: {out_dict['t0']}")
        print(f" t: {out_dict['t'].shape}, dtype={out_dict['t'].dtype}, first 5: {out_dict['t'][:5]}")
        print(f" v_gt: {out_dict['v_gt'].shape}, dtype={out_dict['v_gt'].dtype}")
        print(f" u: {out_dict['u'].shape}, dtype={out_dict['u'].dtype} (gyro first then accel)")

    # Quick sanity checks
    if out_dict["ang_gt"].shape[0] != N or out_dict["p_gt"].shape[0] != N:
        print("Warning: output tensor lengths do not match input CSV row count.")

