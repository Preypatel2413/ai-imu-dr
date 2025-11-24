from utils import stage1_raw, stage2_sync, stage3_orient, stage4_rotate, stage5_pickle

"""
main_preprocessing.py
=====================

Entry point for the full IMU/GPS preprocessing pipeline.

This script runs all processing stages in sequence:

    Stage 1 → raw CSV cleanup  
    Stage 2 → GPS interpolation + IMU synchronization  
    Stage 3 → device orientation (roll, pitch, yaw) calibration  
    Stage 4 → rotate IMU measurements into car frame  
    Stage 5 → serialize final tensors into a .p pickle file

All intermediate CSVs are written inside `raw_data/`, and the final
pickle is written to `data/`.

To run:
    python main_preprocessing.py
"""

def main():

    path = "raw_data/"
    input_file = "input_data.csv"

    stage1_raw.stage1(input_path= path + input_file, output_path= path + "data_stage_1.csv")
    stage2_sync.stage2(input_path= path + "data_stage_1.csv", output_path=path + "data_stage_2.csv")
    roll, pitch, yaw = stage3_orient.stage3(input_path= path + "data_stage_2.csv")

    stage4_rotate.stage4(input_path= path + "data_stage_2.csv", output_path= path + "data_stage_3.csv", yaw_deg= yaw, roll_deg = roll, pitch_deg = pitch)
    stage5_pickle.stage5(input_path= path + "data_stage_3.csv", output_path= "data/2011_09_30_drive_0028_extract.p")


if __name__ == "__main__":
    main()
