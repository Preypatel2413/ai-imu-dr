from utils import stage1_raw, stage2_sync, stage3_orient, stage4_rotate, stage5_pickle, stage1_neo
from os import listdir
from os.path import isfile, join
import os
import subprocess

"""
main_preprocessing.py
=====================

Entry point for the full IMU/GPS preprocessing pipeline.

This script runs all processing stages in sequence:

    Stage 1 -> raw CSV cleanup  
    Stage 2 -> GPS interpolation + IMU synchronization  
    Stage 3 -> device orientation (roll, pitch, yaw) calibration  
    Stage 4 -> rotate IMU measurements into car frame  
    Stage 5 -> serialize final tensors into a .p pickle file

All intermediate CSVs are written inside `raw_data/`, and the final
pickle is written to `data/`.

To run:
    python main_preprocessing.py
"""

def main():
    path = "raw_data/"
    neo_path = "/content/ai-imu-dr/neo-ahrs/src"        #build path of neo-ahrs
    use_neo = True

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    print(onlyfiles)

    if(use_neo):
        os.chdir(neo_path)

    for file in onlyfiles:
        input_file = file
       
        stage1_raw.stage1(input_path= path + input_file, output_path= path + input_file[:-4] +"_data_stage_1.csv")

        input_stage_2 = "_data_stage_1.csv"
        if(use_neo):
            try:
                subprocess.run(["/content/ai-imu-dr/neo-ahrs/src/build/neoahrs", input_file], check=True)
            except subprocess.CalledProcessError as e:
                print("returncode", e.returncode)

            stage1_neo.stage1_neo(input_stage_1= path + input_file[:-4] + "_data_stage_1.csv", input_neo= path + "output_" + input_file, output_file= path + input_file[:-4] + "_data_stage_1_neo.csv")

            input_stage_2 = "_data_stage_1_neo.csv"

        stage2_sync.stage2(input_path= path + input_file[:-4] + input_stage_2, output_path=path + input_file[:-4] + "_data_stage_2.csv")
        roll, pitch, yaw = stage3_orient.stage3(input_path= path + input_file[:-4] + "_data_stage_2.csv")

        stage4_rotate.stage4(input_path= path + input_file[:-4] + "_data_stage_2.csv", output_path= path + input_file[:-4] + "_data_stage_3.csv", yaw_deg= yaw, roll_deg = roll, pitch_deg = pitch)
        stage5_pickle.stage5(input_path= path + input_file[:-4] + "_data_stage_3.csv", output_path= "data/" + input_file[:-4] + ".p")


if __name__ == "__main__":
    main()
