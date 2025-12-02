# I have currently included only some required configurable variables here. 

# Part of main_kitti.py
TRAIN_FILTER = 0                            #should be 1 for training process
CROSS_VALIDATION_SEQUENCES = []             #paths to be used as cross validation sequences
TEST_SEQUENCES = []                         #paths to be used as testing sequences


# Part of train_torch_filter.py
SAVE_TO_DRIVE = False                       #to save checkpoint iekfnets.p on google drive (only useful if google drive is mounted)
DRIVE_PATH = "/content/drive/MyDrive/AI_IMU_DR/checkpoints"         #path to checkpoints folder on drive