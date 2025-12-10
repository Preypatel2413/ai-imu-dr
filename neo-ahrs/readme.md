# neo-ahrs

C++ implementation of AHRS for NeoMatrix.

## Contents

Complete documentation of all functions is yet to be completed. However, the current version includes one-liner documentation to make the review intuitive for Tatyana. This file provides a brief overview of the files involved and their roles in Neo-AHRS.

### Custom Files

0. [data_types](src/libs/neoahrs/data_types.h) - This file consists of the declaration of custom datatypes used throughout the codebase to make interface intuitive. I tried my best to make all important function accept and return those datatype. There is also a datatype for error. However, as of now I have not done error handling which is on my todo list.
1. [neoahrs](src/libs/neoahrs/neoahrs.h) - Implements the architecture of NeoAHRS.
3. [main](src/main.cpp) - Calls functions from [neoahrs](src/libs/neoahrs/neoahrs.h) and performs the computations.
4. [butterworth](src/libs/butterworth/butterworth.h) - Implements a second-order Butterworth filter—both causal and non-causal versions.
5. [madgwick](src/libs/madgwick/madgwick.h) - Standard implementation of the Madgwick filter.
6. [madmadgwick](src/libs/madgwick/madmadgwick.h) - Custom non-causal version of the Madgwick filter that addresses incorrect quaternion initialization, eliminating transient period errors.
7. [magcal](src/libs/magcal/magcal.h) - Magnetometer calibration algorithm that selects one of three fits based on residual error.
8. [attitude](src/libs/attitude/attitude.h) - General utilities for orientation manipulation. Supports Euler angles, quaternions, and rotation matrices.

### Third-Party Libraries

1. [eigen](src/libs/Eigen/) - C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
2. [rapidcsv](src/libs/rapidcsv/rapidcsv.h) - Easy-to-use C++ CSV parser library.

### How to run?

1. Make a directory called `data` just outside this repository directory `neo-ahrs`. It is important that the input Neologs .CSV files are inside the directory `data`. The CSV files can be found [here](https://drive.google.com/drive/folders/1h5J8DNf7p7SMFflW0tmAsimXGsO8u74j?usp=drive_link). If the link is broken, please contact Tatyana.
2. Compilation: `make`
3. The software accepts filename (assuming is is located at directory `data` outside `neo-ahrs`) as argument.
   ```
    ./build/neoahrs.exe 02_05_2025_10_30_11_95e41a68.neolog.csv
   ```
4. The output are stored inside the directory [data](src) with the prefix `output_` appended into the original filename.

### Variables

The ultimate output of the `neoahrs` is hold by the data structure ``neo_output`` with following elements.

| Variable | Unit  | Remark                                            |
| ---------| ----- | --------------------------------------------------|
| `t `     | s     | Master timestamp                                  |
| `ym`     | deg   | Magnetometer yaw                                  |
| `ys`     | deg   | Gyro-accel heading profile starting from 0 deg    |
| `ps`     | deg   | Sensor pitch wrt. NED                             |
| `rs`     | deg   | Sensor roll wrt. NED                              |
| `alx`    | m/s²  | X-axis of gravity compensated acceleration in NED |
| `aly`    | m/s²  | Y-axis of gravity compensated acceleration in NED |
| `alz`    | m/s²  | Z-axis of gravity compensated acceleration in NED |
| `afx`    | m/s²  | X-axis of virtual flat accelerometer              |
| `afy`    | m/s²  | Y-axis of virtual flat accelerometer              |
| `afz`    | m/s²  | Z-axis of virtual flat accelerometer              |
| `gfx`    | rad/s | X-axis of virtual flat gyroscope                  |
| `gfy`    | rad/s | Y-axis of virtual flat gyroscope                  |
| `gfz`    | rad/s | Z-axis of virtual flat gyroscope                  |
