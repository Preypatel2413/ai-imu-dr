#include <string>
#include <chrono>
#include <iostream>

#include "neoahrs.h"

int main(int argc, char* argv[])
{
  // Check if the user provided a filename
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <neolog_filename>" << std::endl;
    return 1;
  }

  // Input Neolog filename from command-line argument
  std::string neolog_filename = argv[1];
  std::string file_path = "../../raw_data/";

  // Start the timer
  auto start = std::chrono::high_resolution_clock::now();
  std::cout << "Processing \"" << neolog_filename << "\"." << std::endl;

  // Select the output from NeoAHRS
  int output_flags = NEOAHRS_OUTPUT_YAW_COMPASS | NEOAHRS_OUTPUT_YAW_SENSOR_GYRO;

  // Construct input relative path and perform NeoAHRS
  std::string input_file = file_path + neolog_filename;
  neo_output output = neoahrs::neoahrs(input_file, output_flags);

  // Construct relative output path and export to CSV
  std::string output_file = "../../raw_data/output_" + neolog_filename;
  neoahrs::export_csv(output_file, output);

  // Display elapsed time
  std::cout << "NeoAHRS Complete!" << std::endl;
  std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
  std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

  return 0;
}
