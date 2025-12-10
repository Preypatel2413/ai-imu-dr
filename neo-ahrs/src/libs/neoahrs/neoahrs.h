// Implementation of the Neo-AHRS architecture.

#ifndef _NEOAHRS_H_
#define _NEOAHRS_H_

#include <vector>

#include "data_types.h"

namespace neoahrs
{
  neo_output neoahrs(std::string &filename, int output_flag);
  neo_output neoahrs(std::vector<neo_log> &data, int output_flag);
  void export_csv(const std::string &filename, const neo_output &output);
}

#endif // neoahrs.h
