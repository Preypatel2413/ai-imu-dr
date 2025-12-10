#ifndef _MADMADGWICK_H_
#define _MADMADGWICK_H_

#include <array>
#include <vector>

#include "data_types.h"

typedef struct
{
  neo_quat q;
  neo_euler ypr;
} madgwick_t;

std::vector<madgwick_t> madmadgwick(
            const std::vector<neo_vec3>& a,
            const std::vector<neo_vec3>& g,
            const std::vector<double>& t,
            double beta);

std::vector<neo_quat> get_quat(const std::vector<madgwick_t> &m);
std::vector<neo_euler> get_euler(const std::vector<madgwick_t> &m);

#endif // madmadgwick.h
