//
// Created by xetql on 8/24/20.
//

#ifndef LB_BRANCH_AND_BOUND_IO_HPP
#define LB_BRANCH_AND_BOUND_IO_HPP

#include "utils.hpp"
#include "simparam.hpp"
#include <fstream>

void save_results(std::string fname,
                  const std::vector<bool>& scenario,
                  const SimParam& p,
                  const std::vector<double>& imp_time,
                  const std::vector<double>& it_max,
                  const std::vector<double>& it_avg);

#endif //LB_BRANCH_AND_BOUND_IO_HPP
