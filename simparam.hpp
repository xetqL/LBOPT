//
// Created by xetql on 17/02/2020.
//

#ifndef LB_BRANCH_AND_BOUND_SIMPARAM_HPP
#define LB_BRANCH_AND_BOUND_SIMPARAM_HPP

#include <functional>
#include "utils.hpp"

struct SimParam {
    /* Initial workload */
    const double W0;
    /* Application workload at each iteration */
    const std::vector<double> W;
    /* Load balancing cost */
    const double C;
    /* Number of iteration to simulate */
    const unsigned int maxI;
    /* Number of processors */
    const unsigned int P;
    /* Workload increase load function */
    const std::function<double(int)> deltaW;
};

#endif //LB_BRANCH_AND_BOUND_SIMPARAM_HPP
