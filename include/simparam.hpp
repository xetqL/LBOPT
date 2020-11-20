//
// Created by xetql on 17/02/2020.
//

#ifndef LB_BRANCH_AND_BOUND_SIMPARAM_HPP
#define LB_BRANCH_AND_BOUND_SIMPARAM_HPP

#include <functional>
#include <ostream>
#include "workload.hpp"

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
    const workload::WorkloadIncreaseRate deltaW;
    std::vector<double> mu {};
    std::vector<double> h  {};
    SimParam(double W0, std::vector<double>& W, double C, unsigned maxI, unsigned P, workload::WorkloadIncreaseRate dW):
    W0(W0), W(W), C(C), maxI(maxI), P(P), deltaW(dW)
    {
        mu = W;
        const auto S = mu.size();
        std::for_each(mu.begin(), mu.end(), [P] ( auto& mu) {mu /= P;});
        h.resize(mu.size());
        for(int i = 0; i < S; ++i) {
            h.at(i) = std::accumulate(mu.begin() + i, mu.end(), 0.0);
        }

    }

    friend std::ostream &operator<<(std::ostream &os, const SimParam &param) {
        os << "W0: " << param.W0 << " C: " << param.C << " maxI: " << param.maxI << " P: " << param.P;
        return os;
    }
};

#endif //LB_BRANCH_AND_BOUND_SIMPARAM_HPP
