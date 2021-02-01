//
// Created by xetql on 17/02/2020.
//

#include "utils.hpp"

double compute_U(double ub, double lb){
    return std::abs(ub - lb);
}

Model operator!(Model a){
    switch(a) {
        case Increasing: return Decreasing;
        case Decreasing: return Increasing;
        default:         return Balanced;
    }
}

State::State(Model m, Workload max, Workload avg, Workload min) : model(model), max(max), avg(avg), min(min) {}

void rebalance(Application& app) {
    app.max = app.avg;
    app.imbalance = 1.;
}

