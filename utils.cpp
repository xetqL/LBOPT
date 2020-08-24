//
// Created by xetql on 17/02/2020.
//

#include "utils.hpp"

double compute_U(double ub, double lb){
    return std::abs(ub - lb);
}
Model operator!(Model a){
    switch(a){
        case Increasing: return Decreasing;
        case Decreasing: return Increasing;
        default:         return Balanced;
    }
}
State::State(Model model, Workload max, Workload avg, Workload min) : model(model), max(max), avg(avg), min(min) {}
void rebalance(State& state) {
    state.model = Balanced;
    state.max   = state.avg;
    state.min   = state.avg;
}