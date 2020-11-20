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

State::State(Model model, Workload max, Workload avg, Workload min) : model(model), max(max), avg(avg), min(min) {}

void rebalance(State& state) {
    state.model = Balanced;
    state.max   = state.avg;
    state.min   = state.avg;
}

void update_workloads(unsigned int iter, unsigned int P, workload::WorkloadIncreaseRate deltaW, State& s){
    double delta = std::visit([iter](auto& wir){return wir(iter);}, deltaW);
    auto&[model, Wmax, Wavg, Wmin] = s;
    switch(model) {
        case Balanced:
            if(delta > 0) {
                Wmax += delta;
                model = Increasing;
            } else if(delta < 0) {
                Wmin += delta;
                model = Decreasing;
            }
            break;
        case Increasing:
            transfer_bound(Wmax, Wmin, Wavg, delta, std::greater_equal<>(), model);
            break;
        case Decreasing:
            transfer_bound(Wmin, Wmax, Wavg, delta, std::less_equal<>(), model);
            break;
    }
    Wmin = std::max(0.0, Wmin);
    Wmax = std::max(0.0, Wmax);
    Wavg = std::max(0.0, Wavg + delta / P);
}
