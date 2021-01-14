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

void update_workloads(unsigned iter, unsigned tau, const std::unique_ptr<workload::Function>& deltaImbalance, Application& s){
    double delta = deltaImbalance->operator()(iter, tau);
    auto&[P, W, max, avg, I, d] = s;

    avg = W.at(iter) / P;

    if(1. <= I+d*delta && I+d*delta <= P){
        I += d*delta;
    } else {
        I += d*delta;
        I = std::min((double) P, I);
        I = std::max(1., I);
        d = -d;
    }

    max = I * avg;
}
