//
// Created by xetql on 10/19/20.
//
#include "workload.hpp"

void update_workloads(unsigned iter, unsigned tau, const ptr_t<workload::Function>& deltaImbalance, Application& s){
    double delta = deltaImbalance->operator()(iter - tau);
    auto&[P, W, max, avg, I, d] = s;
    avg = W.at(iter) / P;

    I += d * delta;
    I = std::min(I, P*avg);
    I = std::max(1., I);
    max = I * avg;
}
