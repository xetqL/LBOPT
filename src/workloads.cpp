//
// Created by xetql on 10/19/20.
//
#include "workload.hpp"
void update_workloads(unsigned int iter, unsigned int P, const ptr_t<workload::Function>& deltaW, State& s){
    double delta = deltaW->operator()(iter);
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

void update_workloads(unsigned last_lb_it, unsigned iter, unsigned int P, const ptr_t<workload::Function>& deltaW, State& s){
    double delta = deltaW->operator()(iter - last_lb_it);
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