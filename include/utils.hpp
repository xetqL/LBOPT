//
// Created by xetql on 17/02/2020.
//

#ifndef LB_BRANCH_AND_BOUND_UTILS_HPP
#define LB_BRANCH_AND_BOUND_UTILS_HPP

#include <vector>
#include <ostream>
#include <iostream>
#include <algorithm>
#include <memory>
#include "workload.hpp"

#define debug(x) std::cout << #x <<"\t"<< x << std::endl;
#define dump(i, x) std::cout << (i) << ": " << (#x) <<"\t" << (x) << std::endl

/* Append value and create new vector */
template<class T>
inline std::vector<T> append(const std::vector<T>& prev, T v){
    std::vector<T> next(prev.begin(), prev.end());
    next.push_back(v);
    return next;
}

/* Append value and create new vector */
template<class T>
inline std::vector<T> _append(std::vector<T> prev, T v){
    prev.push_back(v);
    return std::move(prev);
}

/* Get the latest element of the vector */
template<class T>
inline T latest(const std::vector<T>& prev) {
    return *(prev.end() - 1);
}

template<class T>
inline T get(const std::vector<T>& prev, int i) {
    return *(prev.begin() + i);
}

template<class InputIt, class BinaryOp>
void enumerate(InputIt begin, InputIt end, BinaryOp&& fb){
    auto i = 0;
    while(begin != end){
        fb(i, *begin);
        ++begin;
        ++i;
    }
}

template<class FMath>
double sum(int b, int e, FMath& f){
    double r = 0;
    for(int i = b; i < e; ++i) r += f(i);
    return r;
}

/* Function to compute the application workload at a given iteration */
template<class FMath>
inline double compute_application_workload(double W0, unsigned int i, FMath deltaW) {
   auto r = W0;
   for(int k = 0; k < i; ++k) r = std::max(0.0, r + deltaW(k));
   return r;
}

inline double compute_application_workload(double W0, unsigned int i, workload::WorkloadIncreaseRate deltaW) {
    auto r = W0;
    for(int k = 0; k < i; ++k) r = std::max(0.0, r + std::visit([k](auto& wir){return wir(k);}, deltaW));
    return r;
}

template<class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {

    std::for_each(vec.begin(), vec.end(), [&os](auto val){ os << val << " ";});

    return os;
}
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::shared_ptr<T>& pc) {
    os << *pc;
    return os;
}

double compute_U(double ub, double lb);



using Max = double;
using Avg = double;
using Min = double;
using Workload = double;
enum Model {Increasing=0, Decreasing=1, Balanced=2};

Model operator!(Model a);

struct State {
    State(Model model, Workload max, Workload avg, Workload min);

    Model model  = Balanced;
    Workload max = 0.0,
             avg = 0.0,
             min = 0.0;
};

void rebalance(State& state);

template<class Comp>
void transfer_bound(double& ub, double& lb, double mb, double dt, Comp comp, Model& m) {
    if (comp(ub + dt, mb)) {
        ub += dt;
    } else {
        lb = ub + dt;
        ub = mb;
        m  = !m;
    }
}
void update_workloads(unsigned int iter, unsigned int P, workload::WorkloadIncreaseRate deltaW, State& s);

template<class WorkloadIncreaseFunction>
void update_workloads(unsigned int iter, unsigned int P, WorkloadIncreaseFunction deltaW, State& s){
    double delta = deltaW(iter);
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

#endif //LB_BRANCH_AND_BOUND_UTILS_HPP
