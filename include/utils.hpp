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

#define debug(x) std::cout << #x <<"\t"<< x << std::endl;
#define dump(i, x) std::cout << (i) << ": " << (#x) <<"\t" << (x) << std::endl

template<class T> using ptr_t = std::unique_ptr<T>;

template<typename T>
constexpr auto convert(T&& t) {
    if constexpr (std::is_same<std::remove_cv_t<std::remove_reference_t<T>>, std::string>::value) {
        return std::forward<T>(t).c_str();
    } else {
        return std::forward<T>(t);
    }
}

/**
 * printf like formatting for C++ with std::string
 * Original source: https://stackoverflow.com/a/26221725/11722
 */
template<typename ... Args>
std::string stringFormatInternal(const std::string& format, Args&& ... args)
{
    size_t size = snprintf(nullptr, 0, format.c_str(), std::forward<Args>(args) ...) + 1;
    if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    std::unique_ptr<char[]> buf(new char[size]);
    snprintf(buf.get(), size, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size - 1);
}

template<typename ... Args>
std::string fmt(std::string fmt, Args&& ... args) {
    return stringFormatInternal(fmt, convert(std::forward<Args>(args))...);
}

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
    State(Model m, Workload max, Workload avg, Workload min);
    Model model;
    Workload max = 0.0, avg = 0.0, min=0.;
};

struct Application {
    unsigned P;
    std::vector<Workload> W{};
    Workload max, avg;
    double imbalance;
    int d = 1;
};

void rebalance(Application& app);

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


#endif //LB_BRANCH_AND_BOUND_UTILS_HPP
