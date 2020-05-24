//
// Created by xetql on 17/02/2020.
//

#ifndef LB_BRANCH_AND_BOUND_UTILS_HPP
#define LB_BRANCH_AND_BOUND_UTILS_HPP

#include <vector>
#include <ostream>
#include <algorithm>
#include <memory>

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
    for(int i = b; i < e; ++i) r = std::max(0.0, r+f(i));
    return r;
}

/* Function to compute the application workload at a given iteration */
template<class FMath>
inline double _W(double W0, unsigned int i, FMath deltaW) {
    return W0 + sum(0, i, deltaW);
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
#endif //LB_BRANCH_AND_BOUND_UTILS_HPP
