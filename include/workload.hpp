//
// Created by xetql on 10/19/20.
//

#ifndef LBOPT_WORKLOAD_HPP
#define LBOPT_WORKLOAD_HPP
#include <functional>
#include <cmath>
#include <random>
#include <variant>
namespace workload {

    struct Constant {
        double C;
        double operator()(unsigned x) {
            return C;
        }
    };
    struct Sublinear {
        double a = 1, b = 1, c = 0;
        double operator()(unsigned x) {
            return a/(c+b*x);
        }
    };
    struct Linear {
        double a = 1, b = 0;
        double operator()(unsigned x) {
            return a*x+b;
        }
    };
    struct Quadratic {
        double a = 1, b = 0, c = 0;
        double operator()(unsigned x) {
            return a*x*x+b*x+c;
        }
    };
    struct Exp {
        double i;
        double operator()(unsigned x) {
            return std::exp(x / i);
        }
    };
    struct Log {
        double i;
        double operator()(unsigned x) {
            return std::log(i*x);
        }
    };
    struct Sine {
        double i;
        double operator()(unsigned x) {
            return std::sin(i*x);
        }
    };
    struct Uniform {
        std::vector<double> wir{};
        Uniform(unsigned maxi, double min, double max) {
            std::uniform_real_distribution<double> dist {min, max};
            std::random_device rd{};
            std::mt19937 generator { rd() };
            wir.reserve(maxi);
            for(auto i = 0; i < maxi; ++i)
                wir.push_back(dist(generator));
        }
        double operator()(unsigned x) {
            return wir[x];
        }
    };
    struct Normal  {
        std::vector<double> wir{};
        Normal(unsigned maxi, double mu, double stddev){
            std::normal_distribution<double> dist{mu, stddev};
            std::random_device rd{};
            std::mt19937 generator{0};
            wir.reserve(maxi);
            for(auto i = 0; i < maxi; ++i)
                wir.push_back(dist(generator));
        }
        double operator()(unsigned x) {
            return wir[x];
        }
    };
    struct XorY {
        double a, b;
        double operator()(unsigned x) {
            return x % 2 ? a : b;
        }
    };
    namespace {
        using  __WorkloadIncreaseRate = std::variant<Constant, XorY, Sublinear, Linear, Quadratic, Log, Exp, Sine, Uniform, Normal>;
    }
    using  Perturbator                = std::variant<Uniform, Normal, XorY>;

    struct Perturbation {
        __WorkloadIncreaseRate wir;
        Perturbator alg, sys;
        double operator()(unsigned x) {
            return std::visit([x](auto& w){return w(x); }, wir) + std::visit([x](auto& p){return p(x); }, alg) + std::visit([x](auto& p){return p(x); }, sys);
        }
    };

    using  WorkloadIncreaseRate = std::variant<Constant, XorY, Sublinear, Linear, Quadratic, Log, Exp, Sine, Uniform, Normal, Perturbation>;
}

#endif //LBOPT_WORKLOAD_HPP
