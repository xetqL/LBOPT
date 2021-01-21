//
// Created by xetql on 10/19/20.
//

#ifndef LBOPT_WORKLOAD_HPP
#define LBOPT_WORKLOAD_HPP
#include <functional>
#include <cmath>
#include <random>
#include <variant>
#include <iostream>
namespace workload {

    struct Constant {
        double C;
        double operator()(unsigned x) {
            return C;
        }
    };
    struct Sublinear {
        double a,b,c;

        Sublinear(double a, double b, double c) : a(a), b(b), c(c) {}

        double operator()(unsigned x) {
            return a / (b*x + c);
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
            return a*x*x + b*x +c;
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

    struct GaussianPDF {
        double sigma, mu;
        double operator()(unsigned x) {
            auto y = static_cast<double>(x);
            return std::pow((1. / (sigma * std::sqrt(2. * M_PI))),
                            std::exp(-0.5 * std::pow((y - mu) / sigma, 2.) ) );
        }
    };

    struct SymmetricLinear {
        int i;
        double a, b;
        double operator()(unsigned x) {
            return x > i ? -(a*x + b) : a*x+b;
        }
    };

    using  WorkloadIncreaseRate = std::variant<Constant, XorY, Sublinear, Linear, Quadratic, Log, Exp, Sine, Uniform, Normal, Perturbation, GaussianPDF, SymmetricLinear>;
}

#endif //LBOPT_WORKLOAD_HPP
