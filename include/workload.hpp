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
#include <memory>
#include "utils.hpp"

namespace workload {
    struct Function {
        virtual double operator()(unsigned t) const = 0;
        virtual double operator()(unsigned t, unsigned tau) const {
            return this->operator()(t) - this->operator()(tau);
        }
    };
    struct Constant : public Function {
        double C = 0.0;
        explicit Constant(double C) : C(C){}
        double operator()(unsigned x) const override {
            return C;
        }
        virtual double operator()(unsigned t, unsigned tau) const override {
            return this->operator()(t);
        }
    };
    struct Sublinear : public Function {
        double a, b, c;
        Sublinear(double a, double b, double c) : a(a), b(b), c(c) {}
        double operator() (unsigned x) const override {
            return a / (b * x + c);
        }
    };
    struct Sublinear2 : public Function {
        double a, b, c;
        Sublinear2(double a, double b, double c) : a(a), b(b), c(c) {}
        double operator() (unsigned x) const override {
            return a / (b * (1./x) + c);
        }
    };

    struct Linear : public Function {
        double a = 1, b = 0;
        Linear(double a, double b) : a(a), b(b) {}
        double operator()(unsigned x) const override {
            return a*x+b;
        }
    };

    struct Quadratic : public Function {
        double a = 1, b = 0, c = 0;
        Quadratic(double a, double b, double c) : a(a), b(b), c(c) {}
        double operator()(unsigned x) const override {
            return a * x*x + b*x +c;
        }
    };
    struct Exp : public Function {
        double i = 0;
        double operator()(unsigned x) const override {
            return std::exp(x / i);
        }
    };
    struct Log : public Function {
        double a = 0, b = 0;
        Log(double a, double b) : a(a), b(b) {}
        double operator()(unsigned x) const override {
            return std::log(a*x + b);
        }
    };
    struct Sine : public Function {
        double i = 0;
        double operator()(unsigned x) const override {
            return std::sin(i*x);
        }
    };
    struct Uniform : public Function {
        std::vector<double> wir{};
        Uniform(unsigned maxi, double min, double max) {
            std::uniform_real_distribution<double> dist {min, max};
            std::random_device rd{};
            std::mt19937 generator { rd() };
            wir.reserve(maxi);
            for(auto i = 0; i < maxi; ++i)
                wir.push_back(dist(generator));
        }
        double operator()(unsigned x) const override {
            return wir[x];
        }
    };
    struct Normal : public Function  {
        std::vector<double> wir{};
        Normal(unsigned maxi, double mu, double stddev){
            std::normal_distribution<double> dist{mu, stddev};
            std::random_device rd{};
            std::mt19937 generator{ rd() };
            wir.reserve(maxi);
            for(auto i = 0; i < maxi; ++i)
                wir.push_back(dist(generator));
        }
        double operator()(unsigned x) const override {
            return wir[x];
        }
    };
    struct XorY : public Function {
        double a, b;
        double operator()(unsigned x) const override {
            return x % 2 ? a : b;
        }
    };

    using  Perturbator = std::variant<Uniform, Normal, XorY>;
    struct Perturbation : public Function{
        std::unique_ptr<Function> wir;
        Perturbator alg, sys;
        double operator()(unsigned x) const override {
            return wir->operator()(x) + std::visit([x](auto p){return p(x);}, alg) + std::visit([x](auto p){return p(x);}, sys);
        }
    };
    struct GaussianPDF : public Function {
        double sigma = 0., mu = 0.;
        double operator()(unsigned x) const override {
            auto y = static_cast<double>(x);
            return std::pow((1. / (sigma * std::sqrt(2. * M_PI))),
                            std::exp(-0.5 * std::pow((y - mu) / sigma, 2.) ) );
        }
    };
    struct SymmetricLinear : public Function{
        int i = 0;
        double a = 0, b = 0;

        SymmetricLinear(int i, double a, double b) : i(i), a(a), b(b) {}

        double operator()(unsigned x) const override {
            return x > i ? -(a*x + b) : a*x+b;
        }
    };


    template<class Func>
    struct Symmetric : public Function {
        int i = 0;
        Func f;
        Symmetric(int i, Func f) : i(i), f(std::move(f)) {}
        double operator()(unsigned x) const override {
            return x<i? f(x):-f(x);
        }
        virtual double operator() (unsigned t, unsigned tau) const override {
            if(tau)
                return std::abs(this->operator()(t));
            else
                return this->operator()(t);
        }
    };
    template<>
    struct Symmetric<Constant> : public Function {
        int i = 0;
        Constant f;
        Symmetric(int i, Constant f) : i(i), f(std::move(f)) {}
        double operator()(unsigned x) const override {
            return x<i? f(x):-f(x);
        }
        virtual double operator() (unsigned t, unsigned tau) const override {
            if(tau)
                return std::abs(this->operator()(t));
            else
                return this->operator()(t);
        }
    };
    template<class T, unsigned every>
    struct Repeatable : public Function {
        T repeatedFunc;
        template<class ...Arg>
        explicit Repeatable(Arg... args) : repeatedFunc(args...) {}
        double operator()(unsigned x) const override {
            unsigned l = x / every;
            return repeatedFunc(x - l * every);
        }
        virtual double operator()(unsigned t, unsigned tau) const {
            return this->operator()(t);
        }
    };
    using  WorkloadIncreaseRate = std::variant<Constant, XorY, Sublinear, Linear, Quadratic, Log, Exp, Sine, Uniform, Normal, Perturbation, GaussianPDF, SymmetricLinear>;
}

inline double compute_application_workload(double W0, unsigned int i, const ptr_t<workload::Function>& deltaW) {
    auto r = W0;
    for(unsigned k = 0; k < i; ++k) {
        r = std::max(0.0, r + deltaW->operator()(k));
    }
    return r;
}
void update_workloads(unsigned int iter, unsigned int P, const ptr_t<workload::Function>& deltaW , State& s);
void update_workloads(unsigned last_lb_it, unsigned iter, unsigned int P, const ptr_t<workload::Function>& deltaW, State& s);

#endif //LBOPT_WORKLOAD_HPP
