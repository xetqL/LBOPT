//
// Created by xetql on 8/24/20.
//

#ifndef LB_BRANCH_AND_BOUND_SCENARIO_HPP
#define LB_BRANCH_AND_BOUND_SCENARIO_HPP

#include <vector>
#include "simparam.hpp"
#include "utils.hpp"
#include "io.hpp"

const std::string dir = std::string("./results/");

using Time      = double;
using Scenario  = std::vector<double>;
using Imbalance = std::vector<double>;

std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_freq(SimParam p, int freq);

std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_menon_minus_one(SimParam p);

std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_menon1(SimParam p);

std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_eff(SimParam p);

template<class LBEfficiencyF>
std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_procassini(SimParam p, double desired_speedup, LBEfficiencyF&& getLBEfficiency){
    std::vector<bool> scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    std::vector<double> it_max(p.maxI);
    std::vector<double> it_avg(p.maxI);
    double t_prime = 0;
    double U = 0;
    double cur_eff = 0, lb_eff = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wavg = Wmax;
    double Wmin = Wmax;
    double Tcpu = 0;
    double C = p.C;
    double t;
    Model model = Balanced;
    State state {model, Wmax, Wavg, Wmin};
    for(int iter = 0; iter < p.maxI; ++iter) {
        auto&[model, Wmax, Wavg, Wmin] = state;

        if (t_prime < t) { // trigger load balancing
            U = 0;
            // partitioning and mapping -> load balancing
            rebalance(state);
            Tcpu += C;
            scenario[iter] = true;
        }

        Tcpu += Wmax;
        it_max[iter] = Wmax;
        it_avg[iter] = Wavg;
        U += compute_U(Wmax, Wavg);
        imb_time[iter] = U;

        cur_eff = Wavg / Wmax;
        lb_eff  = getLBEfficiency();
        t_prime = Wmax * cur_eff/lb_eff + C;
        t = Wmax;

        update_workloads(iter, p.P, p.deltaW, state);
    }

    save_results(dir+"proca-solution.txt", scenario, p, imb_time, it_max, it_avg);

    return {Tcpu, scenario, imb_time};
}

std::pair<double, std::vector<double>> compute_tcpu(const std::vector<bool>& scenario, SimParam p, std::string fname);

#endif //LB_BRANCH_AND_BOUND_SCENARIO_HPP
