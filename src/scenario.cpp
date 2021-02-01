//
// Created by xetql on 8/24/20.
//

#include "scenario.hpp"
std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_static(SimParam p){
    std::vector<bool>   scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    std::vector<double> it_max(p.maxI);
    std::vector<double> it_avg(p.maxI);
    double U = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wmin = p.W0 / p.P;
    double Wavg = Wmax;
    double Tcpu = 0;
    double C = p.C;
    Model model = Balanced;
    Application app{p.P, p.W, p.W0 / p.P, p.W0 / p.P, 1.};
    int last_lb_it = 0;

    for(unsigned int iter = 0; iter < p.maxI; ++iter) {
        auto&[P, W0, Wmax, Wavg, Wmin, d] = app;
        // Compute the iteration
        Tcpu += Wmax;
        it_max[iter] = Wmax;
        it_avg[iter] = Wavg;
        // Measure load imbalance
        U += compute_U(Wmax, Wavg);
        // Store cumulative load imbalance
        imb_time[iter] = U;
        // Apply the workload increase rate function
        update_workloads(iter, last_lb_it, p.deltaImbalance, app);
    }
    save_results(fmt("%s/static-solution.txt", dir), scenario, p, imb_time, it_max, it_avg);
    return {Tcpu, scenario, imb_time};
}
std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_freq(SimParam p, int freq){
    std::vector<bool>   scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    std::vector<double> it_max(p.maxI);
    std::vector<double> it_avg(p.maxI);
    double U = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wmin = p.W0 / p.P;
    double Wavg = Wmax;
    double Tcpu = 0;
    double C = p.C;
    Model model = Balanced;
    Application app{p.P, p.W, p.W0 / p.P, p.W0 / p.P, 1.};
    int last_lb_it = 0;

    for(unsigned int iter = 0; iter < p.maxI; ++iter) {
        auto&[P, W0, Wmax, Wavg, Wmin, d] = app;
        // Based on the previous iteration, should I rebalance ?
        if ((iter % freq) == 0 && iter > 0) { // trigger load balancing
            // Reset cumulative LB
            U = 0;
            // partitioning and mapping -> load balancing
            rebalance(app);
            // Account for the LB cost
            Tcpu += C;
            // This iteration is balanced
            scenario[iter] = true;
            last_lb_it = iter;

        }

        // Compute the iteration
        Tcpu += Wmax;
        it_max[iter] = Wmax;
        it_avg[iter] = Wavg;
        // Measure load imbalance
        U += compute_U(Wmax, Wavg);
        // Store cumulative load imbalance
        imb_time[iter] = U;
        // Apply the workload increase rate function
        update_workloads(iter, last_lb_it, p.deltaImbalance, app);
    }

    save_results(fmt("%s/freq-%d-solution.txt", dir, freq), scenario, p, imb_time, it_max, it_avg);

    return {Tcpu, scenario, imb_time};
}
std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_menon1(SimParam p){
    std::vector<bool> scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    std::vector<double> it_max(p.maxI);
    std::vector<double> it_avg(p.maxI);
    double U = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wavg = Wmax;
    double Wmin = Wmax;
    double Tcpu = 0;
    double C = p.C;
    Model model = Balanced;
    Application app{p.P, p.W, p.W0 / p.P, p.W0 / p.P, 1.};
    int last_lb_it = 0;

    for(int iter = 0; iter < p.maxI; ++iter) {
        auto&[P, W0, Wmax, Wavg, Wmin, d] = app;

        // Based on the previous iteration, should I rebalance ?
        if (U > C) { // trigger load balancing
            // Reset cumulative LB
            U = 0;
            // partitioning and mapping -> load balancing
            // partitioning and mapping -> load balancing
            rebalance(app);
            // Account for the LB cost
            Tcpu += C;
            // This iteration is balanced
            scenario[iter] = true;
            last_lb_it = iter;

        }

        // Compute the iteration
        Tcpu += Wmax;
        it_max[iter] = Wmax;
        it_avg[iter] = Wavg;
        // Measure load imbalance
        U += compute_U(Wmax, Wavg);
        // Store cumulative load imbalance
        imb_time[iter] = U;
        // Apply the workload increase rate function
        update_workloads(iter, last_lb_it, p.deltaImbalance, app);
    }
    save_results(dir+"menon-solution.txt", scenario, p, imb_time, it_max, it_avg);
    return {Tcpu, scenario, imb_time};
}
std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_menon_minus_one(SimParam p){
    std::vector<bool> scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    std::vector<double> it_max(p.maxI);
    std::vector<double> it_avg(p.maxI);
    double U = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wavg = Wmax;
    double Wmin = Wmax;
    double Tcpu = 0;
    double C = p.C;
    Model model = Balanced;
    Application app{p.P, p.W, p.W0 / p.P, p.W0 / p.P, 1.};
    int last_lb_it = 0;

    for(int iter = 0; iter < p.maxI; ++iter) {
        auto&[P, W0, Wmax, Wavg, Wmin, d] = app;
        // Based on the previous iteration, should I rebalance ?
        if (U + (Wmax - Wavg) > C) { // trigger load balancing
            // Reset cumulative LB
            U = 0;
            // partitioning and mapping -> load balancing
            rebalance(app);
            // Account for the LB cost
            Tcpu += C;
            // This iteration is balanced
            scenario[iter] = true;
            last_lb_it = iter;

        }
        // Compute the iteration
        Tcpu += Wmax;
        it_max[iter] = Wmax;
        it_avg[iter] = Wavg;
        // Measure load imbalance
        U += compute_U(Wmax, Wavg);
        // Store cumulative load imbalance
        imb_time[iter] = U;
        // Apply the workload increase rate function
        update_workloads(iter, last_lb_it, p.deltaImbalance, app);

    }

    save_results(dir+"menon-solution-1.txt", scenario, p, imb_time, it_max, it_avg);

    return {Tcpu, scenario, imb_time};
}

std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_bastien(SimParam p){
    std::vector<bool> scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    std::vector<double> it_max(p.maxI);
    std::vector<double> it_avg(p.maxI);
    double U = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wavg = Wmax;
    double Wmin = Wmax;
    double Tcpu = 0;
    double C = p.C;
    Application app{p.P, p.W, p.W0 / p.P, p.W0 / p.P, 1.};
    int last_lb_it = 0;
    double last_Ui = 0;
    for(int iter = 0; iter < p.maxI; ++iter) {
        auto&[P, W0, Wmax, Wavg, Wmin, d] = app;
        auto tau = iter-last_lb_it;
        // Based on the previous iteration, should I rebalance ?
        if ((tau*last_Ui) - U > C) { // trigger load balancing
            // Reset cumulative LB
            U = 0;
            // partitioning and mapping -> load balancing
            // partitioning and mapping -> load balancing
            rebalance(app);
            // Account for the LB cost
            Tcpu += C;
            // This iteration is balanced
            scenario[iter] = true;
            last_lb_it = iter;
        }

        // Compute the iteration
        Tcpu += Wmax;
        it_max[iter] = Wmax;
        it_avg[iter] = Wavg;
        // Measure load imbalance
        U += compute_U(Wmax, Wavg);
        last_Ui = compute_U(Wmax, Wavg);
        // Store cumulative load imbalance
        imb_time[iter] = U;
        // Apply the workload increase rate function
        update_workloads(iter, last_lb_it, p.deltaImbalance, app);
    }
    save_results(dir+"bastien-solution.txt", scenario, p, imb_time, it_max, it_avg);
    return {Tcpu, scenario, imb_time};
}

std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_eff(SimParam p){
    std::vector<bool> scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    double U = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wavg = Wmax;
    double Wmin = Wmax;
    double Tcpu = 0;
    double C = p.C;
    Model model = Balanced;
    Application app{p.P, p.W, p.W0 / p.P, p.W0 / p.P, 1.};
    int last_lb_it = 0;

    for(int iter = 0; iter < p.maxI; ++iter) {
        auto&[P, W0, Wmax, Wavg, Wmin, d] = app;
        U += compute_U(Wmax, Wavg);
        imb_time[iter] = U;
        if (U > C) { // trigger load balancing
            U = 0;
            // partitioning and mapping -> load balancing
            rebalance(app);
            Tcpu += C;
            scenario[iter] = true;
            last_lb_it = iter;

        }
        Tcpu += Wmax;

        update_workloads(iter, last_lb_it, p.deltaImbalance, app);
    }
    return {Tcpu, scenario, imb_time};
}

 std::pair<double, std::vector<double>> compute_tcpu(const std::vector<bool>& scenario, SimParam p, std::string fname){
    std::vector<double> imb_time(p.maxI);
    std::vector<double>   it_max(p.maxI);
    std::vector<double>   it_avg(p.maxI);
    double U = 0;
    double Wmax = p.W0 / p.P;
    double Tcpu = 0;
    double C = p.C;

    Application app{p.P, p.W, p.W0 / p.P, p.W0 / p.P, 1.};
    int last_lb_it = 0;

    for(int iter = 0; iter < p.maxI; ++iter) {
        auto&[P, W0, Wmax, Wavg, Wmin, d] = app;
        bool dec = scenario[iter];

        if (dec) { // trigger load balancing && partitioning and mapping -> load balancing
            U = 0;
            rebalance(app);
            Tcpu += C;
            last_lb_it = iter;
        }

        Tcpu += Wmax;

        it_max[iter] = Wmax;
        it_avg[iter] = Wavg;

        U += compute_U(Wmax, Wavg);

        imb_time[iter] = U;

        update_workloads(iter, last_lb_it, p.deltaImbalance, app);
    }

    save_results(fname, scenario, p, imb_time, it_max, it_avg);
    return {Tcpu, imb_time};
}
