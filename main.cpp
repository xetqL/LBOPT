#include <utility>
#include <iostream>
#include <queue>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <zupply.hpp>
#include <fstream>
#include <set>
#include <chrono>
#include <random>

#include "lbnode.hpp"

using namespace zz;

double  __Bfitness(std::vector<bool> s, SimParam p)
{
    double Tcpu=0.;  double Wmin=p.W0/p.P; double Wmax=p.W0/p.P;
    for(int i=0;i<p.maxI;i++){
        if(s[i]) {  // load balancing
            Wmin=((p.P-1)*Wmin + Wmax)/p.P;
            Wmax=Wmin;
            Tcpu+=p.C;
        }
        Tcpu += Wmax;  // add the time of the most loaded proc
        Wmax += p.deltaW(i);
    }
    return Tcpu;
}

std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_freq(SimParam p, int freq){
    std::vector<bool> scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    double U = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wavg = Wmax;
    double Tcpu = 0;
    double C = p.C;
    for(int iter = 0; iter < p.maxI; ++iter) {
        // Based on the previous iteration, should I rebalance ?
        if ((iter % freq) == 0 && iter > 0) { // trigger load balancing
            // Reset cumulative LB
            U = 0;
            // partitioning and mapping -> load balancing
            Wmax = ((P-1) * Wavg + Wmax) / P;
            Wavg = Wmax;
            // Account for the LB cost
            Tcpu += C;
            // This iteration is balanced
            scenario[iter] = true;
        }

        // Compute the iteration
        Tcpu += Wmax;
        // Measure load imbalance
        U += Wmax - Wavg;
        // Store cumulative load imbalance
        imb_time[iter] = U;
        // Apply the workload increase rate function
        double delta = p.deltaW(iter);
        Wavg = std::max(0.0, Wavg + delta / P);
        Wmax = std::max(0.0, Wmax + delta);
        if(Wmax < Wavg) {
            Wmax = Wavg;
        }
    }
    return {Tcpu, scenario, imb_time};
}

std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_menon1(SimParam p){
    std::vector<bool> scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    double U = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wavg = Wmax;
    double Tcpu = 0;
    double C = p.C;
    for(int iter = 0; iter < p.maxI; ++iter) {
        // Based on the previous iteration, should I rebalance ?
        if (U > C) { // trigger load balancing
            // Reset cumulative LB
            U = 0;
            // partitioning and mapping -> load balancing
            Wmax = ((P-1) * Wavg + Wmax) / P;
            Wavg = Wmax;
            // Account for the LB cost
            Tcpu += C;
            // This iteration is balanced
            scenario[iter] = true;
        }

        // Compute the iteration
        Tcpu += Wmax;
        // Measure load imbalance
        U += Wmax - Wavg;
        // Store cumulative load imbalance
        imb_time[iter] = U;
        // Apply the workload increase rate function
        double delta = p.deltaW(iter);
        Wavg = std::max(0.0, Wavg + delta / P);
        Wmax = std::max(0.0, Wmax + delta);
        if(Wmax < Wavg) {
            Wmax = Wavg;
        }
    }
    return {Tcpu, scenario, imb_time};
}

std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_menon_minus_one(SimParam p){
    std::vector<bool> scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    double U = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wavg = Wmax;
    double Tcpu = 0;
    double C = p.C;
    double prev_dw = 0.0;
    for(int iter = 0; iter < p.maxI; ++iter) {
        // Based on the previous iteration, should I rebalance ?
        auto worse_than_before = (Wmax-Wavg) < prev_dw;
        if (U+(Wmax - Wavg) > C) { // trigger load balancing
            // Reset cumulative LB
            U = 0;
            // partitioning and mapping -> load balancing
            Wmax = ((P-1) * Wavg + Wmax) / P;
            Wavg = Wmax;
            // Account for the LB cost
            Tcpu += C;
            // This iteration is balanced
            scenario[iter] = true;
        }

        // Compute the iteration
        Tcpu += Wmax;
        // Measure load imbalance
        U += Wmax - Wavg;
        prev_dw = Wmax-Wavg;
        // Store cumulative load imbalance
        imb_time[iter] = U;
        // Apply the workload increase rate function
        double delta = p.deltaW(iter);
        Wavg = std::max(0.0, Wavg + delta / P);
        Wmax = std::max(0.0, Wmax + delta);
        if(Wmax < Wavg) {
            Wmax = Wavg;
        }
    }
    return {Tcpu, scenario, imb_time};
}

std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_eff(SimParam p){
    std::vector<bool> scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    double U = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wavg = Wmax;
    double Tcpu = 0;
    double C = p.C;
    for(int iter = 0; iter < p.maxI; ++iter) {
        U += Wmax - Wavg;
        imb_time[iter] = U;
        if (U > C) { // trigger load balancing
            U = 0;
            Wmax = ((P-1) * Wavg + Wmax) / P;
            Wavg = Wmax;
            Tcpu += C;
            scenario[iter] = true;
        }
        Tcpu += Wmax;
        double delta = p.deltaW(iter);
        Wavg = std::max(0.0, Wavg + delta / P);
        Wmax = std::max(0.0, Wmax + delta);
        if(Wmax < Wavg) {
            Wmax = Wavg;
        }
    }
    return {Tcpu, scenario, imb_time};
}

std::pair<double, std::vector<double>> compute_tcpu(const std::vector<bool>& scenario, SimParam p){
    std::vector<double> imb_time(p.maxI);
    double U = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wavg = Wmax;
    double Tcpu = 0;
    double C = p.C;
    for(int iter = 0; iter < p.maxI; ++iter) {
        bool dec = scenario[iter];

        if (dec) { // trigger load balancing
            U = 0;
            Wmax = ((P-1) * Wavg + Wmax) / P;
            Wavg = Wmax;
            Tcpu += C;
        }

        Tcpu += Wmax;

        U += Wmax - Wavg;
        imb_time[iter] = U;

        double delta = p.deltaW(iter);
        Wavg = std::max(0.0, Wavg + delta / P);
        Wmax = std::max(0.0, Wmax + delta);
        if(Wmax < Wavg) {
            Wmax = Wavg;
        }
    }
    return {Tcpu, imb_time};
}
template<class LBEfficiencyF>
std::tuple<double, std::vector<bool>, std::vector<double>> create_scenario_procassini(SimParam p, double desired_speedup, LBEfficiencyF&& getLBEfficiency){
    std::vector<bool> scenario(p.maxI);
    std::vector<double> imb_time(p.maxI);
    double t_prime = 0;
    double U = 0;
    double cur_eff = 0, lb_eff = 0;
    double P = p.P;
    double Wmax = p.W0 / p.P;
    double Wavg = Wmax;
    double Tcpu = 0;
    double C = p.C;
    for(int iter = 0; iter < p.maxI; ++iter) {

        Tcpu += Wmax;
        U += Wmax - Wavg;
        imb_time[iter] = U;

        cur_eff = Wavg / Wmax;
        lb_eff  = getLBEfficiency();
        t_prime = Wmax * cur_eff/lb_eff + C;
        auto t = Wmax;

        double delta = p.deltaW(iter);
        Wavg = std::max(0.0, Wavg + delta / P);
        Wmax = std::max(0.0, Wmax + delta);

        if(Wmax < Wavg) {
            Wmax = Wavg;
        }

        if (t_prime < t) { // trigger load balancing
            U = 0;
            Wmax = ((P-1) * Wavg + Wmax) / P;
            Wavg = Wmax;
            Tcpu += C;
            scenario[iter] = true;
        }
    }
    return {Tcpu, scenario, imb_time};
}

int main(int argc, char** argv) {
    using TNode = LBChainedNode;
    std::vector<std::shared_ptr<TNode>> container;
    container.reserve((unsigned long) std::pow(2, 20));

    using PriorityQueue = std::multiset<std::shared_ptr<TNode>, CompareLBChainedNode>;
    PriorityQueue pQueue;
    /* Workload increase rate function, some examples are given below */
    int deltaW_func_id;
    /* Number of solution to get from Branch and Bound */
    int nb_solution_wanted;
    /* Initial workload */
    double W0;
    /* Workload increase load function */
    std::function<double(int)> deltaW;
    /* Load balancing cost */
    double C;
    /* Number of iteration to simulate */
    unsigned int maxI;
    /* Number of processors */
    unsigned int P;
    /* Application workload at each iteration */
    std::vector<double> W;
    std::string app_workload_file;

    cfg::ArgParser parser;
    parser.add_opt_version('V', "version", "0.1");
    parser.add_opt_help('h', "help"); // use -h or --help
    parser.add_opt_value('W', "W0", W0, (double) 0, "Initial workload", "DOUBLE").require(); // require this option

    parser.add_opt_value('d', "deltaW", deltaW_func_id, (int) 0,
                         "Select the workload increase rate function (11):"
                         " (0) dw=1.0,"
                         " (1) dw=5.0/(1.0+0.2*i,"
                         " (2) dw=1.0 + 0.2*i,"
                         " (3) dw=sin(0.2*i)", "INT"); // require this option
    parser.add_opt_value('p', "processor", P, (unsigned int) 0, "Number of processors",
                         "INT").require(); // require this option
    parser.add_opt_value('i', "iteration", maxI, maxI, "Number of iteration to simulate",
                         "INT").require(); // require this option
    parser.add_opt_value('I', "import", app_workload_file, std::string("app_workload_file"),
                         "import application workload file",
                         "STRING"); // require this option
    parser.add_opt_value('C', "lbcost", C, (double) 0, "Load balancing cost",
                         "DOUBLE").require(); // require this option
    parser.add_opt_value('s', "solution", nb_solution_wanted, (int) 1,
                         "Number of output solution from branch and bound", "INT"); // require this option

    bool output;
    auto &verbose = parser.add_opt_flag('v', "verbose", "(1) output iteration by iteration, (2) produce output ",
                                        &output);
    parser.parse(argc, argv);

    if (parser.count_error() > 0) {
        std::cout << parser.get_error() << std::endl;
        std::cout << parser.get_help() << std::endl;
        exit(-1);
    }

    std::random_device dev;
    auto seed = dev();
    std::cout << "The random seed is\t" << seed << std::endl;
    std::mt19937 rng(seed);

    std::uniform_real_distribution<double> uniform_dist(-0.3, 0.5);
    std::uniform_real_distribution<double> uniform_dist_derivative(0, W0 / 20);
    std::normal_distribution<double> normal_dist_derivative(2.0, 10.0);
    std::uniform_real_distribution<double> random_staircase(0.1, 1.1);
    std::normal_distribution<double> normal_dist(0.3, 1.3);

    std::vector<double> customLoad(maxI);
    std::vector<double> uniformLoad(maxI);
    std::vector<double> gaussianLoad(maxI);
    std::vector<double> sinLoad(maxI);
    std::vector<double> cosLoad(maxI);
    std::vector<double> nothing(maxI);
    std::vector<double> minus1(maxI);
    std::vector<double> randomPositiveDerivativeUniform(maxI);
    std::vector<double> randomPositiveDerivativeNormal(maxI);
    std::vector<double> exponential(maxI);

    for (int i = 49; i < maxI; i += 50) {
        customLoad[i] = random_staircase(rng);
    }

    for (int i = 0; i < maxI; i++) {
        uniformLoad[i] = uniform_dist(rng);
        randomPositiveDerivativeUniform[i] = uniform_dist_derivative(rng);
        randomPositiveDerivativeNormal[i] = normal_dist_derivative(rng);
        gaussianLoad[i] = normal_dist(rng);
        sinLoad[i] = 1.0 + std::sin(0.2 * i);
        cosLoad[i] = std::cos(i * 5.0 * M_PI / 180.0) - 0.1;
        nothing[i] = 0.0;
        minus1[i] = i < maxI / 2 ? -1.0 : 1.0;
        exponential[i] = std::exp(i / (maxI * 0.1));
    }

    constexpr int NB_INCREASING_WORKLOAD_F = 13;

    std::function<double(int)> deltaWf[NB_INCREASING_WORKLOAD_F] = {
            [W0, P](int i) { return (double) 0.5 * W0 / P; },
            [](int i) { return (double) 5.0 / (1.0 + 0.2 * i); },
            [](int i) { return (double) i * 0.0001; },
            [&sinLoad](int i) { return sinLoad[i]; },
            [&cosLoad](int i) { return cosLoad[i]; },
            [&uniformLoad](int i) { return uniformLoad[i]; },
            [&gaussianLoad](int i) { return gaussianLoad[i]; },
            [&customLoad](int i) { return (double) customLoad[i]; },
            [&nothing](int i) { return (double) nothing[i]; },
            [&minus1](int i) { return (double) minus1[i]; },
            [&randomPositiveDerivativeUniform](int i) { return (double) randomPositiveDerivativeUniform[i]; },
            [&randomPositiveDerivativeNormal, W0, P](int i) {
                return (double) std::max(1.0, randomPositiveDerivativeNormal[i]);
            },
            [&exponential, W0, P](int i) { return (double) exponential[i]; },
    };

    deltaW_func_id = deltaW_func_id > NB_INCREASING_WORKLOAD_F ? 0 : deltaW_func_id;
    deltaW = deltaWf[deltaW_func_id];

    W.resize(maxI);
    for (unsigned int i = 0; i < maxI; ++i) W[i] = std::max(0.0, _W(W0, i, deltaW));

    std::cout << W << std::endl;
    SimParam param{W0, W, C, maxI, P, deltaW};

    std::shared_ptr<TNode> initNode = std::make_shared<TNode>(&param);
    pQueue.insert(initNode);

    std::vector<std::shared_ptr<TNode>> solutions;
    std::vector<bool> foundYes(maxI, false);

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    /* Branch and bound search with node pruning */
    if (nb_solution_wanted)
        do {
            auto currentNode = *pQueue.begin();
            pQueue.erase(pQueue.begin());
            //Ok, I found a Yes Node for a given depth of the binary tree, no other Yes node at this depth can be better
            if (currentNode->apply_lb) {
                prune_similar_nodes_except(currentNode, pQueue, nb_solution_wanted - 1);
                foundYes[currentNode->iteration] = true;
            }

            if (currentNode->iteration >= maxI - 1) {
                solutions.push_back(currentNode);
            } else {
                auto[yes, no] = currentNode->children();
                //if I did not removed a Yes Node at this iteration, it might be better than what already exist
                if (!foundYes.at(no->iteration))
                    pQueue.insert(yes);
                pQueue.insert(no);
            }
        } while (solutions.size() != nb_solution_wanted);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    std::cout << std::setfill('=') << std::setw(130) << "\n";
    std::cout << "Results from search using branch and bound: (it took " << time_span.count() << " seconds)"
              << std::endl;

    enumerate(solutions.cbegin(), solutions.cend(), [](int i, std::shared_ptr<TNode> val) {
        std::cout << "BaB :" << std::string(13, ' ') << val << " "
                  << std::get<0>(compute_tcpu(get_scenario(val.get()), *val->params)) << std::endl;
    });
    std::cout << std::endl << "Results using mathematical formulation: " << std::endl;
    t1 = high_resolution_clock::now();

    std::cout << std::setfill('=') << std::setw(130) << "\n";

    auto[tmenon, sc1, imb] = create_scenario_menon1(param);

    auto menon_criterion_sol1 = generate_solution_from_scenario(sc1, param);

    std::cout << "U>C: " << std::string(13, ' ') << menon_criterion_sol1 << " " << tmenon << std::endl;

    /* Show the cumulative time (CPU_TIME) of a given solution up to a given iteration */

    if (verbose.get_count() >= 1) {
        std::cout << " Branch and Bound solutions: " << std::endl;
        show_each_iteration(solutions[0]);
        std::cout << " Menon solution  (U>C): " << std::endl;
        show_each_iteration(menon_criterion_sol1);
    }

    std::ofstream fMenon;
    auto menon_criterion_sol = generate_solution_from_scenario(sc1, param);
    fMenon.open("menon-solution.txt");
    fMenon << param.maxI << std::endl;
    fMenon << param.W << std::endl;
    fMenon << std::fixed << std::setprecision(6);
    fMenon << menon_criterion_sol->eval() << std::endl;
    fMenon << imb << std::endl;
    write_data(fMenon, menon_criterion_sol.get(), get_scenario, " ");
    write_data(fMenon, menon_criterion_sol.get(), get_times, " ");
    fMenon << param.C << std::endl;
    fMenon.close();

    {
        std::ofstream fMenon;
        auto[tmenon, sc1, imb] = create_scenario_menon_minus_one(param);
        auto menon_criterion_sol = generate_solution_from_scenario(sc1, param);
        fMenon.open("menon-solution-1.txt");
        fMenon << param.maxI << std::endl;
        fMenon << param.W << std::endl;
        fMenon << std::fixed << std::setprecision(6);
        fMenon << menon_criterion_sol->eval() << std::endl;
        fMenon << imb << std::endl;
        write_data(fMenon, menon_criterion_sol.get(), get_scenario, " ");
        write_data(fMenon, menon_criterion_sol.get(), get_times, " ");
        fMenon << param.C << std::endl;
        fMenon.close();
    }

    {
        std::ofstream fFreq100;
        auto[tfreq, scfreq, imb_freq] = create_scenario_freq(param, 100);
        auto freq_criterion_sol = generate_solution_from_scenario(scfreq, param);
        fFreq100.open("freq-100-solution.txt");
        fFreq100 << param.maxI << std::endl;
        fFreq100 << param.W << std::endl;
        fFreq100 << std::fixed << std::setprecision(6);
        fFreq100 << freq_criterion_sol->eval() << std::endl;
        fFreq100 << imb_freq << std::endl;
        write_data(fFreq100, freq_criterion_sol.get(), get_scenario, " ");
        write_data(fFreq100, freq_criterion_sol.get(), get_times, " ");
        fFreq100 << param.C << std::endl;
        fFreq100.close();
    }
    {
        std::ofstream fProca;
        auto[tproca, sc2, imb_proca] = create_scenario_procassini(param, 0.9, [](){return 0.96;});
        auto proca_criterion_sol = generate_solution_from_scenario(sc2, param);
        fProca.open("proca-solution.txt");
        fProca << param.maxI << std::endl;
        fProca << param.W << std::endl;
        fProca << std::fixed << std::setprecision(6);
        fProca << proca_criterion_sol->eval() << std::endl;
        fProca << imb_proca << std::endl;
        write_data(fProca, proca_criterion_sol.get(), get_scenario, " ");
        write_data(fProca, proca_criterion_sol.get(), get_times, " ");
        fProca << param.C << std::endl;
        fProca.close();
    }
    std::ofstream fOpti;
    for(int i = 0; i < nb_solution_wanted; ++i) {
        auto[t, imb] = compute_tcpu(get_scenario(solutions[i].get()), param);
        fOpti.open("optimal-solution-"+std::to_string(i)+".txt");
        fOpti << param.maxI << std::endl;
        fOpti << param.W << std::endl;
        fOpti << std::fixed << std::setprecision(6);
        fOpti << t << std::endl;
        fOpti << imb << std::endl;
        write_data(fOpti, solutions[i].get(), get_scenario," ");
        write_data(fOpti, solutions[i].get(), get_times," ");
        fOpti << param.C << std::endl;
        fOpti.close();
    }

    return 0;
}
