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
#include <workload.hpp>

#include "lbnode.hpp"
#include "utils.hpp"
#include "io.hpp"
#include "scenario.hpp"

using namespace zz;


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
    workload::WorkloadIncreaseRate deltaW;
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
    std::normal_distribution<double> normal_dist(1.0, 3.0);

    std::vector<double> customLoad(maxI),
                        uniformLoad(maxI),
                        gaussianLoad(maxI),
                        sinLoad(maxI),
                        cosLoad(maxI),
                        nothing(maxI),
                        minus1(maxI),
                        randomPositiveDerivativeUniform(maxI),
                        randomPositiveDerivativeNormal(maxI),
                        exponential(maxI),
                        min_exponential(maxI);

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
        min_exponential[i] = 1 - std::exp(i / (maxI * 1.0));
    }

    constexpr int NB_INCREASING_WORKLOAD_F = 15;

    workload::WorkloadIncreaseRate deltaWf[NB_INCREASING_WORKLOAD_F] = {
            workload::Constant  {0.2* W0 / P},
            workload::Sublinear { 0.1, 0.5, 10.},
            workload::Linear    {-2., 18.},
            workload::Quadratic {-1, 18, 0},
            workload::SymmetricLinear {(int) maxI / 2, 0, 0.2*W0/P},
    };

//    workload::WorkloadIncreaseRate deltaWf[6] = {
//            workload::Constant{0.01 * W0 / P},
//            workload::Constant{-0.01 * W0 / P},
//            workload::Sublinear{5.0, 0.2, 1.0},
//            workload::Uniform(maxI, -0.5, 0.5),
//            workload::Normal{maxI, 0.0, 0.5},
//            workload::Perturbation{
//                workload::Sine{0.012},
//                workload::Uniform{maxI, 10, -10},
//                workload::Normal(maxI, 0.0, 0.0)
//            }
//    };

    deltaW = deltaWf[deltaW_func_id];

    W.resize(maxI);
    for (unsigned int i = 0; i < maxI; ++i) {
        W[i] = std::max(0.0, compute_application_workload(W0, i, deltaW));
    }

    if(verbose.get_count()) std::cout << W << std::endl;

    SimParam param {W0, W, C, maxI, P, deltaW};

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
    std::cout << "Results from search using branch and bound: (it took " << time_span.count() << " seconds)" << std::endl;

    enumerate(solutions.cbegin(), solutions.cend(), [](int i, std::shared_ptr<TNode> val) {
        std::cout << "BaB :" << std::string(13, ' ') << val << " "
                  << std::get<0>(compute_tcpu(get_scenario(val.get()), *val->params, "./results/s1.txt")) << std::endl;
    });

    std::cout << std::endl << "Results using mathematical formulation: " << std::endl;
    t1 = high_resolution_clock::now();

    std::cout << std::setfill('=') << std::setw(130) << "\n";

    auto[tmenon, sc1, imb] = create_scenario_menon1(param);
    auto menon_criterion_sol1 = generate_solution_from_scenario(sc1, param);

    auto[a, sc, c] = create_scenario_bastien(param);

    std::cout << "---- " << std::endl;
    std::cout << param << std::endl;
    std::cout << menon_criterion_sol1->eval() << std::endl;
    std::cout << generate_solution_from_scenario(sc, param)->eval() << std::endl;
    std::cout << solutions[0]->eval() << std::endl;
    std::cout << "---- " << std::endl;
    std::cout << "U>C: "   << std::string(13, ' ') << menon_criterion_sol1->eval() << " " << tmenon << std::endl;
    std::cout << "U+1>C: " << std::string(11, ' ') << generate_solution_from_scenario(sc, param)->eval() << " " << a << std::endl;

    /* Show the cumulative time (CPU_TIME) of a given solution up to a given iteration */

    if (verbose.get_count() >= 1) {
        std::cout << " Branch and Bound solutions: " << std::endl;
        show_each_iteration(solutions[0]);
        std::cout << " Menon solution  (U>C)adsd: " << std::endl;
        show_each_iteration(generate_solution_from_scenario(sc, param));
    }

    //create_scenario_freq(param, 1000);
    create_scenario_freq(param, 100);
    create_scenario_freq(param, 50);
    create_scenario_freq(param, 25);

    create_scenario_static(param);

    create_scenario_procassini(param, 0.9, [](){return 1.0;});

    for(int i = 0; i < nb_solution_wanted; ++i) {
        std::cout << std::get<0>(compute_tcpu(get_scenario(solutions[i].get()), param, "results/optimal-solution-"+std::to_string(i)+".txt")) << std::endl;
        std::cout << solutions[i]->eval() << std::endl;
    }

    return 0;
}
