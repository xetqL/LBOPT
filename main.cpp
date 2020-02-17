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
#include "lbnode.hpp"

using namespace zz;



template<typename T>
std::ostream& operator<<(std::ostream& os, const std::shared_ptr<T>& pc) {
    os << *pc;
    return os;
}

int main(int argc, char** argv) {
    using TNode = LBChainedNode;
    std::vector<std::shared_ptr<TNode>> container;
    container.reserve(std::pow(2, 20));
    //using PriorityQueue = std::priority_queue<std::shared_ptr<TNode>, std::vector<std::shared_ptr<TNode>>, CompareLBChainedNode>;
    using PriorityQueue = std::multiset<std::shared_ptr<TNode>, CompareLBChainedNode>;
    PriorityQueue pQueue;
    /* Workload increase rate function, some examples are given below */
    int deltaW_func_id;
    std::function<double(int)> deltaWf[3] = {
            [](int i){return (double) 1.0;},
            [](int i){return (double) 5.0/(1.0+0.2*i);},
            [](int i){return std::sin((double) i);}
    };

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

    cfg::ArgParser parser;
    parser.add_opt_version('V', "version", "0.1");
    parser.add_opt_help('h', "help"); // use -h or --help
    parser.add_opt_value('W', "W0",W0, (double) 0, "Initial workload", "DOUBLE").require(); // require this option
    parser.add_opt_value('d', "deltaW", deltaW_func_id, (int) 0, "Select the workload increase rate function: default (0): dw=1, (1) dw=i, (2) dw=sin(i)", "INT"); // require this option
    parser.add_opt_value('p', "processor",  P, (unsigned int) 0, "Number of processors", "INT").require(); // require this option
    parser.add_opt_value('i', "iteration",  maxI, (unsigned int) maxI, "Number of iteration to simulate", "INT").require(); // require this option
    parser.add_opt_value('C', "lbcost",     C, (double) 0, "Load balancing cost", "DOUBLE").require(); // require this option
    parser.add_opt_value('s', "solution", nb_solution_wanted, (int) 1, "Number of output solution from branch and bound, default = 1", "INT"); // require this option

    bool output;
    parser.add_opt_flag('o', "out", "Produce output (yes/no)", &output);

    parser.parse(argc, argv);

    if (parser.count_error() > 0) {
        std::cout << parser.get_error() << std::endl;
        std::cout << parser.get_help() << std::endl;
        exit(-1);
    }

    deltaW_func_id = deltaW_func_id > 2 ? 0 : deltaW_func_id;

    deltaW = deltaWf[deltaW_func_id];

    W.resize(maxI+1);
    for(int i = 0; i < maxI; ++i) W[i] = _W(W0, i, deltaW);

    SimParam param {W0, W, C, maxI, P, deltaW};

    std::shared_ptr<TNode> initNode = std::make_shared<TNode>(0, 0, 0, false, nullptr, &param);
    pQueue.insert(initNode);

    std::vector<std::shared_ptr<TNode>> solutions;
    std::vector<bool> foundYes (maxI, false);
    do {
        std::shared_ptr<TNode> currentNode = *pQueue.begin();
        //std::cout << currentNode.iteration << " " << (currentNode.pnode != nullptr ? currentNode.pnode->iteration:0) << std::endl;
        pQueue.erase(pQueue.begin());

        //Ok, I found a Yes Node for a given depth of the binary tree, no other Yes node at this depth can be better atm
        //remove them
        if(currentNode->apply_lb){
            auto it = pQueue.begin();
            const auto end = pQueue.end();
            while(it != end) {
                // copy the current iterator then increment it
                auto current = it++;
                auto node = *current;
                if(node->iteration == currentNode->iteration && node->apply_lb){
                    pQueue.erase(current);
                }
            }
            foundYes[currentNode->iteration] = true;
        }

        if(currentNode->iteration >= maxI-1) {
            solutions.push_back(currentNode);
        } else {
            auto [ yes, no ] = currentNode->children();
            //if I did not removed a Yes Node at this iteration, it might be better than what already exist

            if(!foundYes.at(no->iteration))
                pQueue.insert(yes);
            pQueue.insert(no);
        }
        //std::cout << pQueue.size() << std::endl;
    } while(solutions.size() < nb_solution_wanted);

    std::cout << std::setfill('=') << std::setw(130) << "\n";
    std::cout << "Results from search using branch and bound: " << std::endl ;

    enumerate(solutions.cbegin(), solutions.cend(), [](int i, std::shared_ptr<TNode> val){
        std::cout << std::right << std::setw (22) << std::fixed << std::setprecision(5) << std::setfill(' ')
        << "Solution " + std::to_string(i) << ": " << std::string(13, ' ') << val << std::endl;
    });

    std::cout << std::endl << "Results using mathematical formulation: "<<std::endl;
    auto average_increase_load = (W[maxI-1] - W[0]) / (maxI-1);
    std::vector<bool> apply_lb1(maxI, false);
    auto fLB1 = [&](int cnt){return cnt*std::sqrt(2*C / (average_increase_load*(1-1.0/P)));};
    for(int i = 1, lb = std::round(fLB1(i)); lb <= maxI; ++i, lb = std::round(fLB1(i))) apply_lb1[lb-1] = true;
    std::cout << std::left << std::setw (24) << std::fixed << std::setprecision(5) << std::setfill(' ')
              << "sqrt[2*C/dW*(1-1.0/P)]: "<<"Tau=" << fLB1(1) << "; "
              << eval(LBNode{0, 0, {0}, apply_lb1, &param}) << std::endl;

    std::vector<bool> apply_lb2(maxI, false);
    auto fLB2 = [&](int cnt){return cnt*std::sqrt(2*C / (average_increase_load*(1+1.0/P)));};
    for(int i = 1, lb = std::round(fLB2(i)); lb <= maxI; ++i, lb = std::round(fLB2(i))) apply_lb2[lb-1] = true;
    std::cout << std::left << std::setw (24) << std::fixed << std::setprecision(5) << std::setfill(' ')
              << "sqrt[2*C/dW*(1+1.0/P)]: "<<"Tau=" << fLB2(1) << "; "
              << eval(LBNode{0, 0, {0}, apply_lb2, &param}) << std::endl;

    std::vector<bool> apply_lb3(maxI, false);
    auto fLB3 = [&](int cnt){return cnt*std::sqrt(2*C / (average_increase_load));};
    for(int i = 1, lb = std::round(fLB3(i)); lb <= maxI; ++i, lb = std::round(fLB3(i))) apply_lb3[lb-1] = true;
    std::cout << std::left << std::setw (24) << std::fixed << std::setprecision(5) << std::setfill(' ')
              <<  "sqrt(2*C/dW):  "<<"Tau=" << fLB3(1) <<"; "
              << eval(LBNode{0, 0, {0}, apply_lb3, &param}) << std::endl;

    std::cout << std::setfill('=') << std::setw(130) << "\n";

    std::vector<bool> s(maxI);
    std::vector<bool> apply_lb4 = solutions[0]->get_scenario(s, solutions[0].get());


    std::cout << std::left << std::setw (24) << std::fixed << std::setprecision(5) << std::setfill(' ')
              <<  "Verification:  "<<"Tau=" <<"; "
              << eval(LBNode{0, 0, {0}, apply_lb4, &param}) << std::endl;

    std::cout << std::setfill('=') << std::setw(130) << "\n";


    if(output){
        std::ofstream fCpuTime;
        for(int i = 0; i < nb_solution_wanted; ++i){
            fCpuTime.open("optimal-cpu-time-"+std::to_string(i)+".txt");
            fCpuTime << solutions[i]->get_cpu_time();
            fCpuTime.close();
        }
    }

    /* Show the cumulative time (CPU_TIME) of a given solution until a given iteration */
    show_each_iteration(solutions[0], maxI);

    return 0;
}