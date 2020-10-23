//
// Created by xetql on 17/02/2020.
//

#include <iostream>
#include "lbnode.hpp"
#include "utils.hpp"

std::shared_ptr<LBChainedNode> generate_solution_from_scenario(const std::vector<bool>& scenario, SimParam p) {
    std::shared_ptr<LBChainedNode> node = std::make_shared<LBChainedNode>(&p);
    for(auto it = (std::begin(scenario)+1); it != std::end(scenario); it++) {
        node = make_next(node.get(), *it);
    }
    return node;
}

void reverse(std::shared_ptr<LBChainedNode>& node){
    std::shared_ptr<LBChainedNode>
            prev = nullptr,
            curr = node,
            next = node->pnode;
    while(curr != nullptr) {
        curr->pnode = prev;
        prev = curr;
        curr = next;
        if(next != nullptr)
            next = next->pnode;
    }
    node = prev;
}

void show_each_iteration(std::shared_ptr<LBChainedNode> n) {
    reverse(n);
    for(int i = 0; i < n->params->maxI; ++i) {
        auto head = n;
        int readable_index = i+1;
        std::cout << "iteration " << std::right << std::setfill(' ') << std::setw(3) << (readable_index) << " -> { ";
        std::cout << head->get_decision() << " ";
        for(int j = 0; j < i; ++j) {
            head = head->pnode;
            if(head == nullptr) break;
            std::cout << head->get_decision() << " ";
        }
        std::cout << "} = " << head->eval() << std::endl;
    }
    reverse(n);
}



std::shared_ptr<LBChainedNode> make_next(LBChainedNode* parent, bool nextDecision){
    auto P = parent->params->P;

    auto nextIteration = parent->iteration+1;

    auto previous_lb_call = parent->apply_lb ? parent->iteration : parent->prev_lb;

    State s = parent->s;
    update_workloads(parent->iteration, P, parent->params->deltaW, s);

    if(nextDecision) {
        rebalance(s);
    }

    return std::make_shared<LBChainedNode>(nextIteration, previous_lb_call, parent->eval(), s, nextDecision, parent->shared_from_this(), parent->params);
}

std::vector<int>    get_lb_iterations(LBChainedNode* parent) {
    std::vector<int> ret;
    auto tmp = parent;
    while(tmp != nullptr) {
        if(tmp->apply_lb)
            ret.push_back(tmp->iteration);
        tmp = tmp->pnode.get();
    }
    std::reverse(ret.begin(), ret.end());
    return ret;
}
std::vector<bool>   get_scenario(const LBChainedNode* n) {
    std::vector<bool> scenario;
    const LBChainedNode* curr = n;
    while(curr != nullptr) {
        scenario.push_back(curr->get_decision());
        curr = curr->pnode.get();
    }
    std::reverse(scenario.begin(), scenario.end());
    return scenario;
}
std::vector<double> get_times(const LBChainedNode* n) {
    std::vector<double> cpu_times;
    auto curr = n;
    while(curr != nullptr) {
        cpu_times.push_back(curr->eval());
        curr = curr->pnode.get();
    }
    std::reverse(cpu_times.begin(), cpu_times.end());
    return cpu_times;
}
std::vector<double> get_imbalance_time(const LBChainedNode* n) {
    std::vector<double> cpu_times;
    auto curr = n;
    while(curr != nullptr){
        cpu_times.push_back(curr->s.max - curr->s.avg);
        curr = curr->pnode.get();
    }
    std::reverse(cpu_times.begin(), cpu_times.end());
    return cpu_times;
}
std::vector<double> get_cumulative_imbalance_time(const LBChainedNode* n){
    auto scenario = get_scenario(n);
    auto imb_time = get_imbalance_time(n);
    std::vector<double> cum_imb_time(imb_time.size());
    auto size = scenario.size();
    double U = 0;
    for(int i = 0; i < size; ++i) {
        if(scenario[i]) U = imb_time[i];
        else U += imb_time[i];
        cum_imb_time[i] = U;
    }
    return cum_imb_time;
}
std::vector<double> get_cumulative_time(const LBChainedNode* n){
    auto scenario = get_scenario(n);
    auto times = get_times(n);
    std::vector<double> cum_times(times.size());
    auto size = scenario.size();
    double U = 0;
    for(int i = 0; i < size; ++i) {
        if(scenario[i])
            U = 0;
        else
            U += times[i];
        cum_times[i] = U;
    }
    return cum_times;
}

/* The node must be evaluated once */
void write_cpu_times(std::ostream& str, LBChainedNode* node, const char* separator = "\n") {
    std::vector<double> times = get_times(node);
    std::for_each(times.cbegin(), times.cend(), [&](auto val){str << val << separator;});
    str << std::endl;
}

std::ostream &operator<<(std::ostream &os, const LBChainedNode &node) {
    //std::vector<bool> scenario = get_value<bool>(&node, [](const LBChainedNode* n){return n->get_decision();} );
    std::vector<bool> scenario = get_scenario(&node);
    os << "iteration " << std::right << std::setfill(' ') << std::setw(3) << (scenario.size()) << " -> { ";
    std::for_each(scenario.begin(), scenario.end(), [&](auto val){ os << val << " ";});
    os << "} = " << node.eval();
    return os;
}