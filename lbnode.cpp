//
// Created by xetql on 17/02/2020.
//

#include <iostream>
#include "lbnode.hpp"

LBNode eval(LBNode&& n) {
    int maxI = n.apply_lb.size() - 1;
    LBNode tmp = n;
    while(tmp.iteration < maxI) {
        tmp = tmp.next();
    }
    return tmp;
}
LBNode eval(LBNode&  n, int until) {
    int maxI = std::min((int) n.apply_lb.size() - 1, until);
    LBNode tmp = n;
    while(tmp.iteration < maxI) {
        tmp = tmp.next();
    }
    return tmp;
}
LBNode eval(LBNode&& n, int until) {
    int maxI = std::min((int) n.apply_lb.size() - 1, until);
    LBNode tmp = n;
    while(tmp.iteration < maxI) {
        tmp = tmp.next();
    }
    return tmp;
}
LBChainedNode get_node_at(LBChainedNode&& head, int at){
    int maxI = std::min((int) head.iteration, at);
    maxI = maxI > 0 ? maxI : 0;
    if(head.iteration == at){
        return std::forward<LBChainedNode>(head);
    }

    LBChainedNode tmp = head;
    while(tmp.iteration != maxI){
        tmp = *tmp.pnode;
    }
    return tmp;
}
std::shared_ptr<LBChainedNode> eval(std::vector<bool> scenario, SimParam p) {
    std::shared_ptr<LBChainedNode> node = std::make_shared<LBChainedNode>(0, 0, 0, 0, scenario.at(0), nullptr, &p);
    auto size = scenario.size();
    while(node->iteration+1 < size) {
        node = node->next(scenario);
    }
    return node;
}

// reverse a linked list
void reverse(std::shared_ptr<LBChainedNode>& node){
    std::shared_ptr<LBChainedNode>
            prev = nullptr,
            curr = node,
            next = node->pnode;
    while(curr != nullptr) {
        curr->pnode = prev;
        prev = curr;
        curr = next;
        if(next != nullptr) next = next->pnode;
    }
    node = prev;
}

/* Show the cumulative time (CPU_TIME) of a given solution until a given iteration */
void show_each_iteration(std::shared_ptr<const LBChainedNode> n, int until) {
    auto head = n;
    for(int i = 0; i < until-1; ++i) {
        std::cout << "iteration " << std::right << std::setfill(' ') << std::setw(3) << (i+1) << " -> { ";
        for(int j = 0; j < i+1; ++j){
            std::cout << n->apply_lb << " ";
            if(n->pnode == nullptr) break;
            n = n->pnode;

        }
        std::cout << "} = " << n->cpu_time << std::endl;
        n = head;
    }
    std::cout << "iteration " << std::right << std::setfill(' ') << std::setw(3) << (until) << " -> { ";
    for(int j = 0; j < until; ++j){
        std::cout << n->apply_lb << " ";
        if(n->pnode == nullptr) break;
        n = n->pnode;
    }
    std::cout << "} = " << n->eval() << std::endl;
}

void show_each_iteration(LBNode& n, int until) {
    for(int i = 0; i < until; ++i)
        std::cout << eval(LBNode{0,0,{0}, n.apply_lb, n.params}, i) << std::endl;
}

std::vector<int> get_lb_iterations(const std::shared_ptr<LBChainedNode> n) {
    std::vector<int> ret;
    auto tmp = n;
    while(tmp != nullptr){
        if(tmp->apply_lb)
            ret.push_back(tmp->iteration);
        tmp = tmp->pnode;
    }
    return ret;
}