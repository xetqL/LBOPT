//
// Created by xetql on 17/02/2020.
//

#ifndef LB_BRANCH_AND_BOUND_LBNODE_HPP
#define LB_BRANCH_AND_BOUND_LBNODE_HPP

#include <memory>
#include <vector>
#include "simparam.hpp"
#include "utils.hpp"
#include <sstream>
#include <iomanip>
#include <cassert>

/*Use LBChainedNode as it consumes less memory*/
struct LBChainedNode : std::enable_shared_from_this<LBChainedNode> {

    unsigned int iteration = 0;
    unsigned int prev_lb   = 0;
    double cpu_time = 0;
    mutable double Wmax     = 0;
    bool   apply_lb = false;
    std::shared_ptr<LBChainedNode> pnode;
    const SimParam * params;

    LBChainedNode(unsigned int iteration, unsigned int prevLb, double cpuTime, double Wmax, bool applyLb,
                  std::shared_ptr<LBChainedNode> pnode, const SimParam * params) : iteration(iteration), prev_lb(prevLb),
                                                          cpu_time(cpuTime), Wmax(Wmax), apply_lb(applyLb), pnode(std::move(pnode)), params(params)
    {

    }

    /* Functions to evaluate the next computing time given the current state */
    inline double eval() const { return eval(iteration); }
    inline double eval2() const { return eval(iteration); }

    inline double eval(unsigned int iteration) const {
        return apply_lb ?
               cpu_time + (params->W.at(iteration) / params->P) + params->C :
               cpu_time + (params->W.at(prev_lb)   / params->P) + sum(prev_lb, iteration, params->deltaW);
    }
    inline double eval2(unsigned int iteration) const {
        if(apply_lb) {
            Wmax = (params->W.at(iteration) / params->P);
            return cpu_time + Wmax + params->C;
        } else {
            Wmax += params->deltaW(iteration);
            return cpu_time + Wmax;
        }
    }

    std::shared_ptr<LBChainedNode> next(const std::vector<bool>& apply_lb) {
        return std::make_shared<LBChainedNode>(iteration + 1,
                                               get(apply_lb, iteration) ? iteration : prev_lb,
                                               eval(iteration),
                                               Wmax,
                                               get(apply_lb, iteration + 1),
                                               this->shared_from_this(),
                                               params);
    }

    /* Get the possible children that may appear after the current scenario (apply_lb) */
    inline std::pair<std::shared_ptr<LBChainedNode>, std::shared_ptr<LBChainedNode>> children() {
        return
                {
                    std::make_shared<LBChainedNode>(
                            iteration+1,
                            apply_lb ? iteration : prev_lb,
                            eval2(),
                            Wmax,
                            true,
                            this->shared_from_this(),
                            params),
                    std::make_shared<LBChainedNode>(
                            iteration+1,
                            apply_lb ? iteration : prev_lb,
                            eval2(),
                            Wmax,
                            false,
                            this->shared_from_this(),
                            params)
                };
    }

    std::vector<bool> get_scenario(std::vector<bool>& scenario, const LBChainedNode* n) const {
        if(n == nullptr) {
            return scenario;
        } else {
            scenario.at(n->iteration) = (n->apply_lb);
            return get_scenario(scenario, n->pnode.get());
        }
    }

    std::vector<double> get_times(std::vector<double>& scenario, const LBChainedNode* n) const {
        scenario.at(n->iteration) = n->eval();
        if(n->pnode == nullptr) {
            return scenario;
        } else {
            return get_times(scenario, n->pnode.get());
        }
    }

    friend std::ostream &operator<<(std::ostream &os, const LBChainedNode &node) {
        std::vector<bool> scenario(node.iteration+1);
        scenario = node.get_scenario(scenario, &node);
        os << "iteration " << std::right << std::setfill(' ') << std::setw(3) << (node.iteration+1) << " -> { ";
        std::for_each(scenario.begin(), scenario.begin() + node.iteration + 1, [&](auto val){ os << val << " ";});
        os << "} = " << node.eval2();
        return os;
    }

    /* The node must be evaluated once */
    std::string get_cpu_time() {
        std::ostringstream str;
        std::vector<double> scenario(iteration+1);
        scenario = get_times(scenario, this);
        std::for_each(scenario.cbegin(), scenario.cend(), [&](auto val){str << std::to_string(val) << std::endl;});
        //str << std::to_string(eval());
        return str.str();
    }

};

/* Structure of LBNode, it contains the full CPU_Time history and the LB history*/
struct LBNode {
    unsigned int iteration = 0;
    unsigned int prev_lb   = 0;
    std::vector<double> cpu_time;
    std::vector<bool>   apply_lb;
    const SimParam *params;
    LBNode(unsigned int iteration, unsigned int prevLb, std::vector<double> cpuTime, std::vector<bool> applyLb,
            const SimParam * params):
        iteration(iteration),
        prev_lb(prevLb),
        cpu_time(std::move(cpuTime)),
        apply_lb(std::move(applyLb)),
        params(params) {}
    /* Functions to evaluate the next computing time given the current scenario */
    inline double eval() const { return eval(iteration); }
    inline double eval(int i) const {
        auto v = get(apply_lb, i) ?
                 get(cpu_time, i) + (params->W[iteration] / params->P) + params->C:
                 get(cpu_time, i) + (params->W[prev_lb]   / params->P) + sum(prev_lb, iteration, params->deltaW);
        return v;
    }

    /* Get the next node given the current scenario (apply_lb) */
    LBNode next() const {
        return {iteration+1,
                get(apply_lb, iteration) ? iteration : prev_lb,
                append(cpu_time, eval(iteration)),
                apply_lb,
                params};
    }
    /* Get the possible children that may appear after the current scenario (apply_lb) */
    inline std::pair<LBNode, LBNode> children() {
        return
                {
                        {iteration+1, latest(apply_lb) ? iteration : prev_lb, append(cpu_time, eval()), append((apply_lb), true), params},
                        {iteration+1, latest(apply_lb) ? iteration : prev_lb, _append(std::move(cpu_time), eval()), _append(std::move(apply_lb), false), params}
                };
    }

    friend std::ostream &operator<<(std::ostream &os, const LBNode &node) {
        os << "iteration " << std::right << std::setfill(' ') << std::setw(3) << (node.iteration+1) << " -> { ";
        std::for_each(node.apply_lb.begin(), node.apply_lb.begin() + node.iteration + 1, [&](auto val){ os << val << " ";});
        os << "} = " << node.eval();
        return os;
    }

    /* The node must be evaluated once */
    std::string get_cpu_time() {
        std::ostringstream str;
        std::for_each(cpu_time.cbegin(), cpu_time.cend(), [&](auto val){str << std::to_string(val) << std::endl;});
        str << std::to_string(eval());
        return str.str();
    }

    bool operator<(const LBNode &rhs)  const {
        return this->eval() < rhs.eval();
    }
    bool operator>(const LBNode &rhs)  const {
        return rhs.eval() < this->eval();
    }
    bool operator<=(const LBNode &rhs) const {
        return !(rhs.eval() < this->eval());
    }
    bool operator>=(const LBNode &rhs) const {
        return !(this->eval() < rhs.eval());

    }
};

LBNode eval(LBNode&& n);
LBNode eval(LBNode&  n, int until);
LBNode eval(LBNode&& n, int until);
LBChainedNode get_node_at(LBChainedNode&& head, int at);
std::shared_ptr<LBChainedNode> eval(std::vector<bool> scenario);
void reverse(std::shared_ptr<LBChainedNode>& node);

/* Show the cumulative time (CPU_TIME) of a given solution until a given iteration */
void show_each_iteration(std::shared_ptr<const LBChainedNode> n, int until);
void show_each_iteration(LBNode& n, int until);

std::vector<int> get_lb_iterations(std::shared_ptr<LBChainedNode> n);

/* Structure for comparing LBNode */
struct CompareLBNode {
    bool operator()(const LBNode &a, const LBNode &b){
        return a.eval() > b.eval();
    }
};

struct CompareLBChainedNode {
    inline bool operator()(const std::shared_ptr<LBChainedNode>& a, const std::shared_ptr<LBChainedNode>& b) const {
        return a->eval2() < b->eval2();
    }
};

template<class Container>
void prune_similar_nodes(const std::shared_ptr<LBChainedNode>& n, Container& c){
    auto it = c.begin();
    const auto end = c.end();
    while(it != end) {
        auto current = it++; // copy the current iterator then increment it
        auto node = *current;
        if(node->iteration == n->iteration && node->apply_lb) {
            c.erase(current);
        }
    }
}

#endif //LB_BRANCH_AND_BOUND_LBNODE_HPP
