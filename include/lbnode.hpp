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

struct LBChainedNode;
std::shared_ptr<LBChainedNode> make_next(LBChainedNode* parent, bool nextDecision);

std::vector<bool>   get_scenario(const LBChainedNode* n);
std::vector<double> get_times(const LBChainedNode* n);
std::vector<double> get_imbalance_time(const LBChainedNode* n);
std::vector<double> get_cumulative_imbalance_time(const LBChainedNode* n);
std::vector<double> get_cumulative_time(const LBChainedNode* n);
void write_cpu_times(std::ostream&, LBChainedNode*, const char*);

/*Use LBChainedNode as it consumes less memory*/
struct LBChainedNode : std::enable_shared_from_this<LBChainedNode> {
    const unsigned int iteration, prev_lb;
    const double cpu_time;
    const bool   apply_lb = false;

    std::shared_ptr<LBChainedNode> pnode;
    const SimParam * params;
    const Application app;
    LBChainedNode(unsigned int iteration, unsigned int prevLb, double cpuTime, Application app, bool applyLb,
                  std::shared_ptr<LBChainedNode> pnode, const SimParam * params) :
                  iteration(iteration),
                  prev_lb(prevLb), cpu_time(cpuTime), app(std::move(app)), apply_lb(applyLb), pnode(std::move(pnode)), params(params) {
    }
public:
    explicit LBChainedNode(const SimParam* params):
        iteration(0), prev_lb(0), cpu_time(0),
        app(Application{params->P, params->W, params->W0 / params->P,params->W0 / params->P, 1.}),
        apply_lb(false), params(params) {}

    double predict() const {
        return eval() + params->h.at(iteration);
    }

    double eval()  const {
        double w;
        w = cpu_time + app.max;
        if (apply_lb)
            w += params->C;
        return w;
    }

    /* Get the possible children that may appear after the current scenario (apply_lb) */
    inline std::pair<std::shared_ptr<LBChainedNode>, std::shared_ptr<LBChainedNode>> children() {
        return {make_next(this, true),make_next(this, false) };
    }

    bool get_decision() const { return apply_lb; }

    friend std::shared_ptr<LBChainedNode> make_next(LBChainedNode* parent, bool nextDecision);
    friend std::vector<bool> get_scenario(const LBChainedNode* n);
    friend std::vector<double>  get_times(const LBChainedNode* n);
    friend std::vector<double> get_imbalance_time(const LBChainedNode* n);
    friend std::vector<double> get_cumulative_imbalance_time(const LBChainedNode* n);
    friend std::vector<double> get_cumulative_time(const LBChainedNode* n);
    friend void write_cpu_times(std::ostream&, LBChainedNode*, const char*);
};

template<class T, class GetDataF>
std::vector<T> get_value(const LBChainedNode* n, GetDataF&& f) {
    std::vector<T> scenario;
    scenario.reserve(n->iteration+1);
    const LBChainedNode* curr = n;
    while(curr != nullptr) {
        scenario.push_back(f(curr));
        curr = curr->pnode.get();
    }
    std::reverse(scenario.begin(), scenario.end());
    return scenario;
}

std::ostream &operator<<(std::ostream &os, const LBChainedNode &node);

std::shared_ptr<LBChainedNode> generate_solution_from_scenario(const std::vector<bool>& scenario, SimParam p);

void reverse(std::shared_ptr<LBChainedNode>& node);

/* Show the cumulative time (CPU_TIME) of a given solution until a given iteration */
void show_each_iteration(std::shared_ptr<LBChainedNode> n);

std::vector<int> get_lb_iterations(LBChainedNode* parent);

struct CompareLBChainedNode {
    inline bool operator()(const std::shared_ptr<LBChainedNode>& a, const std::shared_ptr<LBChainedNode>& b) const {
        return a->predict() < b->predict();
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

template<class Container>
void prune_similar_nodes_except(const std::shared_ptr<LBChainedNode>& n, Container& c, int except){
    auto it = c.begin();
    const auto end = c.end();
    while(it != end) {
        auto current = it++; // copy the current iterator then increment it
        auto node = *current;
        if(node->iteration == n->iteration && node->apply_lb) {
            if(!except)
                c.erase(current);
            else
                except--;
        }
    }
}

template<class DataGetter, class NodeType>
void write_data(std::ostream& str, NodeType* node, DataGetter getField, const char* separator = "\n"){
    auto data = getField(node);
    std::for_each(data.cbegin(), data.cend(), [&](auto val){ str << val << separator;});
    str << std::endl;
}

#endif //LB_BRANCH_AND_BOUND_LBNODE_HPP
