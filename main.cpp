#include <iostream>
#include <queue>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>
#include <functional>
#include <iomanip>

double W0;
std::function<double(int)> deltaW;
double C;
unsigned int maxI;
unsigned int P;

template<class FMath>
double sum(int b, int e, FMath& f){
    double r = 0;
    for(int i = b; i < e; ++i) r += f(i);
    return r;
}

/* Function to compute the application workload at a given iteration */
inline double _W(unsigned int i) {
    return W0 + sum(0, i, deltaW);
}
/* Application workload at each iteration */
std::vector<double> W;

/* Append value and create new vector */
template<class T>
inline std::vector<T> append(const std::vector<T>& prev, T v){
    std::vector<T> next(prev.begin(), prev.end());
    next.push_back(v);
    return next;
}

/* Get the latest element of the vector */
template<class T>
inline T latest(const std::vector<T>& prev) {
    return *(prev.end() - 1);
}
template<class T>
inline T get(const std::vector<T>& prev, int i){
    return *(prev.begin() + i);
}

/* Structure of LBNode, it contains the full CPU_Time history and the LB history*/
struct LBNode {
    unsigned int iteration = 0;
    unsigned int prev_lb   = 0;
    std::vector<double> cpu_time;
    std::vector<bool>   apply_lb;

    /* Functions to evaluate the next computing time given the current scenario */
    double eval() const { return eval(iteration); }
    double eval(int i) const {
        auto v = get(apply_lb, i) ?
                 get(cpu_time, i) + (W[iteration+1] / P) + C:
                 get(cpu_time, i) + (W[prev_lb == 0 ? 0 : prev_lb+1] / P) + deltaW(iteration) * (iteration - prev_lb);
        return v;
    }

    /* Get the next node given the current scenario (apply_lb) */
    LBNode next() const {
        return {iteration+1,
                get(apply_lb, iteration) ? iteration : prev_lb,
                append(cpu_time, eval(iteration)),
                apply_lb};
    }
    /* Get the possible children that may appear after the current scenario (apply_lb) */
    inline std::pair<LBNode, LBNode> children() {
        return
        {
            {iteration+1, latest(apply_lb) ? iteration : prev_lb, append(cpu_time, eval()), append(apply_lb, true)},
            {iteration+1, latest(apply_lb) ? iteration : prev_lb, append(cpu_time, eval()), append(apply_lb, false)}
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
};

/* Structure for comparing LBNode */
struct CompareLBNode {
    bool operator()(LBNode &a, LBNode &b){
        return a.eval() > b.eval();
    }
};

LBNode eval(LBNode&& n) {
    int maxI = n.apply_lb.size() - 1;
    LBNode tmp = n;
    while(tmp.iteration < maxI) {
        tmp = tmp.next();
    }
    return tmp;
}
LBNode eval(LBNode& n,  int until) {
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

/* Enumerate (python-like) function, i.e., apply a binary function on a pair (int, T) where int is the id and T a
 * data to be processed. This function may or may not change the data via non-const iterator.*/
template<class InputIt, class BinaryOp>
void enumerate(InputIt begin, InputIt end, BinaryOp&& fb){
    auto i = 0;
    while(begin != end){
        fb(i, *begin);
        ++begin;
        ++i;
    }
}

/* Show the cumulative time (CPU_TIME) of a given solution until a given iteration */
void show_each_iteration(LBNode& n, int until){
    for(int i = 0; i < until; ++i)
        std::cout << eval(LBNode{0,0,{0}, n.apply_lb}, i) << std::endl;
}

int main(int argc, char** argv) {
    using PriorityQueue = std::priority_queue<LBNode, std::vector<LBNode>, CompareLBNode>;
    PriorityQueue pQueue;

    /* Initial workload */
    W0     = 10;

    /* Number of processors */
    P      = 10;

    /* Workload increase rate function, some examples are given below */
    deltaW = [](int i){return 1.0;};    // constant -> CPU_time is quadratic
    //deltaW = [](int i){return i;};   // linear
    // deltaW = [](int i){return std::sin((double) i);};

    /* Load balancing cost */
    C      = 25;

    /* Number of iteration to simulate */
    maxI   = std::atoi(argv[1]);

    W.resize(maxI);
    for(int i = 0; i < maxI; ++i) W[i] = _W(i);
    int nb_solution_wanted = std::atoi(argv[2]);

    LBNode initNode {0, 0, {0}, {false}};
    pQueue.push(initNode);
    std::vector<LBNode> solutions;

    do {
        LBNode currentNode = pQueue.top();
        pQueue.pop();
        if(currentNode.iteration >= maxI - 1) {
            solutions.push_back(currentNode);
        } else {
            auto [ l, r ] = currentNode.children();
            pQueue.push(l);
            pQueue.push(r);
        }
    } while(solutions.size() < nb_solution_wanted);
    std::cout << std::setfill('=') << std::setw(130) << "\n";
    std::cout << "Results from search using branch and bound: " << std::endl ;
    enumerate(solutions.cbegin(), solutions.cend(), [&](int i, auto val){
        std::cout << std::right << std::setw (22) << std::fixed << std::setprecision(5) << std::setfill(' ')
        << "Solution " + std::to_string(i) << ": " << std::string(13, ' ') << val << std::endl;
    });
    std::cout << std::endl << "Results using mathematical formulation: "<<std::endl;
    auto average_increase_load = (W[maxI-1] - W[0]) / maxI;

    std::vector<bool> apply_lb1(maxI, false);
    auto fLB1 = [&](int cnt){return cnt*std::sqrt(2*C / (average_increase_load*(1-1.0/P)));};
    for(int i = 1, lb = std::round(fLB1(i)); lb <= maxI; ++i, lb = std::round(fLB1(i))) apply_lb1[lb-1] = true;
    std::cout << std::left << std::setw (24) << std::fixed << std::setprecision(5) << std::setfill(' ')
              << "sqrt[2*C/dW*(1-1.0/P)]: "<<"Tau=" << fLB1(1) << "; "
              << eval(LBNode{0, 0, {0}, apply_lb1}) << std::endl;

    std::vector<bool> apply_lb2(maxI, false);
    auto fLB2 = [&](int cnt){return cnt*std::sqrt(2*C / (average_increase_load*(1+1.0/P)));};
    for(int i = 1, lb = std::round(fLB2(i)); lb <= maxI; ++i, lb = std::round(fLB2(i))) apply_lb2[lb-1] = true;
    std::cout << std::left << std::setw (24) << std::fixed << std::setprecision(5) << std::setfill(' ')
              << "sqrt[2*C/dW*(1+1.0/P)]: "<<"Tau=" << fLB2(1) << "; "
              << eval(LBNode{0, 0, {0}, apply_lb2}) << std::endl;

    std::vector<bool> apply_lb3(maxI, false);
    auto fLB3 = [&](int cnt){return cnt*std::sqrt(2*C / (average_increase_load));};
    for(int i = 1, lb = std::round(fLB3(i)); lb <= maxI; ++i, lb = std::round(fLB3(i))) apply_lb3[lb-1] = true;
    std::cout << std::left << std::setw (24) << std::fixed << std::setprecision(5) << std::setfill(' ')
              <<  "sqrt(2*C/dW):  "<<"Tau=" << fLB3(1) <<"; "
              << eval(LBNode{0, 0, {0}, apply_lb3}) << std::endl;

    std::cout << std::setfill('=') << std::setw(130) << "\n";

    std::ofstream fCpuTime;
    for(int i = 0; i < nb_solution_wanted; ++i){
        fCpuTime.open("optimal-cpu-time-"+std::to_string(i)+".txt");
        fCpuTime << solutions[i].get_cpu_time();
        fCpuTime.close();
    }

    /* Show the cumulative time (CPU_TIME) of a given solution until a given iteration */
    show_each_iteration(solutions[0], maxI);

    return 0;
}