#include <iostream>
#include <queue>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>

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

    double eval() const {
        return latest(apply_lb) ?
               latest(cpu_time) + (W[iteration] / P) + C:
               latest(cpu_time) + (W[prev_lb]   / P) + deltaW(iteration) * (iteration-prev_lb);
    }

    double eval(int i) const {
        return get(apply_lb, i) ?
               get(cpu_time, i) + (W[iteration] / P) + C:
               get(cpu_time, i) + (W[prev_lb]   / P) + deltaW(iteration) * (iteration-prev_lb);
    }

    LBNode next(){
        return {iteration+1, get(apply_lb, iteration) ? iteration : prev_lb, append(cpu_time, eval(iteration)), apply_lb};
    }

    inline std::pair<LBNode, LBNode> children() {
        return {
            {iteration+1, latest(apply_lb) ? iteration : prev_lb, append(cpu_time, eval()), append(apply_lb, true)},
            {iteration+1, latest(apply_lb) ? iteration : prev_lb, append(cpu_time, eval()), append(apply_lb, false)}
        };
    }

    friend std::ostream &operator<<(std::ostream &os, const LBNode &node) {
        os << "iteration " << node.iteration << " -> { ";
        std::for_each(node.apply_lb.cbegin(), node.apply_lb.cend(), [&](auto val){os<< val << " ";});
        os << "} = " << node.eval();
        return os;
    }

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

int main(int argc, char** argv) {
    using PriorityQueue = std::priority_queue<LBNode, std::vector<LBNode>, CompareLBNode>;
    PriorityQueue pQueue;

    /* Initial workload */
    W0     = 10;

    /* Number of processors */
    P      = 10;

    /* Workload increase rate function, some examples are given below */
    deltaW = [](int i){return 1.0;};    // constant -> CPU_time is quadratic
    // deltaW = [](int i){return i;};   // linear
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
        if(currentNode.iteration >= maxI) {
            solutions.push_back(currentNode);
        } else {
            auto [ l, r ] = currentNode.children();
            pQueue.push(l);
            pQueue.push(r);
        }
    } while(solutions.size() < nb_solution_wanted);

    std::for_each(solutions.cbegin(), solutions.cend(), [&](auto val){std::cout << val << std::endl;});

    auto average_increase_load = (W[maxI-1] - W[0]) / maxI;

    std::vector<bool> apply_lb1(maxI, false);
    auto fLB1 = [&](int cnt){return std::sqrt(2*C / (average_increase_load*(1-1.0/P)));};
    for(int i = 1, lb = (int) fLB1(i); lb <= maxI; ++i, lb = (int) fLB1(i)) apply_lb1[lb] = true;
    LBNode n1{0, 0, {0}, apply_lb1};
    std::cout << n1 << std::endl;

    std::vector<bool> apply_lb2(maxI, false);
    auto fLB2 = [&](int cnt){return std::sqrt(2*C / (average_increase_load*(1+1.0/P)));};
    for(int i = 1, lb = (int) fLB2(i); lb <= maxI; ++i, lb = (int) fLB2(i)) apply_lb2[lb] = true;

    std::vector<bool> apply_lb3(maxI, false);
    auto fLB3 = [&](int cnt){return std::sqrt(2*C / (average_increase_load));};
    for(int i = 1, lb = (int) fLB3(i); lb <= maxI; ++i, lb = (int) fLB3(i)) apply_lb3[lb] = true;


    std::ofstream fCpuTime;
    for(int i = 0; i < nb_solution_wanted; ++i){
        fCpuTime.open("optimal-cpu-time-"+std::to_string(i)+".txt");
        fCpuTime << solutions[i].get_cpu_time();
        fCpuTime.close();
    }

    //

    return 0;
}