//
// Created by xetql on 8/24/20.
//

#include "io.hpp"
#include "lbnode.hpp"

void save_results(std::string fname,
                  const std::vector<bool>& scenario,
                  const SimParam& p,
                  const std::vector<double>& imb_time,
                  const std::vector<double>& it_max,
                  const std::vector<double>& it_avg) {
    std::ofstream f;
    auto sol = generate_solution_from_scenario(scenario, p);
    f.open(fname);
    f << p.maxI << std::endl;
    f << p.W << std::endl;
    f << std::fixed << std::setprecision(6);
    f << sol->eval() << std::endl;
    f << imb_time   << std::endl;
    write_data(f, sol.get(), get_scenario, " ");
    write_data(f, sol.get(), get_times, " ");
    f << p.C << std::endl;
    f << it_max << std::endl;
    f << it_avg << std::endl;
    f.close();
}
