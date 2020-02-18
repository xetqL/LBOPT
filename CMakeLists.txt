cmake_minimum_required(VERSION 3.14)

project(LB_branch_and_bound)

set(CMAKE_CXX_STANDARD 17)
find_package (Threads)

add_executable(lb_exhaustive_search main.cpp zupply/src/zupply.cpp lbnode.cpp lbnode.hpp simparam.hpp utils.cpp utils.hpp)

target_link_libraries (lb_exhaustive_search ${CMAKE_THREAD_LIBS_INIT})

target_include_directories(lb_exhaustive_search PUBLIC zupply/src)