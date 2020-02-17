cmake_minimum_required(VERSION 3.14)

project(LB_branch_and_bound)

set(CMAKE_CXX_STANDARD 17)
find_package (Threads)

add_executable(LB_branch_and_bound main.cpp zupply/src/zupply.cpp)

target_link_libraries (LB_branch_and_bound ${CMAKE_THREAD_LIBS_INIT})
target_include_directories(LB_branch_and_bound PUBLIC zupply/src)