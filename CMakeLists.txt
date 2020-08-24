cmake_minimum_required(VERSION 3.12)

project(LB_branch_and_bound)

#set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 17)

find_package (Threads)

add_executable(lb main.cpp zupply/src/zupply.cpp lbnode.cpp lbnode.hpp simparam.hpp utils.cpp utils.hpp io.cpp io.hpp scenario.cpp scenario.hpp)

target_link_libraries (lb ${CMAKE_THREAD_LIBS_INIT})

target_include_directories(lb PUBLIC zupply/src)
