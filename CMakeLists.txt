cmake_minimum_required(VERSION 3.10)

project(LBOPT)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package (Threads)

file(GLOB SRC src/*.cpp)

add_executable(lb ${SRC} zupply/src/zupply.cpp )

target_link_libraries (lb ${CMAKE_THREAD_LIBS_INIT})

target_include_directories(lb PUBLIC zupply/src include)
