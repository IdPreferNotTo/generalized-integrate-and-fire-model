cmake_minimum_required(VERSION 3.10)
project(generalized-integrate-and-fire)

set(CMAKE_CXX_STANDARD 11)

find_package(Boost 1.49.0 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(generalized-integrate-and-fire neuron_models.cpp neuron_models.h linear_response.cpp linear_response.h main.cpp)

target_link_libraries(generalized-integrate-and-fire ${Boost_LIBRARIES})