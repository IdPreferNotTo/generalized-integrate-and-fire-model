cmake_minimum_required(VERSION 3.13)
project(SCH13)

set(CMAKE_CXX_STANDARD 14)

find_package(Boost 1.49.0 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(SCH13 main.cpp neuron_models.cpp neuron_models.h linear_response.cpp linear_response.h)

target_link_libraries(SCH13 ${Boost_LIBRARIES})