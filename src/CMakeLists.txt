cmake_minimum_required(VERSION 3.1)
project(PDE_solver VERSION 0.1.0 LANGUAGES CXX )
add_executable(pde_solver main.cpp numerical.cpp utils.cpp inputoutput.cpp)

target_include_directories(pde_solver PRIVATE /usr/include/python3.8)
target_link_libraries(pde_solver -lpython3.8)