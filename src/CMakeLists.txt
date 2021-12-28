cmake_minimum_required(VERSION 3.1)
project(PDE_solver VERSION 0.2.0 LANGUAGES CXX )
add_executable(pde_solver main.cpp Initiation.cpp Domain.cpp EllipticPDE.cpp Mesh.cpp PostProcessing.cpp)

target_include_directories(pde_solver PRIVATE /usr/include/python3.8)
target_link_libraries(pde_solver -lpython3.8)
