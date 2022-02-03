cmake_minimum_required(VERSION 3.1)
project(PDE_solver VERSION 0.2.0 LANGUAGES CXX )
add_executable(pde_solver main.cpp Initiation.cpp Domain.cpp EllipticPDE.cpp Mesh.cpp Solver.cpp PostProcessing.cpp)

find_package(Eigen3 3 REQUIRED NO_MODULE)
target_include_directories(pde_solver PRIVATE /usr/include/python3.8 /usr/local/include/eigen)
target_link_libraries(pde_solver -lpython3.8 Eigen3::Eigen)

target_compile_options(pde_solver PUBLIC -fno-omit-frame-pointer)
set (CMAKE_CXX_FLAGS "-O0 -pg")