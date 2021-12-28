#include <iostream>

#include "headers/Utilities.h"
#include "headers/Initiation.h"
#include "headers/EllipticPDE.h"
#include "headers/Domain.h"
#include "headers/Mesh.h"

void initiateProgram () {
    // Print  welcome message and problem description

    cout << "PDE Solver started!\n";
    cout << "The current version of the solver is capable of solving the following 2D heat transfer problem. For a rectangular domain with customized ";
    cout << "width and height in the cartesian coordinate system (i.e. x and y), the steady-state temperature values for the internal nodes are evaluated. ";
    cout << "Only Dirichlet boundary conditions are allowed for homogeneous materials (only a single material for the whole domain is expected). ";
    cout << "PS:  the bottom-left corner of the domain conincides with the origion of the coordinate system.\n\n";
    cout << "Now, let's start with defining the problem values!\n\n";
}


int main () {
    initiateProgram();

    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();
    
    vector<string> m = {"* Please enter desired coordinate system (Cartesian = 1, Polar = 2): ",
                        "* Run test case? (1 or 0): ",
                        "* Would you like to solve for the homogenous case (1 or 0): "};
    
    bool is_test_case = static_cast<bool>(utils->requestInput('i', 0, 1, m[1]));
    
    CoordinateSystem coordinate_sys;
    if(is_test_case) coordinate_sys = CoordinateSystem::Cartesian;
    else coordinate_sys = static_cast<CoordinateSystem>(utils->requestInput('i', 1, 2, m[0]));
    
    bool is_homog;
    if (!is_test_case) is_homog = static_cast<bool>(utils->requestInput('i', 0, 1, m[2]));
    else is_homog = true;

    std::unique_ptr<Initiation> Heat_2D = std::make_unique<EllipticPDE>(coordinate_sys, is_test_case, is_homog);
    
    std::unique_ptr<Domain> domain = std::make_unique<Domain>(Heat_2D);
    
    if(!is_homog) Heat_2D->setInHomogeneous(domain);

    Heat_2D->setBCs();
    std::unique_ptr<Mesh> mesh = std::make_unique<Mesh>(Heat_2D, domain);

    utils->print_matrix(mesh->getMesh());

    // TOD: Solve
    // TODO: post-processing
}