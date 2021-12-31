#include <iostream>

#include "./headers/Utilities.h"
#include "./headers/Initiation.h"
#include "./headers/EllipticPDE.h"
#include "./headers/Domain.h"
#include "./headers/Mesh.h"
#include "./headers/Solver.h"
#include "./headers/PostProcessing.h"

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
// Why make this a pointer and why unique?
// Maybe using the "i" char is not the best alternative? It's only used for integers
    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();
    
// Changed the name to message to make it clearer
    vector<string> message = {"* Please enter desired coordinate system (Cartesian = 1, Polar = 2): ",
                        "* Run test case? (1 or 0): ",
                        "* Would you like to solve for the homogenous case (1 or 0): ",
                        "* Please enter desired solution method (LU = 0, Gauss-Seidel = 1): "};
    
    bool is_test_case = static_cast<bool>(utils->requestInput('i', 0, 1, message[1]));
    
    CoordinateSystem coordinate_sys;
    if(is_test_case) coordinate_sys = CoordinateSystem::Cartesian;
    else coordinate_sys = static_cast<CoordinateSystem>(utils->requestInput('i', 1, 2, message[0]));
    
    bool is_homog;
    if (!is_test_case) is_homog = static_cast<bool>(utils->requestInput('i', 0, 1, message[2]));
    else is_homog = true;

    std::unique_ptr<Initiation> Heat_2D = std::make_unique<EllipticPDE>(coordinate_sys, is_test_case, is_homog);
    
    std::unique_ptr<Domain> domain = std::make_unique<Domain>(Heat_2D);
    
    if(!is_homog) Heat_2D->setInHomogeneous(domain);

    Heat_2D->setBCs();
    std::unique_ptr<Mesh> mesh = std::make_unique<Mesh>(Heat_2D, domain);

    utils->print_matrix(mesh->getMesh());

    // TOD: Solve
    bool solution_method = static_cast<bool>(utils->requestInput('i', 0, 1, message[3]));

    // if (!solution_method) {
    //     std::unique_ptr<Solver> solver = std::make_unique<LU>(is_test_case);
    // }
    // else {
    //     std::unique_ptr<Solver> solver = std::make_unique<Seidel>(is_test_case);
    // }
    
    std::unique_ptr<Solver> solver = std::make_unique<LU>(is_test_case);

    solver->solve(Heat_2D, mesh);
    solver->setError(Heat_2D, mesh, domain);
    // std::cout << "printing solution \n" << std::endl;
    // utils->print_matrix(solver->getSolution());

    // std::cout << "printing error \n" << std::endl;
    // utils->print_vector(solver->getError());
    
    // TODO: post-processing
    std::unique_ptr<PostProcessing> postproc = std::make_unique<PostProcessing>(Heat_2D,mesh,solver);
    postproc->printError();
    postproc->exportResult();
    postproc->plotResult();
}