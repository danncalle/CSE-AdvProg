#include <iostream>
#include <chrono>

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
    cout << "The current version (v2.0) of the solver is intented for the numerical solving of the 2D heat equation on numerous domains.\n";
    cout << "Due to the structure of the implementation, the capabilities and features of the solver can be extended ";
    cout << "by adding more child classes to the currently defined base classes.";
    cout << "The solver will be guiding you to specify your heat problem.\n Note: If the test case is requested, the problem definition";
    cout << " will be restrictied to a rectangular domain with equal boundary conditions for all sides except for the top side. ";
    cout << "The solution will also be limited to the homogenous case.\n\n";
    cout << "Now, let's start with defining the problem values!\n\n";
}


int main () {
    // Print out the welcome messages and program direction messages
    initiateProgram();

    // using some functionalities from the Utilities class
    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();
    
    // Storing the messages for the request to come
    vector<string> message = {"* Please enter desired coordinate system (Cartesian = 1, Polar = 2): ",
                        "\n* Run test case? (1 or 0): ",
                        "\n* Would you like to solve for the homogenous case (1 or 0): ",
                        "\n* Please enter desired solution method (LU = 0, LU sparse = 1, Gauss-Seidel = 2): "};
    
    // Define if test case is required
    bool is_test_case = static_cast<bool>(utils->requestInput('i', 0, 1, message[1]));
    
    // Define coordinate system
    // Test case is limited to the Cartesian coordinate system only
    CoordinateSystem coordinate_sys;
    if(is_test_case) coordinate_sys = CoordinateSystem::Cartesian;
    else coordinate_sys = static_cast<CoordinateSystem>(utils->requestInput('i', 1, 2, message[0]));
    
    // Define Laplace or Poisson's Equation
    // Test case AND polar case are limited to the Homogeneous case only for simplicity
    bool is_homog;
    if (!(is_test_case || coordinate_sys == CoordinateSystem::Polar)) is_homog = static_cast<bool>(utils->requestInput('i', 0, 1, message[2]));
    else is_homog = true;

    // By default, initialize Elliptic PDE problem.
    // To extend this solver, the user can specifiy the type of PDE to be solved
    std::unique_ptr<Initiation> Heat_2D = std::make_unique<EllipticPDE>(coordinate_sys, is_test_case, is_homog);
    
    // Initialize the Domain
    std::unique_ptr<Domain> domain = std::make_unique<Domain>(Heat_2D);
    
    // If inhomogenity is specified, then define the respective values
    if(!is_homog) Heat_2D->setInHomogeneous(domain);

    // Define the types and values of the boundary conditions
    Heat_2D->setBCs();

    // Create the mesh based on the specified domain 
    std::unique_ptr<Mesh> mesh = std::make_unique<Mesh>(Heat_2D, domain);

    // utils->print_matrix(mesh->getMesh());

    // Select the solution method
    int solution_method;
        // classic LU solver only implemented for the polar case
    if(Heat_2D->getCoordinateSystem() != CoordinateSystem::Polar) {
        solution_method = utils->requestInput('i', 0, 2, message[3]);
    }
    else {
        solution_method = 0;
    }
    
    
    // Create the solver object
    std::unique_ptr<Solver> solver;
    if (solution_method == 0) {
        solver = std::make_unique<LU_direct>(is_test_case);
    }
    else if (solution_method == 1) {
        solver = std::make_unique<LU_sparse>(is_test_case);
    }
    else {
        solver = std::make_unique<Seidel>(is_test_case);
    }
    
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    
    start = std::chrono::system_clock::now();
    solver->solve(Heat_2D, mesh);
    solver->setError(Heat_2D, mesh, domain);
    end = std::chrono::system_clock::now();
    
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "** Solver Elapsed time: " << elapsed_time.count() << " ms" << std::endl;


    start = std::chrono::system_clock::now();
    // Post-process the solution: printing error if available, saving results and plotting the solution.
    std::unique_ptr<PostProcessing> postproc = std::make_unique<PostProcessing>(Heat_2D,mesh,solver);
    postproc->printError();
    postproc->exportResult();

    // don't plot for the oval case
    if (!(domain->getShape() == Shape::Oval)) {
        postproc->plotResult();
    }
    end = std::chrono::system_clock::now();
    elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "** Post-processing Elapsed time: " << elapsed_time.count() << " ms" << std::endl;
}