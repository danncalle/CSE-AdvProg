#include <iostream>
#include <fstream>

#include "headers/PostProcessing.h"

namespace plt = matplotlibcpp;

using std::cout;
using std::endl;

PostProcessing::PostProcessing(const std::unique_ptr<Initiation> & pde_type, const std::unique_ptr<Mesh> & mesh, const std::unique_ptr<Solver> & solver)
: _pde_type(pde_type.get()), _mesh(mesh->getMesh()), _sol(solver->getSolution())
{}


void PostProcessing::exportResult() {
    // Save all results collected from all previous computations to a csv file

    // get the size of the system
    const size_t rows = _mesh.size();
    const size_t mesh_cols = _mesh[0].size();

    // Open file stream and create output "results.csv" file
    std::ofstream resultfile;
    resultfile.open ("../results/results.csv");

    // Check if resultfile connection is ready, otherwise stop this function 
    if(!resultfile) {
        cout << "Couldn't establish connection to file system. Saving aborted!" << endl;
        return;
    }


    // Table headers
    resultfile << "Node number,BC (y/n),X,Y,T_numerical";
 
    if (_pde_type->isTestCase()) {
        resultfile << ",T_analytical,Error\n";
    }
    else { resultfile << "\n"; }

    // Actual data
    for (size_t i = 0; i < rows; i++) {
        // node number
        resultfile << i+1;
        // mesh details
        for (size_t j = 0; j < mesh_cols; j++) {
            resultfile << "," << _mesh[i][j];
        }
        // temperature values from the numerical computation
        resultfile << "," << _sol[i][0];
        
        // if test case is requested, output the analytical solution and the error percentages
        if (_pde_type->isTestCase()) {
            resultfile << "," << _sol[i][1] << "," << _sol[i][2];
        }
        // End row
        resultfile << "\n";
    }

    // close the file stream and ouput termination message
    resultfile.close();
    cout << "\n-------------\nComputation done!\n";
    cout << "Results can be found in results/results.csv file!\n";
}

