#include <iostream>
#include <fstream>

#include "./headers/PostProcessing.h"

namespace plt = matplotlibcpp;

using std::cout;
using std::endl;

PostProcessing::PostProcessing(const std::unique_ptr<Initiation> & pde_type, const std::unique_ptr<Mesh> & mesh, const std::unique_ptr<Solver> & solver)
: _pde_type(pde_type.get()), _mesh(mesh.get()), _sol(solver->getSolution())
{
    if(_pde_type->isTestCase()) {
        int N = _sol[0].size();
        _error = vector<double> (N, 0);

        for (int i = 0; i < N; i++) {
            if (abs(_sol[1][i]) > pow(10,-14)) {
                _error[i] = abs(_sol[1][i] - _sol[0][i])/abs(_sol[1][i]);
            }
            else {
                _error[i] = 0;
            }
        }    
    }
}


void PostProcessing::exportResult() {
    // @ Save all results collected from all previous computations to a csv file
    const auto& mesh = _mesh->getMesh();

    // get the size of the system
    const size_t rows = mesh.size();
    const size_t mesh_cols = mesh[0].size();

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
            resultfile << "," << mesh[i][j];
        }
        // temperature values from the numerical computation
        resultfile << "," << _sol[0][i];
        
        // if test case is requested, output the analytical solution and the error percentages
        if (_pde_type->isTestCase()) {
            resultfile << "," << _sol[1][i] << "," << _error[i];
        }
        // End row
        resultfile << "\n";
    }

    // close the file stream and ouput termination message
    resultfile.close();
    cout << "\n-------------\nSaving is done!\n";
    cout << "Results can be found in results/results.csv file!\n";
}

void PostProcessing::printError(){
    // @ Print a summary of the error of the numerical against analytical solution if test case is true.
    if(_pde_type->isTestCase()){
        double avg = accumulate(_error.begin(), _error.end(), 0.0) / _error.size();
        cout << "\n-------------------- Error report for Unit test case: ---------------- \n" << endl;
        cout << "The maximum error is: " << *std::max_element(_error.begin(), _error.end()) << endl;
        cout << "The minimum error is: " << *std::min_element(_error.begin(), _error.end()) << endl;
        cout << "The average error is: " << avg << endl;
    }
    else {
        cout << "\n-------------------- Not a test case, no error output ---------------- \n" << endl;
    }
}

void PostProcessing::plotResult(){
    // @ Using matplotlibcpp, open a GUI to plot the results of the simulated temperature values on a 3D contour plot
    
    // initialize the x, y, and z matrices (and their corresponding rows)
    const array<double,2> step = _mesh->getStepSize();
    const array<int,2> nodes = _mesh->getNumNodes();
    const auto& mesh = _mesh->getTransMesh();

    size_t counter = 0;
    
    vector<vector<double>> x, y, z;
    x.reserve(nodes[1]);
    y.reserve(nodes[1]);
    z.reserve(nodes[1]);

    vector<double> x_row (nodes[0], 0), y_row (nodes[0], 0), z_row (nodes[0], 0);
    
    for (int i = 0; i < nodes[1];  i++) {
        for (int j = 0; j < nodes[0]; j++) {
            // evaluate the current x and y coordinates and the equivalent temperature value, then add to the x, y, and z rows
            x_row[j] = mesh[1][counter];
            y_row[j] = mesh[2][counter];      
            z_row[j] = _sol[0][counter];
            
            counter++;
        }
        // add these rows to the respective matrices
        x.push_back(x_row);
        y.push_back(y_row);
        z.push_back(z_row);
    }


    // repeat the first row of the solution to close the circle for the polar case 
    if(_pde_type->getCoordinateSystem() == CoordinateSystem::Polar) {
        counter = 0;
        x_row = {}, y_row = {}, z_row = {};
        for (int j = 0; j < nodes[0]; j++) {
            // evaluate the current x and y coordinates and the equivalent temperature value, then add to the x, y, and z rows
            
            x_row[counter] = mesh[1][counter];
            y_row[counter] = mesh[2][counter];      
            z_row[counter] = _sol[0][counter];
            
            counter++;
        }

        x.push_back(x_row);
        y.push_back(y_row);
        z.push_back(z_row);
    }
    

    

    // plot and show the generated surface
    plt::plot_surface(x, y, z);
    plt::show();

}
