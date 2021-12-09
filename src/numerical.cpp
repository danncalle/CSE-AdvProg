#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <fstream>
#include <ios>
#include <limits>
#include <string>
#include <tuple>
#include <utility>

#include "./headers/numerical.h"
#include "./headers/utils.h"

using std::vector;
using std::cout;
using std::endl;
using std::cin;
using std::accumulate;

vector <double> analytical_sol (vector<vector<double>>& mesh, double W, double H, double T1, double T2) {
    // Calculates the analytical solution of the unit test case based on the desired mesh and BCs

    size_t no_of_nodes = mesh.size();
    vector<double> temperature (no_of_nodes);
    double f_sum = 0;

    for (int i = 0; i<no_of_nodes; i++) {
        // directly use the input BC values for the boundary nodes
        // case distinction applied for the corner nodes (top cornes nodes should belong to the left and right sides)
        if (mesh[i][0] && mesh[i][2] == H && mesh[i][1] != 0 && mesh[i][1] != W) {
            temperature[i] = T2;
        }
        else if (mesh[i][0]) {
            temperature[i] = T1;
        }
        // if not a boundary node, use the closed-form formula to calculate the exact temperature value
        else {
           f_sum = 0;
            for (int n = 1; n<= 200; n++) {
                f_sum += (pow(-1,n+1) +1)/n * sin(n*M_PI*mesh[i][1]/W) * sinh(n*M_PI*mesh[i][2]/W) / sinh(n*M_PI*H/W);
            }
            temperature[i] = (2/M_PI)*f_sum*(T2-T1)+T1; 
        }   
    }
    

    return temperature;
}

vector<vector <double>> generate_mesh (double x, double y, size_t n_x, size_t n_y, int coordinate_sys) {
    //  Generate a mesh based on the major dimensions of the domain and the number of nodes in each direction

    // calculating the step size and total number of nodes
    double dx = x/(n_x-1);
    double dy = y/(n_y-1);
    size_t total_nodes = n_x * n_y;

    // Initializing the rows and columns of the mesh matrix
    vector<double> rows (3);
    vector<vector<double>> mesh (total_nodes, rows);

    int current_node = 0;
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            
            // Specify if the node is at the boundary or not (if yes, then directly used the x or y values)
            // Specify the coordinates of the current node
            if (i == 0) {
                mesh[current_node][0] = 1;
                mesh[current_node][1] = 0;
                mesh[current_node][2] = i*dy;
            }
            else if (i == n_x - 1) {
                mesh[current_node][0] = 1;
                mesh[current_node][1] = x;
                mesh[current_node][2] = i*dy;
            }
            else if (j == 0) {
                mesh[current_node][0] = 1;
                mesh[current_node][1] = j*dx;
                mesh[current_node][2] = 0;
            }
            else if (j == n_y - 1) {
                mesh[current_node][0] = 1;
                mesh[current_node][1] = j*dx;
                mesh[current_node][2] = y;
            }
            else {
                mesh[current_node][0] = 0;
                mesh[current_node][1] = j*dx;
                mesh[current_node][2] = i*dy;
            }

            current_node++;
        }
    }

    return mesh;
}

vector<vector<double>> fdm_mesh(size_t n_x, size_t n_y, int coordinate_sys) {
    
    size_t no_of_elems = n_x*n_y;
    vector<double> rows (no_of_elems,0.0);
    vector<vector<double>> fdm_matrix (no_of_elems,rows);
    
    int row;

// Inner points of the mesh (not on boundaries)
    for (int i=1; i < n_x-1; i++) {
        for (int j=1; j < n_y-1; j++) {
            row = i + j*n_x;
            fdm_matrix[row][row] = -4;
            fdm_matrix[row][row-1] = 1;
            fdm_matrix[row][row+1] = 1;
            fdm_matrix[row][row+n_x] = 1;
            fdm_matrix[row][row-n_x] = 1;
        }
    }

// Lower boundary of mesh
    int j = 0;
    for (int i=0; i<n_x;i++) {
        row = i + j*n_x;
        fdm_matrix[row][row] = 1;
    }

// Left boundary of mesh
    int i = 0;
    for (int j=0; j<n_y;j++) {
        row = i + j*n_x;
        fdm_matrix[row][row] = 1;
    }

// Right boundary of mesh
    i = n_x-1;
    for (int j=0; j<n_y;j++) {
        row = i + j*n_x;
        fdm_matrix[row][row] = 1;
    }

// Upper boundary of mesh
    j = n_y-1;
    for (int i=0; i<n_x;i++) {
        row = i + j*n_x;
        fdm_matrix[row][row] = 1;
    }

    return fdm_matrix;
}

// Receives a vector of boundary conditions where in case 
// of cartesian coordinates, has the shape: [Lower BC, Left BC, Right BC, Upper BC]
// TODO: describe the shape in case of different coordinates system.

vector<double> fdm_rhs(size_t n_x, size_t n_y, int coordinate_sys,const vector<double>& bc) {
    
    size_t no_of_elems = n_x*n_y;
    vector<double> b (no_of_elems);

    int row;

// Inner points of the mesh are 0 for vector b.

// Lower boundary condition
    int j = 0;
    for (int i=0; i<n_x;i++) {
        row = i + j*n_x;
        b[row] = bc[0];
    }

// Left boundary condition
    int i = 0;
    for (int j=0; j<n_y;j++) {
        row = i + j*n_x;
        b[row] = bc[1];
    }  


// Right boundary condition
    i = n_x-1;
    for (int j=0; j<n_y;j++) {
        row = i + j*n_x;
        b[row] = bc[2];
    }  

// Upper boundary condition
    j = n_y-1;
    for (int i=0; i<n_x;i++) {
        row = i + j*n_x;
        b[row] = bc[3];
    }
    return b;
}

void LU_factorization(const vector<vector<double>>& M, vector<vector<double>>& L, vector<vector<double>>& U) {
    const int N = M.size();
    double sum;
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < i; k++) {
            
            sum = 0.0;
            for (int j = 0; j < k; j++) {
                sum = sum + L[i][j]*U[j][k];
            }
            L[i][k] = (M[i][k] - sum)/U[k][k];
        }

        L[i][i] = 1;

        for (int k = 0; k < N; k++) {
            sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum = sum + L[i][j]*U[j][k];
            }
            U[i][k] = M[i][k] - sum;
        }
        
    }
    return;
}

vector<vector<double>> test_LU(const vector<vector<double>>& L, const vector<vector<double>>& U) {

    int N = L.size();
    vector<double> v (N,0);
    vector<vector<double>> test(N, v); 
    double tmp;

    for (int i=0; i<N; i++) {

        for (int j=0; j<N; j++) {
        
        	tmp = 0;
        
            for (int k=0; k<N; k++) {
                
                tmp = tmp + L[i][k]*U[k][j];
            
            }

            test[i][j] = tmp;
        }
    }

    return test;
}

vector<double> backward_sustitution(const vector<vector<double>>& M, const vector<double>& b) {

    double sum;
    int N = M.size();
    vector<double> y (N,0);

    for (int i = N-1; i >=0; i--) {
        
        sum = 0;

        for (int j = N - 1; j > i; j--) {
            sum = sum + M[i][j]*y[j];
        }

        y[i] = (b[i] - sum)/M[i][i];
    }

    return y;

}

vector<double> forward_sustitution(const vector<vector<double>>& M, const vector<double>& b) {
    
    double sum;
    int N = M.size();
    vector<double> y (N,0);

    for (int i = 0; i < N; i++) {
        
        sum = 0;

        for (int j = 0; j < i; j++) {
            sum = sum + M[i][j]*y[j];
        }

        y[i] = (b[i] - sum)/M[i][i];
    }

    return y;
}

vector<double> finite_diff_4_pts(size_t n_x, size_t n_y, int coordinate_sys, const vector<double>& bc) {

// Creating the finite difference system for which Mx = b
    int N = n_x*n_y;
    vector<vector<double>> M = fdm_mesh(n_x, n_y, coordinate_sys);
    vector<double> b = fdm_rhs(n_x, n_y, coordinate_sys, bc);
    
// Creating the skeleton of the L and U matrices 
    vector<double> v (N);
    vector<vector<double>> L (N, v);
    vector<vector<double>> U (N, v);
    
// Solution of the system: LU factorization
    LU_factorization(M,L,U);
    vector<double> y = forward_sustitution(L,b);
    vector<double> x = backward_sustitution(U,y);
    
    return x;
}

vector<vector<double>> test_fdm_error(const vector<double>& x_analytical, const vector<double>& x_fdm) {
    int N = x_analytical.size();
    vector<double> v (N, 0);
    vector<vector<double>> error_output (2, v);
    double avg;

    for (int i = 0; i < N; i++) {
        error_output[0][i] = abs(x_analytical[i] - x_fdm[i]);
        if (abs(x_analytical[i]) > pow(10,-14)) {
            error_output[1][i] = abs(x_analytical[i] - x_fdm[i])/abs(x_analytical[i]);
        }
        else {
            error_output[1][i] = 0;
        }
    }

    avg = accumulate(error_output[1].begin(), error_output[1].end(), 0.0) / N;

    cout << "\n-------------------- Error report for Unit test case: ---------------- \n" << endl;
    cout << "The maximum error is: " << *std::max_element(error_output[1].begin(), error_output[1].end()) << endl;
    cout << "The minimum error is: " << *std::min_element(error_output[1].begin(), error_output[1].end()) << endl;
    cout << "The average error is: " << avg << endl;

    return error_output;
}
