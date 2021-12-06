#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <tuple>

#include "utils.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cin;


template <typename T> void print_vector(const T& v) {
    for (const auto& elem: v) {
        cout << elem << " ";
    }
    cout << endl;
}

template <typename T> void print_matrix(const T& m) {
    for (const auto& vec : m ) {
        print_vector(vec);
    }
}


std::tuple<vector<double>, vector<int>, bool, vector<double>> requestInputs () {
    // initializing the input values
    vector<double> dimensions = {0,0};
    vector<int> nodes = {0,0};
    bool isTestCase = false;
    int NotempsRequired = 4;
    vector<double> boundary_conditions = {0,0,0,0};

    // variables used for the automation of checking the input correctness (type and size)
    bool type_error = true;
    bool condition = true;

    // holding the message for each group of inputs (for efficiency purposes)
    vector<string> messages = {"Width (must be greater than 0): ", "Height (must be greater than 0): ",
                    "Nodes in X direction (must be integer and greater than 1): ", 
                    "Nodes in Y direction (must be integer and greater than 1): "};

    // ask for the first couple of inputs:  the dimensions of the domain and the nodes in each direction
    cout << "* Please specify the size of the rectangular domain (Width x Height)\n";
    for(int i = 0; i < 4; i++) {
        if (i == 2) { 
            cout << ">> Width: " << dimensions[0] << ", Height: " << dimensions[1] << endl;
            cout << "* Please specify the number of nodes in each direction\n"; 
        }

        condition  = true;
        type_error = true;
        while (condition || type_error) {
            cout << messages[i];

            if (i <= 1) { cin >> dimensions[i]; }
            else { cin >> nodes[i-2]; }
            
            type_error = cin.fail();
            if (i <= 1) { condition = (dimensions[i] <= 0); }
            else { condition = (nodes[i-2] <= 1); }
            
            // clear current input
            cin.clear();
            //clear buffer before taking a new line
            cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
    cout << ">> nodes in X: " << nodes[0] << ", in Y: " << nodes[1] << endl;
    
    // ask for third input: run a unit test case?
    type_error = true;
    while (type_error) {
        cout << "* Run unit test? (1/0): ";
        cin >> isTestCase;
        type_error = cin.fail();
        cin.clear();
        cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if (isTestCase) {
        cout << "You've chosen to run the test case. In this case, there are only 2 temperature values that are required. ";
        cout << "The first temperature value is for all sides of the domain except the top one. The second value is ";
        cout << "for the top side only.\n";
        NotempsRequired = 2;
    } else {

    }

    // ask for the temperature at the boundaries that are required
    // TODO implement asking for each temperature condition (Lower, Left, Right, Upper) depending on Coordinate System
    vector<string> messages_bc = {"Lower boundary - ", "Left boundary - ",
                    "Right boundary - ", 
                    "Upper boundary - "};


    for (int i = 0; i < NotempsRequired; i++) {
        type_error = true;
        while (type_error) {
            if (NotempsRequired == 4) {
                cout << messages_bc[i];
            }
            cout << "Temperature T" << i+1 << " (in degree Celsius): ";

            if (NotempsRequired == 2 && i == 0) {
                cin >> boundary_conditions[0];
                boundary_conditions[1] = boundary_conditions[0];
                boundary_conditions[2] = boundary_conditions[0];
            }
            else if (NotempsRequired == 2 && i == 1) {
                cin >> boundary_conditions[3];
            }
            else {
                cin >> boundary_conditions[i];
            }
            type_error = cin.fail();
            cin.clear();
            cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }

    return std::make_tuple(dimensions, nodes, isTestCase, boundary_conditions);
}

vector<vector<double>> vec_2_mat(size_t n_x, size_t n_y, const vector<double>& x) {
// This function converts a vector to a 2D matrix.
    
    vector<double> v (n_x,0);
    vector<vector<double>> M (n_y, v);
    return M;

}

void initiateProgram () {
    cout << "PDE Solver started!\n";
    cout << "The current version of the solver is capable of solving the following 2D heat transfer problem. For a rectangular domain with customized ";
    cout << "width and height in the cartesian coordinate system (i.e. x and y), the steady-state temperature values for the internal nodes are evaluated. ";
    cout << "Only Dirichlet boundary conditions are allowed for homogeneous materials (only a single material for the whole domain is expected). ";
    cout << "Boundary conditions are defined as follows.\n";
    cout << "Left side: Temperature T1\n";
    cout << "Top side: Temperature T2\n";
    cout << "Right side: Temperature T3\n";
    cout << "Bottom side: Temperature T4\n";
    cout << "PS:  the left-bottom corner of the domain conincides with the origion of the coordinate system.\n\n";
    cout << "Now, let's start with defining the problem values!\n\n";
}
