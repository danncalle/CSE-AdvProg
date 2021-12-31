#include "./headers/Solver.h"

Solver::Solver(SolverType solver_type, bool test_case):
__solver_type(solver_type), __test_case(test_case) {};
        
vector<double> Solver::getError() const {return __error;};
vector<vector<double>> Solver::getSolution() const {return __solution;};

LU::LU(bool test_case):
Solver(SolverType::Numerical, test_case) {}

void LU::setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){
        // Right hand side of the linear system Ax = b, in includes the boundary conditions (and the heat generation term in future implementation).
        // Receives a vector of boundary conditions where in case 
        // of cartesian coordinates, has the shape: [Lower BC, Left BC, Right BC, Upper BC]
        array<int,2> nodes = mesh->getNumNodes();
        size_t no_of_elems = nodes[0]*nodes[1];
        
        vector<double> bc (4,0.0);

        if(pde->isTestCase()){
            bc[0] = pde->getBcs()[0];
            bc[1] = bc[0];
            bc[2] = bc[0];
            bc[3] = pde->getBcs()[1];
        }
        else {
            bc = pde->getBcs();
        }

        __b = vector<double> (no_of_elems, 0.0);
        int row;

    // Inner points of the mesh are 0 for vector b when no internal heat generation is present.

    // Upper boundary condition
        int j =  nodes[1]-1;
        for (int i=0; i< nodes[0];i++) {
            row = i + j* nodes[0];
            __b[row] = bc[3];
        }

    // Bottom boundary condition
        j = 0;
        for (int i=0; i< nodes[0];i++) {
            row = i + j* nodes[0];
            __b[row] = bc[0];
        }

    // Left boundary condition
        int i = 0;
        for (int j=0; j< nodes[1];j++) {
            row = i + j* nodes[0];
            __b[row] = bc[1];
        }  


    // Right boundary condition
        i =  nodes[0]-1;
        for (int j=0; j< nodes[1];j++) {
            row = i + j* nodes[0];
            __b[row] = bc[2];
        }  

        return;
    }
    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar) {
        return;
    }
};

void LU::setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){
        // Creates the matrix of the finite differences system of equations based on 
        // the number of nodes and coordinate system.

        //Calculates size of the matrix
        array<int,2> nodes = mesh->getNumNodes();
        size_t no_of_elems = nodes[0]*nodes[1];

        double dx = mesh->getStepSize()[0];
        double dy = mesh->getStepSize()[1];
        vector<double> rows (no_of_elems, 0.0);
        __M = vector<vector<double>> (no_of_elems,rows);
        __L = vector<vector<double>> (no_of_elems,rows);
        __U = vector<vector<double>> (no_of_elems,rows);

        int row;

    // Filling the elements of the matrix
    // Inner points of the mesh (not on boundaries)
        for (int i=1; i < nodes[0]-1; i++) {
            for (int j=1; j < nodes[1]-1; j++) {
                row = i + j*nodes[0];
                __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[row][row-1] = pow(dy, 2);
                __M[row][row+1] = pow(dy, 2);
                __M[row][row+nodes[0]] = pow(dx, 2);
                __M[row][row-nodes[0]] = pow(dx, 2);
            }
        }

    // Lower boundary of mesh
        int j = 0;
        for (int i=0; i<nodes[0];i++) {
            row = i + j*nodes[0];
            __M[row][row] = 1;
        }

    // Left boundary of mesh
        int i = 0;
        for (int j=0; j<nodes[1];j++) {
            row = i + j*nodes[0];
            __M[row][row] = 1;
        }

    // Right boundary of mesh
        i = nodes[0]-1;
        for (int j=0; j<nodes[1];j++) {
            row = i + j*nodes[0];
            __M[row][row] = 1;
        }

    // Upper boundary of mesh
        j = nodes[1]-1;
        for (int i=0; i<nodes[0];i++) {
            row = i + j*nodes[0];
            __M[row][row] = 1;
        }

        return;
    }
    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar) {
        return;
    }   
};

void LU::LUFactorization() {
// Perform the LU factorization of a matrix M.

    const int N = __M.size();
    double sum;

// Calculate elements of L and then elements of U 
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < i; k++) {
            
            sum = 0.0;
            for (int j = 0; j < k; j++) {
                sum = sum + __L[i][j]*__U[j][k];
            }
            __L[i][k] = (__M[i][k] - sum)/__U[k][k];
        }

        __L[i][i] = 1;

        for (int k = 0; k < N; k++) {
            sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum = sum + __L[i][j]*__U[j][k];
            }
            __U[i][k] = __M[i][k] - sum;
        }
        
    }
    return;
}

vector<double> LU::backwardSubstitution(const vector<double>& rhs){
    
    double sum;
    int N = __U.size();
    vector<double> y (N,0);

    for (int i = N-1; i >=0; i--) {
        
        sum = 0;

        for (int j = N - 1; j > i; j--) {
            sum = sum + __U[i][j]*y[j];
        }

        y[i] = (rhs[i] - sum)/__U[i][i];
    }    
    return y;
}

vector<double> LU::forwardSubstitution(const vector<double>& rhs){
// Takes a lower triangular matrix M and vector b and solves the system Mx = b via forward substitution.

    double sum;
    int N = __L.size();
    vector<double> y (N,0);

    for (int i = 0; i < N; i++) {
        
        sum = 0;

        for (int j = 0; j < i; j++) {
            sum = sum + __L[i][j]*y[j];
        }

        y[i] = (rhs[i] - sum)/__L[i][i];
    }

    return y;
}

void LU::solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh) {
// Solve a 2D poisson's equation using the 4 point stencil finite differences method.

// Creating the finite difference system for which Mx = b
    LU::setupMatrix(pde, mesh);
    LU::setupRhs(pde, mesh);

// Solution of the system: LU factorization
    LUFactorization();
    vector<double> y = forwardSubstitution(__b);
    __solution[0] = backwardSubstitution(y);

    return;
}

void LU::setError(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain){
    if(!__test_case){
        int N = __solution[0].size();
        __error = vector<double> (N, 0);
        return;
    }
    else {
        std::unique_ptr<Analytic> analytic = std::make_unique<Analytic>();
        analytic->solve(pde, mesh, domain);
        __solution[1] = analytic->getAnalyticSol();

        int N = __solution[0].size();
        __error = vector<double> (N, 0);
        double avg;

        for (int i = 0; i < N; i++) {
            if (abs(__solution[1][i]) > pow(10,-14)) {
                __error[i] = abs(__solution[1][i] - __solution[0][i])/abs(__solution[1][i]);
            }
            else {
                __error[i] = 0;
            }
        }

        // avg = accumulate(__error.begin(), __error.end(), 0.0) / N;

        // cout << "\n-------------------- Error report for Unit test case: ---------------- \n" << endl;
        // cout << "The maximum error is: " << *std::max_element(__error.begin(), __error.end()) << endl;
        // cout << "The minimum error is: " << *std::min_element(__error.begin(), __error.end()) << endl;
        // cout << "The average error is: " << avg << endl;

        return;
    }
}

Analytic::Analytic():
Solver(SolverType::Analytical, true) {}

void Analytic::solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain){
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){
    // Calculates the analytical solution of the unit test case based on the desired mesh and BCs

        size_t no_of_nodes = mesh->getTotalNodes();
        double T1 = pde->getBcs()[0];
        double T2 = pde->getBcs()[1];
        double W = domain->getDimensions()[0];
        double H = domain->getDimensions()[1];

        __analytic_sol = vector<double>(no_of_nodes,0);
        vector<vector<double>> mesh_matrix = mesh->getMesh();
        double f_sum = 0;

        for (int i = 0; i<no_of_nodes; i++) {
            // directly use the input BC values for the boundary nodes
            // case distinction applied for the corner nodes (top cornes nodes should belong to the left and right sides)
            if (mesh_matrix[i][0] && (mesh_matrix[i][2] == 0 || mesh_matrix[i][1] == 0 || mesh_matrix[i][1] == W)) {
                __analytic_sol[i] = T1;
            }
            else if (mesh_matrix[i][0] && mesh_matrix[i][2] == H && mesh_matrix[i][1] != 0 && mesh_matrix[i][1] != W) {
                __analytic_sol[i] = T2;
            }
            // if not a boundary node, use the closed-form formula to calculate the exact __solution value
            else {
            f_sum = 0;
                for (int n = 1; n<= 200; n++) {
                    f_sum += (pow(-1,n+1) +1)/n * sin(n*M_PI*mesh_matrix[i][1]/W) * sinh(n*M_PI*mesh_matrix[i][2]/W) / sinh(n*M_PI*H/W);
                }
                __analytic_sol[i] = (2/M_PI)*f_sum*(T2-T1)+T1; 
            }   
        }

        return;
    }
    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar){
        return;
    }

};

vector<double> Analytic::getAnalyticSol() const {return __analytic_sol;}
void Analytic::setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){
    return;
};
void Analytic::setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){
    return;
};
void Analytic::setError(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain){
    return;
};
void Analytic::solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){
    return;
};
