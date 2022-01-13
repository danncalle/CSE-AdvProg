#include "./headers/Solver.h"

Solver::Solver(SolverType solver_type, bool test_case):
__solver_type(solver_type), __test_case(test_case) {
    cout << "Solver initiated correctly\n";
};
        
vector<double> Solver::getError() const {return __error;};
vector<vector<double>> Solver::getSolution() const {return __solution;};

void Solver::solve_analytical(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain){
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){
    // Calculates the analytical solution of the unit test case based on the desired mesh and BCs

        size_t no_of_nodes = mesh->getTotalNodes();
        double T1 = pde->getBcs()[0];
        double T2 = pde->getBcs()[1];
        double W = domain->getDimensions()[0];
        double H = domain->getDimensions()[1];

        __solution[1] = vector<double>(no_of_nodes,0);
        vector<vector<double>> mesh_matrix = mesh->getMesh();
        double f_sum = 0;

        for (int i = 0; i<no_of_nodes; i++) {
            // directly use the input BC values for the boundary nodes
            // case distinction applied for the corner nodes (top cornes nodes should belong to the left and right sides)
            if (mesh_matrix[i][0] && (mesh_matrix[i][2] == 0 || mesh_matrix[i][1] == 0 || mesh_matrix[i][1] == W)) {
                __solution[1][i] = T1;
            }
            else if (mesh_matrix[i][0] && mesh_matrix[i][2] == H && mesh_matrix[i][1] != 0 && mesh_matrix[i][1] != W) {
                __solution[1][i] = T2;
            }
            // if not a boundary node, use the closed-form formula to calculate the exact __solution value
            else {
            f_sum = 0;
                for (int n = 1; n<= 200; n++) {
                    f_sum += (pow(-1,n+1) +1)/n * sin(n*M_PI*mesh_matrix[i][1]/W) * sinh(n*M_PI*mesh_matrix[i][2]/W) / sinh(n*M_PI*H/W);
                }
                __solution[1][i] = (2/M_PI)*f_sum*(T2-T1)+T1; 
            }   
        }

        return;
    }
    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar){
        return;
    }

};

void Solver::setError(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain){
    if(!__test_case){
        int N = __solution[0].size();
        __error = vector<double> (N, 0);
        return;
    }
    else {
        solve_analytical(pde, mesh, domain);
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

LU_direct::LU_direct(const bool test_case):
Solver(SolverType::LU_dir, test_case) {
    cout << "LU_direct solver initiated correctly\n";
}

void LU_direct::setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){
        // Right hand side of the linear system Ax = b, in includes the boundary conditions (and the heat generation term in future implementation).
        // Receives a vector of boundary conditions where in case 
        // of cartesian coordinates, has the shape: [Lower BC, Left BC, Right BC, Upper BC]
        array<int,2> nodes = mesh->getNumNodes();
        size_t no_of_elems = nodes[0]*nodes[1];
        vector<double> bc (4,0.0);
        vector<BoundaryTypes> bc_types = pde->getBcTypes();
        __b = vector<double> (no_of_elems, 0.0);
        double heat_ss = pde->getHeatSs();
        int row;

        // for (auto &elem: bc_types){
        //     cout << static_cast<std::underlying_type<BoundaryTypes>::type>(elem)<< endl;
        // }

        if(pde->isTestCase()){
            bc[0] = pde->getBcs()[0];
            bc[1] = bc[0];
            bc[2] = bc[0];
            bc[3] = pde->getBcs()[1];
        }
        else {
            bc = pde->getBcs();
        }

    // Inner points of the mesh are 0 for vector b when no internal heat generation is present.

    // Upper boundary condition
        int j =  nodes[1]-1;
        if (bc_types[3] == BoundaryTypes::Dirichlet) {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b[row] = bc[3];
            }
        }
        else {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b[row] = bc[3]*mesh->getStepSize()[1];
            }            
        }

    // Bottom boundary condition
        j = 0;
        if (bc_types[0] == BoundaryTypes::Dirichlet) {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b[row] = bc[0];
            }
        }
        else {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b[row] = bc[0]*mesh->getStepSize()[1];
            }
        }

    // Left boundary condition
        int i = 0;
        if (bc_types[1] == BoundaryTypes::Dirichlet) {
            for (int j=1; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b[row] = bc[1];
            }  
        }
        else {
            for (int j=0; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b[row] = bc[1]*mesh->getStepSize()[0];
            }              
        }


    // Right boundary condition
        i =  nodes[0]-1;
        if (bc_types[2] == BoundaryTypes::Dirichlet) {
            for (int j=0; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b[row] = bc[2];
            }  
        }
        else {
            for (int j=0; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b[row] = bc[2]*mesh->getStepSize()[0];
            }              
        }

        if (heat_ss!=0){ // Accounting only for point heat generation/sink
            double k = pde->getK();
            array<double,2> step = mesh->getStepSize();
            array<int,2> heat_ss_node_loc;
            heat_ss_node_loc[0] = (int) round(pde->getHeatSsLoc()[0]/step[0]);
            heat_ss_node_loc[1] = (int) round(pde->getHeatSsLoc()[1]/step[1]);
            __b[heat_ss_node_loc[0] + heat_ss_node_loc[1]*nodes[0]] += - heat_ss/k * pow(step[0]*step[1], 2);
        }

        return;
    }
    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar) {
        return;
    }
};

void LU_direct::setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){
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
        vector<BoundaryTypes> bc_types = pde->getBcTypes();

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
        if (bc_types[0] == BoundaryTypes::Dirichlet) {
            for (int i=0; i<nodes[0];i++) {
                row = i + j*nodes[0];
                __M[row][row] = 1;
            }
        }
        else {
            // corner node at i,j = 0,0 
            if (bc_types[1] == BoundaryTypes::Neumann) {
                __M[0][0] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[0][1] = 2 * pow(dy, 2);
                __M[0][nodes[0]] = 2 * pow(dx, 2);
            }
            else {
                __M[0][0] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[0][1] = pow(dy, 2);
                __M[0][nodes[0]] = 2 * pow(dx, 2);
            }
            // corner node at i,j = Nx-1,0 
            if (bc_types[2] == BoundaryTypes::Neumann) {
                __M[nodes[0]-1][nodes[0]-1] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[nodes[0]-1][nodes[0]-2] = 2 * pow(dy, 2);
                __M[nodes[0]-1][2*nodes[0]-1] = 2 * pow(dx, 2);
            }
            else {
                __M[nodes[0]-1][nodes[0]-1] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[nodes[0]-1][nodes[0]-2] = pow(dy, 2);
                __M[nodes[0]-1][2*nodes[0]-1] = 2 * pow(dx, 2);
            }
            // non-corner nodes
            for (int i=1; i<nodes[0]-1;i++) {
                row = i + j*nodes[0];
                __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[row][row-1] = pow(dy, 2);
                __M[row][row+1] = pow(dy, 2);
                __M[row][row+nodes[0]] = 2 * pow(dx, 2);
            }
        }

    // Left boundary of mesh
        int i = 0;
        if (bc_types[1] == BoundaryTypes::Dirichlet) {
            for (int j=0; j<nodes[1];j++) {
                row = i + j*nodes[0];
                __M[row][row] = 1;
            }
        }
        else {
            // corner node at i,j = 0,Ny-1
            row = (nodes[1]-1)*nodes[0];
            if (bc_types[3] == BoundaryTypes::Neumann) {
                __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[row][row+1] = 2 * pow(dy, 2);
                __M[row][row - nodes[0]] = 2 * pow(dx, 2);
            }
            else {
                __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[row][row+1] = pow(dy, 2);
                __M[row][row - nodes[0]] = 2 * pow(dx, 2);
            }
            // non-corner nodes
            for (int j=1; j<nodes[1]-1;j++) {
                row = i + j*nodes[0];
                __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[row][row-1] = pow(dy, 2);
                __M[row][row+1] = pow(dy, 2);
                __M[row][row-nodes[0]] = 2 * pow(dx, 2);
            }
        }

    // Right boundary of mesh
        i = nodes[0]-1;
        if (bc_types[2] == BoundaryTypes::Dirichlet) {
            for (int j=0; j<nodes[1];j++) {
                row = i + j*nodes[0];
                __M[row][row] = 1;
            }
        }
        else {
            // corner node at i,j = Nx-1,Ny-1
            row = nodes[0]*nodes[1] - 1;
            if (bc_types[3] == BoundaryTypes::Neumann) {
                __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[row][row-1] = 2 * pow(dy, 2);
                __M[row][row - nodes[0]] = 2 * pow(dx, 2);
            }
            else {
                __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[row][row-1] = pow(dy, 2);
                __M[row][row - nodes[0]] = 2 * pow(dx, 2);
            }
            // non-corner nodes
            for (int j=1; j<nodes[1]-1;j++) {
                row = i + j*nodes[0];
                __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[row][row-1] = pow(dy, 2);
                __M[row][row+1] = pow(dy, 2);
                __M[row][row-nodes[0]] = 2 * pow(dx, 2);
            }
        }

    // Upper boundary of mesh
        j = nodes[1]-1;
        if (bc_types[3] == BoundaryTypes::Dirichlet) {
            for (int i=0; i<nodes[0];i++) {
                row = i + j*nodes[0];
                __M[row][row] = 1;
            }
        }
        else {
            // all corner nodes were already assigned
            // non-corner nodes
            for (int i=1; i<nodes[0]-1;i++) {
                row = i + j*nodes[0];
                __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                __M[row][row-1] = pow(dy, 2);
                __M[row][row+1] = pow(dy, 2);
                __M[row][row-nodes[0]] = 2 * pow(dx, 2);
            }
        }
        return;
    }
    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar) {
        return;
    }   
};

void LU_direct::LUFactorization() {
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

vector<double> LU_direct::backwardSubstitution(const vector<double>& rhs){
    
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

vector<double> LU_direct::forwardSubstitution(const vector<double>& rhs){
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

void LU_direct::solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh) {
// Solve a 2D poisson's equation using the 4 point stencil finite differences method.
    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();
// Creating the finite difference system for which Mx = b
    LU_direct::setupMatrix(pde, mesh);
    LU_direct::setupRhs(pde, mesh);

    cout<< "Matrix M\n";
    utils->print_matrix(__M);
    cout<< "vector b\n";
    utils->print_vector(__b);


// Solution of the system: LU factorization
    LUFactorization();
    vector<double> y = forwardSubstitution(__b);
    __solution[0] = backwardSubstitution(y);

    cout<< "solution T\n";
    utils->print_vector(__solution[0]);

    return;
}

void LU_direct::writeLU() {
    // Save all results collected from all previous computations to a csv file

    // get the size of the system
    const size_t rows = __L.size();
    const size_t cols = __L[0].size();

    // Open file stream and create output "L_direct.csv" and "U_direct" files
    std::ofstream L;
    L.open ("../results/L_direct.csv");

    std::ofstream U;
    U.open ("../results/U_direct.csv");

    // Check if resultfile connection is ready, otherwise stop this function 
    if(!L || !U) {
        cout << "Couldn't establish connection to file system. Saving aborted!" << endl;
        return;

    }
    
    for (int i=0; i<rows; i++) {
        for (int j=0; j<cols; j++) {
            if(i == j){
                    L << __L[i][j] << "\n";            
                    U << __U[i][j] << ",";
            }
            else if(i < j) {
                U << __U[i][j] << ",";
                if (j==cols-1) U << __U[i][j] << "\n";
            }
            else if(i > j) {
                L << __L[i][j] << ",";
            }
        }
    }
    L.close();
    U.close();
    std::cout << "\n-------------\nL and U matrixes were successfully stored in the results folder\n";
    return;
};

void LU_direct::loadLU() {
    // Save all results collected from all previous computations to a csv file
    // Open file stream and create output "L_direct.csv" and "U_direct" files
    std::ifstream L;
    L.open ("../results/L_direct.csv");

    std::ifstream U;
    U.open ("../results/U_direct.csv");

    size_t N = __b.size();
    vector<double> rows (N, 0.0);
    __L = vector<vector<double>> (N,rows);
    __U = vector<vector<double>> (N,rows);

    int i=0, j=0;
    vector<string> row;
    string line, word, temp;
    std::stringstream ss;

    if(!L.is_open() || !U.is_open()) {
        cout << "Couldn't establish connection to file system. Loading aborted!" << endl;
        return;
    }

    while (L >> temp) {
        row.clear();
        ss.clear();

// Reading a line and convert it into stringstream (splits it into words)
        std::getline(L,line,'\n');
        ss.str(line);

// Asign each part of the line to the variable word and converts it into a double
        while (std::getline(ss, word, ',')){
            __L[i][j] = std::stod(word);
            j++;
        }
        i++;
        j=0;
    }
    // Restarting counters to fill matrix
    i=0,j=0;
    while (U >> temp) {
        row.clear();
        ss.clear();

        std::getline(U,line,'\n');
        ss.str(line);

        while (std::getline(ss, word, ',')){
            __L[i][j] = std::stod(word);
            j++;
        }
        i++;
        j=0;
    }

    return;
};

LU_sparse::LU_sparse(bool test_case):
Solver(SolverType::LU_sp, test_case){
    cout << "LU_sparse solver initiated correctly\n";
};

void LU_sparse::setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){
        // Right hand side of the linear system Ax = b, in includes the boundary conditions (and the heat generation term in future implementation).
        // Receives a vector of boundary conditions where in case 
        // of cartesian coordinates, has the shape: [Lower BC, Left BC, Right BC, Upper BC]
        array<int,2> nodes = mesh->getNumNodes();
        size_t no_of_elems = nodes[0]*nodes[1];
        vector<double> bc (4,0.0);
        vector<BoundaryTypes> bc_types = pde->getBcTypes();
        __b = VectorXd (no_of_elems);
        double heat_ss = pde->getHeatSs();
        int row;

        // for (auto &elem: bc_types){
        //     cout << static_cast<std::underlying_type<BoundaryTypes>::type>(elem)<< endl;
        // }

        if(pde->isTestCase()){
            bc[0] = pde->getBcs()[0];
            bc[1] = bc[0];
            bc[2] = bc[0];
            bc[3] = pde->getBcs()[1];
        }
        else {
            bc = pde->getBcs();
        }

    // Inner points of the mesh are 0 for vector b when no internal heat generation is present.

    // Upper boundary condition
        int j =  nodes[1]-1;
        if (bc_types[3] == BoundaryTypes::Dirichlet) {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b.coeffRef(row) = bc[3];
            }
        }
        else {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b.coeffRef(row) = bc[3]*mesh->getStepSize()[1];
            }            
        }

    // Bottom boundary condition
        j = 0;
        if (bc_types[0] == BoundaryTypes::Dirichlet) {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b.coeffRef(row) = bc[0];
            }
        }
        else {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b.coeffRef(row) = bc[0]*mesh->getStepSize()[1];
            }
        }

    // Left boundary condition
        int i = 0;
        if (bc_types[1] == BoundaryTypes::Dirichlet) {
            for (int j=1; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b.coeffRef(row) = bc[1];
            }  
        }
        else {
            for (int j=0; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b.coeffRef(row) = bc[1]*mesh->getStepSize()[0];
            }              
        }


    // Right boundary condition
        i =  nodes[0]-1;
        if (bc_types[2] == BoundaryTypes::Dirichlet) {
            for (int j=0; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b.coeffRef(row) = bc[2];
            }  
        }
        else {
            for (int j=0; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b.coeffRef(row) = bc[2]*mesh->getStepSize()[0];
            }              
        }

        if (!pde->isHomogeneous()){ // Accounting only for point heat generation/sink
            double k = pde->getK();
            array<double,2> step = mesh->getStepSize();
            array<int,2> heat_ss_node_loc;
            heat_ss_node_loc[0] = (int) round(pde->getHeatSsLoc()[0]/step[0]);
            heat_ss_node_loc[1] = (int) round(pde->getHeatSsLoc()[1]/step[1]);
            __b.coeffRef(heat_ss_node_loc[0] + heat_ss_node_loc[1]*nodes[0]) += - heat_ss/k * pow(step[0]*step[1], 2);
        }

        return;
    }
    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar) {
        return;
    }
};

void LU_sparse::setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){
        // Creates the matrix of the finite differences system of equations based on 
        // the number of nodes and coordinate system.

        //Calculates size of the matrix
        array<int,2> nodes = mesh->getNumNodes();
        size_t no_of_elems = nodes[0]*nodes[1];

        double dx = mesh->getStepSize()[0];
        double dy = mesh->getStepSize()[1];
        vector<double> rows (no_of_elems, 0.0);
        __M = SparseMatrix<double> (no_of_elems,no_of_elems);
        vector<BoundaryTypes> bc_types = pde->getBcTypes();

        int row;

    // Filling the elements of the matrix
    // Inner points of the mesh (not on boundaries)
        for (int i=1; i < nodes[0]-1; i++) {
            for (int j=1; j < nodes[1]-1; j++) {
                row = i + j*nodes[0];
                __M.coeffRef(row,row) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(row,row-1) = pow(dy, 2);
                __M.coeffRef(row,row+1) = pow(dy, 2);
                __M.coeffRef(row,row+nodes[0]) = pow(dx, 2);
                __M.coeffRef(row,row-nodes[0]) = pow(dx, 2);
            }
        }

    // Lower boundary of mesh
        int j = 0;
        if (bc_types[0] == BoundaryTypes::Dirichlet) {
            for (int i=0; i<nodes[0];i++) {
                row = i + j*nodes[0];
                __M.coeffRef(row,row) = 1;
            }
        }
        else {
            // corner node at i,j = 0,0 
            if (bc_types[1] == BoundaryTypes::Neumann) {
                __M.coeffRef(0,0) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(0,1) = 2 * pow(dy, 2);
                __M.coeffRef(0,nodes[0]) = 2 * pow(dx, 2);
            }
            else {
                __M.coeffRef(0,0) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(0,1) = pow(dy, 2);
                __M.coeffRef(0,nodes[0]) = 2 * pow(dx, 2);
            }
            // corner node at i,j = Nx-1,0 
            if (bc_types[2] == BoundaryTypes::Neumann) {
                __M.coeffRef(nodes[0]-1,nodes[0]-1) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(nodes[0]-1,nodes[0]-2) = 2 * pow(dy, 2);
                __M.coeffRef(nodes[0]-1,2*nodes[0]-1) = 2 * pow(dx, 2);
            }
            else {
                __M.coeffRef(nodes[0]-1,nodes[0]-1) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(nodes[0]-1,nodes[0]-2) = pow(dy, 2);
                __M.coeffRef(nodes[0]-1,2*nodes[0]-1) = 2 * pow(dx, 2);
            }
            // non-corner nodes
            for (int i=1; i<nodes[0]-1;i++) {
                row = i + j*nodes[0];
                __M.coeffRef(row,row) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(row,row-1) = pow(dy, 2);
                __M.coeffRef(row,row+1) = pow(dy, 2);
                __M.coeffRef(row,row+nodes[0]) = 2 * pow(dx, 2);
            }
        }

    // Left boundary of mesh
        int i = 0;
        if (bc_types[1] == BoundaryTypes::Dirichlet) {
            for (int j=0; j<nodes[1];j++) {
                row = i + j*nodes[0];
                __M.coeffRef(row,row) = 1;
            }
        }
        else {
            // corner node at i,j = 0,Ny-1
            row = (nodes[1]-1)*nodes[0];
            if (bc_types[3] == BoundaryTypes::Neumann) {
                __M.coeffRef(row,row) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(row,row+1) = 2 * pow(dy, 2);
                __M.coeffRef(row,row - nodes[0]) = 2 * pow(dx, 2);
            }
            else {
                __M.coeffRef(row,row) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(row,row+1) = pow(dy, 2);
                __M.coeffRef(row,row - nodes[0]) = 2 * pow(dx, 2);
            }
            // non-corner nodes
            for (int j=1; j<nodes[1]-1;j++) {
                row = i + j*nodes[0];
                __M.coeffRef(row,row) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(row,row-1) = pow(dy, 2);
                __M.coeffRef(row,row+1) = pow(dy, 2);
                __M.coeffRef(row,row-nodes[0]) = 2 * pow(dx, 2);
            }
        }

    // Right boundary of mesh
        i = nodes[0]-1;
        if (bc_types[2] == BoundaryTypes::Dirichlet) {
            for (int j=0; j<nodes[1];j++) {
                row = i + j*nodes[0];
                __M.coeffRef(row,row) = 1;
            }
        }
        else {
            // corner node at i,j = Nx-1,Ny-1
            row = nodes[0]*nodes[1] - 1;
            if (bc_types[3] == BoundaryTypes::Neumann) {
                __M.coeffRef(row,row) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(row,row-1) = 2 * pow(dy, 2);
                __M.coeffRef(row,row - nodes[0]) = 2 * pow(dx, 2);
            }
            else {
                __M.coeffRef(row,row) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(row,row-1) = pow(dy, 2);
                __M.coeffRef(row,row - nodes[0]) = 2 * pow(dx, 2);
            }
            // non-corner nodes
            for (int j=1; j<nodes[1]-1;j++) {
                row = i + j*nodes[0];
                __M.coeffRef(row,row) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(row,row-1) = pow(dy, 2);
                __M.coeffRef(row,row+1) = pow(dy, 2);
                __M.coeffRef(row,row-nodes[0]) = 2 * pow(dx, 2);
            }
        }

    // Upper boundary of mesh
        j = nodes[1]-1;
        if (bc_types[3] == BoundaryTypes::Dirichlet) {
            for (int i=0; i<nodes[0];i++) {
                row = i + j*nodes[0];
                __M.coeffRef(row,row) = 1;
            }
        }
        else {
            // all corner nodes were already assigned
            // non-corner nodes
            for (int i=1; i<nodes[0]-1;i++) {
                row = i + j*nodes[0];
                __M.coeffRef(row,row) = -2 * (pow(dx, 2) + pow(dy, 2));
                __M.coeffRef(row,row-1) = pow(dy, 2);
                __M.coeffRef(row,row+1) = pow(dy, 2);
                __M.coeffRef(row,row-nodes[0]) = 2 * pow(dx, 2);
            }
        }
        return;
    }
    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar) {
        return;
    }   
};

void LU_sparse::solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){
// Solve a 2D poisson's equation using the 4 point stencil finite differences method.
// Creating the finite difference system for which Mx = b
    LU_sparse::setupMatrix(pde, mesh);
    LU_sparse::setupRhs(pde, mesh);
    // Eigen::SimplicialLDLT<SparseMatrix<double>> solver; //Did not work with the Cholesky decomposition
    Eigen::SparseLU<SparseMatrix<double>> solver;
    solver.compute(__M);

    if (solver.info()!=Eigen::ComputationInfo::Success){
        cout << "Decomposition failed" << endl;
    }

    // cout << "Matrix L: \n" << solver.matrixL();
    // cout << "Matrix U: \n" << solver.matrixU();

    VectorXd sol = solver.solve(__b);

    if (solver.info()!=Eigen::ComputationInfo::Success){
        cout << "Solution failed" << endl;
    }
    // cout<< "Matrix M\n" << __M << endl;
    // cout<< "vector b\n" << __b << endl;

// Converting Eigen::vector to std::vector 
    // cout<<sol<<endl;
    __solution[0] = vector<double> (sol.size(),0.0);
    for (int i=0; i<sol.size(); i++){
        __solution[0][i] = sol(i);
    }

// Storing decomposition

    return;    
};

Seidel::Seidel(bool test_case):
Solver(SolverType::Gauss_Seidel, test_case){
    cout << "Gauss-Seidel solver initiated correctly\n";
};

void Seidel::setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){
        // Right hand side of the linear system Ax = b, in includes the boundary conditions (and the heat generation term in future implementation).
        // Receives a vector of boundary conditions where in case 
        // of cartesian coordinates, has the shape: [Lower BC, Left BC, Right BC, Upper BC]
        array<int,2> nodes = mesh->getNumNodes();
        size_t no_of_elems = nodes[0]*nodes[1];
        vector<double> bc (4,0.0);
        vector<BoundaryTypes> bc_types = pde->getBcTypes();
        __b = vector<double> (no_of_elems, 0.0);
        double heat_ss = pde->getHeatSs();
        int row;

        // for (auto &elem: bc_types){
        //     cout << static_cast<std::underlying_type<BoundaryTypes>::type>(elem)<< endl;
        // }

        if(pde->isTestCase()){
            bc[0] = pde->getBcs()[0];
            bc[1] = bc[0];
            bc[2] = bc[0];
            bc[3] = pde->getBcs()[1];
        }
        else {
            bc = pde->getBcs();
        }

    // Inner points of the mesh are 0 for vector b when no internal heat generation is present.

    // Upper boundary condition
        int j =  nodes[1]-1;
        if (bc_types[3] == BoundaryTypes::Dirichlet) {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b[row] = bc[3];
            }
        }
        else {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b[row] = bc[3]*mesh->getStepSize()[1];
            }            
        }

    // Bottom boundary condition
        j = 0;
        if (bc_types[0] == BoundaryTypes::Dirichlet) {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b[row] = bc[0];
            }
        }
        else {
            for (int i=0; i< nodes[0];i++) {
                row = i + j* nodes[0];
                __b[row] = bc[0]*mesh->getStepSize()[1];
            }
        }

    // Left boundary condition
        int i = 0;
        if (bc_types[1] == BoundaryTypes::Dirichlet) {
            for (int j=1; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b[row] = bc[1];
            }  
        }
        else {
            for (int j=0; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b[row] = bc[1]*mesh->getStepSize()[0];
            }              
        }


    // Right boundary condition
        i =  nodes[0]-1;
        if (bc_types[2] == BoundaryTypes::Dirichlet) {
            for (int j=0; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b[row] = bc[2];
            }  
        }
        else {
            for (int j=0; j< nodes[1];j++) {
                row = i + j* nodes[0];
                __b[row] = bc[2]*mesh->getStepSize()[0];
            }              
        }

        if (heat_ss!=0){ // Accounting only for point heat generation/sink
            double k = pde->getK();
            array<double,2> step = mesh->getStepSize();
            array<int,2> heat_ss_node_loc;
            heat_ss_node_loc[0] = (int) round(pde->getHeatSsLoc()[0]/step[0]);
            heat_ss_node_loc[1] = (int) round(pde->getHeatSsLoc()[1]/step[1]);
            __b[heat_ss_node_loc[0] + heat_ss_node_loc[1]*nodes[0]] += - heat_ss/k * pow(step[0]*step[1], 2);
        }

        return;
    }
    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar) {
        return;
    }
};

void Seidel::setupSolver(){
    vector<string> message;
    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();

    message = {"- Maximum number of iterations for the solver (must be integer, >0): ",
                "- Residual error tolerance to stop the solver (must be real, >0): ",
                "- Relaxation parameter omega (must be real, >1): "};

    __max_iter = utils->requestInput('i', 1, static_cast<int>(INFINITY), message[0]);
    __tol = utils->requestInput('d', 1*pow(10,-16), static_cast<double>(INFINITY), message[1]);
    __omega = utils->requestInput('d', 1.0, static_cast<double>(INFINITY), message[2]);
    return;
}

void Seidel::solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh) {
// Solve a 2D poisson's equation using the 4 point stencil finite differences method.
// Creating the finite difference system for which Mx = b
    Seidel::setupRhs(pde, mesh);
    Seidel::setupSolver();
    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();
    utils->print_vector(__b);
    
// Solution of the system: 
    bool flag = true;
    vector<bool> dirichlet_neigh (4,false);
    double res_norm;
    double tmp;
    double tmp2;
    double den = 1.0;
    int row;
    size_t iter = 0;

    array<int,2> nodes = mesh->getNumNodes();
    vector<double> tmp_sol(mesh->getTotalNodes() + 2*nodes[0], 0.0);
    array<double,2> step_inv;

    std::transform(mesh->getStepSize().begin(), mesh->getStepSize().end(),
                step_inv.begin(),
                [](double &x){return 1/pow(x,2);});

    utils->print_vector(step_inv);
    
    while (flag){
        res_norm = 0.0;
        for(int j=0; j < nodes[1];j++){
            for(int i=0; i < nodes[0];i++){
                row = i + j*nodes[0] + nodes[0];

                if(j==0 && pde->getBcTypes()[0]==BoundaryTypes::Dirichlet) {
                    tmp_sol[row] = __b[row - nodes[0]];
                    continue;
                }
                else if (i==0 && pde->getBcTypes()[1]==BoundaryTypes::Dirichlet)
                {
                    tmp_sol[row] = __b[row - nodes[0]];
                    continue;
                }
                else if (i==nodes[0]-1 && pde->getBcTypes()[2]==BoundaryTypes::Dirichlet)
                {
                    tmp_sol[row] = __b[row - nodes[0]];
                    continue;
                }                
                else if (j==nodes[1]-1 && pde->getBcTypes()[3]==BoundaryTypes::Dirichlet)
                {
                    tmp_sol[row] = __b[row - nodes[0]];
                    continue;
                }
                else {
                    den = -2*(step_inv[0] + step_inv[1]);
                    tmp2 = step_inv[0]*(tmp_sol[row-1] + tmp_sol[row+1])
                    + step_inv[1]*(tmp_sol[row-nodes[0]] + tmp_sol[row+nodes[0]]);
                }

                tmp_sol[row] = (1.0 - __omega)*tmp_sol[row] + __omega*(__b[row - nodes[0]] - tmp2)/den;
            
            }
        }

        for(int j=0; j < nodes[1];j++){
            for(int i=0; i < nodes[0];i++){
                row = i + j*nodes[0] + nodes[0];

                if(j==0 && pde->getBcTypes()[0]==BoundaryTypes::Dirichlet) {
                    continue;
                }
                else if (i==0 && pde->getBcTypes()[1]==BoundaryTypes::Dirichlet)
                {
                    continue;
                }
                else if (i==nodes[0]-1 && pde->getBcTypes()[2]==BoundaryTypes::Dirichlet)
                {
                    continue;
                }        
                else if (j==nodes[1]-1 && pde->getBcTypes()[3]==BoundaryTypes::Dirichlet)
                {
                    continue;
                }
                tmp = __b[row - nodes[0]] 
                    - (step_inv[0]*tmp_sol[row-1] + step_inv[0]*tmp_sol[row+1]
                    + step_inv[1]*tmp_sol[row-nodes[0]] + step_inv[1]*tmp_sol[row+nodes[0]]
                    + tmp_sol[row]*den);

                res_norm += pow(tmp,2);
            }
        }

        res_norm = sqrt(res_norm);
        iter += 1;
        if (res_norm < __tol || iter == __max_iter) {
            flag = false;
            cout << "res_norm = " << res_norm;
            cout << "\niter = " << iter << endl;
        }
        if (iter%10 == 0) cout << "iteration " << iter << " done! \n";

    }
    __solution[0] = vector<double>(tmp_sol.begin() + nodes[0],
                    tmp_sol.begin()+(mesh->getTotalNodes() + nodes[0]));

    return;
}