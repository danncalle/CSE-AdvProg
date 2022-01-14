#include "./headers/Solver.h"

/* Default constructor for a general Solver object.
*/
Solver::Solver(SolverType solver_type, bool test_case):
__solver_type(solver_type), __test_case(test_case) {
    cout << "Solver initiated correctly\n";
};

/* Get method for error parameter of the solver class. Only takes values when test_case is true.
*/
vector<double> Solver::getError() const {return __error;};

/* Get method for the solution vector of the solver class. When test_case is true, it contains both the
* numerical and analytical solutions over the meshed domain.
*/
vector<vector<double>> Solver::getSolution() const {return __solution;};

/* Method that solves analytical case over the meshed domain. Is only used when test_case is true.
*
*/
void Solver::solve_analytical(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain){

// Implementation for a test case in cartesian coordinates. Refer to ReadMe file for more information.
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){

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

// Placeholder for test case in polar coordinates.
    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar){
        return;
    }

};

/* Set method to calculate the relative error of the numerical solution compared to the analytical solution. Is empty when test_case is false.
*/
void Solver::setError(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain){
    if(!__test_case){
        __error = {};
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

        return;
    }
}

/* Default constructor for the LU_direct class. 
*/
LU_direct::LU_direct(const bool test_case):
Solver(SolverType::LU_dir, test_case) {
    cout << "LU_direct solver initiated correctly\n";
}

/* Method to set up the right-hand side of the Finite Differences linear system Ax = b. This methods 
* sets the b. A 5-point stencil is used for the system. Current implementation accounts only for
* squared or rectangular shapes in cartesian coordinates. Therefore, only 4 boundaries are considered.
*/
void LU_direct::setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){

// Implementation for the cartesian coordinate case. 
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){

        int row;
        int t;
        std::list<array<int,2>> bc_order_idx {};

        double heat_ss = pde->getHeatSs();
        array<int,2> nodes = mesh->getNumNodes();
        size_t no_of_elems = nodes[0]*nodes[1];
        __b = vector<double> (no_of_elems, 0.0);

        vector<double> bc (4,0.0);
        vector<BoundaryTypes> bc_types = pde->getBcTypes();

// Set the order of implementation for boundary conditions. When two connecting lines of the
// domain have one Neumann and one Dirichlet conditions respectively, the Dirichlet has 
// predominance in the shared node. When there are two continuous Dirichlet conditions, the 
// order is not relevant, because the corner node is not considered for the 5-point stencil.
        if(pde->isTestCase()){
            bc[0] = pde->getBcs()[0];
            bc[1] = bc[0];
            bc[2] = bc[0];
            bc[3] = pde->getBcs()[1];
            bc_order_idx = {{3,nodes[1]-1},{0,0},{1,0},{2,nodes[0]-1}};
        }
        else {
            bc = pde->getBcs();
            for (int i=0;i<bc.size(); i++) {

                if (i==0 || i==1) t = 0;
                else if (i==2) t = nodes[0]-1;
                else if (i==3) t = nodes[1]-1;

                if (bc_types[i] == BoundaryTypes::Dirichlet){
                    bc_order_idx.push_back({i,t});
                }
                else if (bc_types[i] == BoundaryTypes::Neumann){
                    bc_order_idx.push_front({i,t});
                }
            }
        }

// Set the values for the vector b according to the respective boundary condition in 
// the predefined order.

        for (auto &elem: bc_order_idx) {
            if(elem[0] == 0 || elem[0] == 3){
                int j =  elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int i=0; i< nodes[0];i++) {
                        row = i + j* nodes[0];
                        __b[row] = bc[elem[0]];
                    }
                }
                else {
                    for (int i=0; i< nodes[0];i++) {
                        row = i + j* nodes[0];
                        __b[row] = 2*bc[elem[0]]*mesh->getStepSize()[1]*pow(mesh->getStepSize()[0],2);
                    }            
                }                
            }
            else if(elem[0] == 1 || elem[0] == 2) {
                int i = elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int j=0; j< nodes[1];j++) {
                        row = i + j* nodes[0];
                        __b[row] = bc[elem[0]];
                    }  
                }
                else {
                    for (int j=0; j< nodes[1];j++) {
                        row = i + j* nodes[0];
                        __b[row] = 2*bc[elem[0]]*mesh->getStepSize()[0]*pow(mesh->getStepSize()[0],2);
                    }              
                }                
            }
        }

// Accounting only for point heat generation/sink. Internal nodes with no heat generation/sink or
// Neumann boundary have a zero entry in the b vector. A positive value input by the user accounts for 
// heat extraction of the system. 
        if (!pde->isHomogeneous()){ 
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

/* Method to set up the left-hand side of the Finite Differences linear system Ax = b. This methods 
* sets the matrix A. A 5-point stencil is used for the system. Current implementation accounts only for
* squared or rectangular shapes in cartesian coordinates. Therefore, only 4 boundaries are considered and 
* the matrix has at maximum 5 non-zero entries per row. It is a band matrix with width 2*N_x.
*/
void LU_direct::setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){

// Implementation for the cartesian coordinate case. 
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){

        int row;
        int t;

        array<int,2> nodes = mesh->getNumNodes();
        size_t no_of_elems = nodes[0]*nodes[1];
        vector<BoundaryTypes> bc_types = pde->getBcTypes();
        std::list<array<int,2>> bc_order_idx {};

        vector<double> rows (no_of_elems, 0.0);
        __M = vector<vector<double>> (no_of_elems,rows);
        __L = vector<vector<double>> (no_of_elems,rows);
        __U = vector<vector<double>> (no_of_elems,rows);


        double dx = mesh->getStepSize()[0];
        double dy = mesh->getStepSize()[1];

// Set the order of implementation for boundary conditions. When two connecting lines of the
// domain have one Neumann and one Dirichlet conditions respectively, the Dirichlet has 
// predominance in the shared node. When there are two continuous Dirichlet conditions, the 
// order is not relevant, because the corner node is not considered for the 5-point stencil.
        if(pde->isTestCase()){
            bc_order_idx = {{3,nodes[1]-1},{0,0},{1,0},{2,nodes[0]-1}};

        }
        else {
            for (int i=0;i<pde->getBcs().size(); i++) {

                if (i==0 || i==1) t = 0;
                else if (i==2) t = nodes[0]-1;
                else if (i==3) t = nodes[1]-1;

                if (bc_types[i] == BoundaryTypes::Dirichlet){
                    bc_order_idx.push_back({i,t});
                }
                else if (bc_types[i] == BoundaryTypes::Neumann){
                    bc_order_idx.push_front({i,t});
                }
            }
        }

// Filling the elements of the matrix
// Rows corresponding to inner points of the mesh (not on boundaries)
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

// Rows corresponding to boundary nodes according to the respective boundary condition in 
// the predefined order.
        for (auto &elem: bc_order_idx) {

//Lower boundary 
            if(elem[0] == 0){
                int j =  elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int i=0; i< nodes[0];i++) {
                        row = i + j* nodes[0];
                        __M[row] = vector<double> (no_of_elems,0.0);
                        __M[row][row] = 1;
                    }
                }
                else {
                    for (int i=1; i< nodes[0]-1;i++) {
                        row = i + j* nodes[0];
                        __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M[row][row-1] = pow(dy, 2);
                        __M[row][row+1] = pow(dy, 2);
                        __M[row][row+nodes[0]] = 2 * pow(dx, 2);
                    }

// Lower Right corner
                    if (bc_types[2] == BoundaryTypes::Neumann){
                        __M[nodes[0]-1][nodes[0]-1] = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M[nodes[0]-1][nodes[0]-2] = 2 * pow(dy, 2);
                        __M[nodes[0]-1][2*nodes[0]-1] = 2 * pow(dx, 2);                        
                    }
                }                
            }

// Upper boundary
            if(elem[0] == 3){
                int j =  elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int i=0; i< nodes[0];i++) {
                        row = i + j* nodes[0];
                        __M[row] = vector<double> (no_of_elems,0.0);
                        __M[row][row] = 1;
                    }
                }
                else {
                    for (int i=1; i< nodes[0]-1;i++) {
                        row = i + j* nodes[0];
                        __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M[row][row-1] = pow(dy, 2);
                        __M[row][row+1] = pow(dy, 2);
                        __M[row][row-nodes[0]] = 2 * pow(dx, 2);
                    }    

// Upper Left corner
                    if (bc_types[1] == BoundaryTypes::Neumann){
                        __M[(nodes[1]-1)*nodes[0]][(nodes[1]-1)*nodes[0]] = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M[(nodes[1]-1)*nodes[0]][(nodes[1]-1)*nodes[0]+1] = 2 * pow(dy, 2);
                        __M[(nodes[1]-1)*nodes[0]][nodes[1]*nodes[0] - 2*nodes[0]] = 2 * pow(dx, 2);                       
                    }
                }                
            }

// Left boundary
            else if(elem[0] == 1) {
                int i = elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int j=0; j< nodes[1];j++) {
                        row = i + j* nodes[0];
                        __M[row] = vector<double> (no_of_elems,0.0);
                        __M[row][row] = 1;
                    }  
                }
                else {
                    for (int j=1; j< nodes[1]-1;j++) {
                        row = i + j* nodes[0];
                        __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M[row][row-1] = pow(dy, 2);
                        __M[row][row+1] = pow(dy, 2);
                        __M[row][row-nodes[0]] = 2 * pow(dx, 2);
                    }     

// Lower Left corner
                    if (bc_types[0] == BoundaryTypes::Neumann){
                        __M[0][0] = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M[0][1] = 2 * pow(dy, 2);
                        __M[0][nodes[0]] = 2 * pow(dx, 2);
                    }         
                }                
            }

// Right boundary
            else if(elem[0] == 2) {
                int i = elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int j=0; j< nodes[1];j++) {
                        row = i + j* nodes[0];
                        __M[row] = vector<double> (no_of_elems,0.0);
                        __M[row][row] = 1;
                    }  
                }
                else {
                    for (int j=1; j< nodes[1]-1;j++) {
                        row = i + j* nodes[0];
                        __M[row][row] = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M[row][row-1] = pow(dy, 2);
                        __M[row][row+1] = pow(dy, 2);
                        __M[row][row-nodes[0]] = 2 * pow(dx, 2);
                    }   

// Upper Right corner
                    if (bc_types[3] == BoundaryTypes::Neumann){
                        __M[nodes[0]*nodes[1] - 1][nodes[0]*nodes[1] - 1] = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M[nodes[0]*nodes[1] - 1][nodes[0]*nodes[1] - 2] = 2 * pow(dy, 2);
                        __M[row][nodes[0]*(nodes[1] - 1) - 1] = 2 * pow(dx, 2);
                    }         
                }                
            }
        }

        return;
    }

    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar) {
        return;
    }   
};

/* Perform classical LU factorization on the matrix of the Finite Differences system.
* This approach takes O(N^3) operations to be completed since it considers the full matrix.
*/
void LU_direct::LUFactorization() {

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

/* Solve an upper-triangular linear system via backward substitution.
*/
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

/* Solve a lower-triangular linear system via forward substitution.
*/
vector<double> LU_direct::forwardSubstitution(const vector<double>& rhs){

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

/* Solve the Finite Differences linear system by first setting the system itself, i.e.,
* A and b in Ax=b. Then performs the LU factorization and solves the system of equations via 
* forward and backward substitution.
*/
void LU_direct::solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh) {

    LU_direct::setupMatrix(pde, mesh);
    LU_direct::setupRhs(pde, mesh);

    LUFactorization();
    vector<double> y = forwardSubstitution(__b);
    __solution[0] = backwardSubstitution(y);

    cout << "Solution finished" << endl;

    return;
}

/* Method to allow the user to create a file with the L and U matrices from the
* LU factorization in a csv format.
*/
void LU_direct::writeLU() {

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

/* Method to allow the user to load a file with the L and U matrices from the
* LU factorization from a csv format. It then assigns these values to the L and U matrices of 
* the solver object.
*/
void LU_direct::loadLU() {

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

// Check if resultfile connection is ready, otherwise stop this function 
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

/* Default constructor for the LU_sparse class. This class uses the Eigen library methods and data
* structures to solve the PDE.
*/
LU_sparse::LU_sparse(bool test_case):
Solver(SolverType::LU_sp, test_case){
    cout << "LU_sparse solver initiated correctly\n";
};

/* Method to set up the right-hand side of the Finite Differences linear system Ax = b. This methods 
* sets the b. A 5-point stencil is used for the system. Current implementation accounts only for
* squared or rectangular shapes in cartesian coordinates. Therefore, only 4 boundaries are considered.
*/
void LU_sparse::setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){

// Implementation for the cartesian coordinate case. 
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){

        int row;
        int t;
        std::list<array<int,2>> bc_order_idx {};

        double heat_ss = pde->getHeatSs();
        array<int,2> nodes = mesh->getNumNodes();
        size_t no_of_elems = nodes[0]*nodes[1];
        __b = VectorXd (no_of_elems);

        vector<double> bc (4,0.0);
        vector<BoundaryTypes> bc_types = pde->getBcTypes();

// Set the order of implementation for boundary conditions. When two connecting lines of the
// domain have one Neumann and one Dirichlet conditions respectively, the Dirichlet has 
// predominance in the shared node. When there are two continuous Dirichlet conditions, the 
// order is not relevant, because the corner node is not considered for the 5-point stencil.
        if(pde->isTestCase()){
            bc[0] = pde->getBcs()[0];
            bc[1] = bc[0];
            bc[2] = bc[0];
            bc[3] = pde->getBcs()[1];
            bc_order_idx = {{3,nodes[1]-1},{0,0},{1,0},{2,nodes[0]-1}};
        }
        else {
            bc = pde->getBcs();
            for (int i=0;i<bc.size(); i++) {

                if (i==0 || i==1) t = 0;
                else if (i==2) t = nodes[0]-1;
                else if (i==3) t = nodes[1]-1;

                if (bc_types[i] == BoundaryTypes::Dirichlet){
                    bc_order_idx.push_back({i,t});
                }
                else if (bc_types[i] == BoundaryTypes::Neumann){
                    bc_order_idx.push_front({i,t});
                }
            }
        }

// Set the values for the vector b according to the respective boundary condition in 
// the predefined order.
        for (auto &elem: bc_order_idx) {
            if(elem[0] == 0 || elem[0] == 3){
                int j =  elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int i=0; i< nodes[0];i++) {
                        row = i + j* nodes[0];
                        __b.coeffRef(row) = bc[elem[0]];
                    }
                }
                else {
                    for (int i=0; i< nodes[0];i++) {
                        row = i + j* nodes[0];
                        __b.coeffRef(row) = 2*bc[elem[0]]*mesh->getStepSize()[1]*pow(mesh->getStepSize()[0],2);
                    }            
                }                
            }
            else if(elem[0] == 1 || elem[0] == 2) {
                int i = elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int j=0; j< nodes[1];j++) {
                        row = i + j* nodes[0];
                        __b.coeffRef(row) = bc[elem[0]];
                    }  
                }
                else {
                    for (int j=0; j< nodes[1];j++) {
                        row = i + j* nodes[0];
                        __b.coeffRef(row) = 2*bc[elem[0]]*mesh->getStepSize()[0]*pow(mesh->getStepSize()[0],2);
                    }              
                }                
            }
        }

// Accounting only for point heat generation/sink. Internal nodes with no heat generation/sink or
// Neumann boundary have a zero entry in the b vector. A positive value input by the user accounts for 
// heat extraction of the system.
        if (!pde->isHomogeneous()){ 
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

/* Method to set up the left-hand side of the Finite Differences linear system Ax = b. This methods 
* sets the matrix A. A 5-point stencil is used for the system. Current implementation accounts only for
* squared or rectangular shapes in cartesian coordinates. Therefore, only 4 boundaries are considered and 
* the matrix has at maximum 5 non-zero entries per row. It is a band matrix with width 2*N_x.
*/
void LU_sparse::setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){

// Implementation for the cartesian coordinate case. 
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){


        int row;
        int t;

        array<int,2> nodes = mesh->getNumNodes();
        size_t no_of_elems = nodes[0]*nodes[1];
        vector<BoundaryTypes> bc_types = pde->getBcTypes();
        std::list<array<int,2>> bc_order_idx {};

        __M = SparseMatrix<double> (no_of_elems,no_of_elems);

        double dx = mesh->getStepSize()[0];
        double dy = mesh->getStepSize()[1];

// Set the order of implementation for boundary conditions. When two connecting lines of the
// domain have one Neumann and one Dirichlet conditions respectively, the Dirichlet has 
// predominance in the shared node. When there are two continuous Dirichlet conditions, the 
// order is not relevant, because the corner node is not considered for the 5-point stencil.
        if(pde->isTestCase()){
            bc_order_idx = {{3,nodes[1]-1},{0,0},{1,0},{2,nodes[0]-1}};

        }
        else {
            for (int i=0;i<pde->getBcs().size(); i++) {

                if (i==0 || i==1) t = 0;
                else if (i==2) t = nodes[0]-1;
                else if (i==3) t = nodes[1]-1;

                if (bc_types[i] == BoundaryTypes::Dirichlet){
                    bc_order_idx.push_back({i,t});
                }
                else if (bc_types[i] == BoundaryTypes::Neumann){
                    bc_order_idx.push_front({i,t});
                }
            }
        }


// Filling the elements of the matrix
// Rows corresponding to inner points of the mesh (not on boundaries)
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

// Rows corresponding to boundary nodes according to the respective boundary condition in 
// the predefined order.
        for (auto &elem: bc_order_idx) {

// Lower boundary 
            if(elem[0] == 0){
                int j =  elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int i=0; i< nodes[0];i++) {
                        row = i + j* nodes[0];
                        __M.row(row) *= 0.0;
                        __M.coeffRef(row, row) = 1;
                    }
                }
                else {
                    for (int i=1; i< nodes[0]-1;i++) {
                        row = i + j* nodes[0];
                        __M.coeffRef(row, row) = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M.coeffRef(row, row-1) = pow(dy, 2);
                        __M.coeffRef(row, row+1) = pow(dy, 2);
                        __M.coeffRef(row, row+nodes[0]) = 2 * pow(dx, 2);
                    }

// Lower Right corner
                    if (bc_types[2] == BoundaryTypes::Neumann){
                        __M.coeffRef(nodes[0]-1, nodes[0]-1) = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M.coeffRef(nodes[0]-1, nodes[0]-2) = 2 * pow(dy, 2);
                        __M.coeffRef(nodes[0]-1, 2*nodes[0]-1) = 2 * pow(dx, 2);                        
                    }
                }                
            }

// Upper boundary
            if(elem[0] == 3){
                int j =  elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int i=0; i< nodes[0];i++) {
                        row = i + j* nodes[0];
                        __M.row(row) *= 0.0;
                        __M.coeffRef(row, row) = 1;
                    }
                }
                else {
                    for (int i=1; i< nodes[0]-1;i++) {
                        row = i + j* nodes[0];
                        __M.coeffRef(row, row) = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M.coeffRef(row, row-1) = pow(dy, 2);
                        __M.coeffRef(row, row+1) = pow(dy, 2);
                        __M.coeffRef(row, row-nodes[0]) = 2 * pow(dx, 2);
                    }    

// Upper Left corner
                    if (bc_types[1] == BoundaryTypes::Neumann){
                        __M.coeffRef((nodes[1]-1)*nodes[0], (nodes[1]-1)*nodes[0]) = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M.coeffRef((nodes[1]-1)*nodes[0], (nodes[1]-1)*nodes[0]+1) = 2 * pow(dy, 2);
                        __M.coeffRef((nodes[1]-1)*nodes[0], nodes[1]*nodes[0] - 2*nodes[0]) = 2 * pow(dx, 2);                       
                    }
                }                
            }

// Left boundary
            else if(elem[0] == 1) {
                int i = elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int j=0; j< nodes[1];j++) {
                        row = i + j* nodes[0];
                        __M.row(row) *= 0.0;
                        __M.coeffRef(row, row) = 1;
                    }  
                }
                else {
                    for (int j=1; j< nodes[1]-1;j++) {
                        row = i + j* nodes[0];
                        __M.coeffRef(row, row) = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M.coeffRef(row, row-1) = pow(dy, 2);
                        __M.coeffRef(row, row+1) = pow(dy, 2);
                        __M.coeffRef(row, row-nodes[0]) = 2 * pow(dx, 2);
                    }   
  
// Lower Left corner
                    if (bc_types[0] == BoundaryTypes::Neumann){
                        __M.coeffRef(0, 0) = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M.coeffRef(0, 1) = 2 * pow(dy, 2);
                        __M.coeffRef(0, nodes[0]) = 2 * pow(dx, 2);
                    }         
                }                
            }

// Right boundary
            else if(elem[0] == 2) {
                int i = elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int j=0; j< nodes[1];j++) {
                        row = i + j* nodes[0];
                        __M.row(row) *= 0.0;
                        __M.coeffRef(row, row) = 1;
                    }  
                }
                else {
                    for (int j=1; j< nodes[1]-1;j++) {
                        row = i + j* nodes[0];
                        __M.coeffRef(row, row) = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M.coeffRef(row, row-1) = pow(dy, 2);
                        __M.coeffRef(row, row+1) = pow(dy, 2);
                        __M.coeffRef(row, row-nodes[0]) = 2 * pow(dx, 2);
                    }     

// Upper Right corner
                    if (bc_types[3] == BoundaryTypes::Neumann){
                        __M.coeffRef(nodes[0]*nodes[1] - 1, nodes[0]*nodes[1] - 1) = -2 * (pow(dx, 2) + pow(dy, 2));
                        __M.coeffRef(nodes[0]*nodes[1] - 1, nodes[0]*nodes[1] - 2) = 2 * pow(dy, 2);
                        __M.coeffRef(row, nodes[0]*(nodes[1] - 1) - 1) = 2 * pow(dx, 2);
                    }         
                }                
            }
        }

        return;
    }

    else if(pde->getCoordinateSystem() == CoordinateSystem::Polar) {
        return;
    }   
};

/* Solve the Finite Differences linear system by first setting the system itself, i.e.,
* A and b in Ax=b. Then performs the sparse LU factorization method and solves the system of equations. 
* According to Eigen references, the system is solved in O(N^3/2) operations.
*/
void LU_sparse::solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){

    LU_sparse::setupMatrix(pde, mesh);
    LU_sparse::setupRhs(pde, mesh);

    Eigen::SparseLU<SparseMatrix<double>> solver;
    solver.compute(__M);

    if (solver.info()!=Eigen::ComputationInfo::Success){
        cout << "Decomposition failed" << endl;
    }

    VectorXd sol = solver.solve(__b);

    if (solver.info()!=Eigen::ComputationInfo::Success){
        cout << "Solution failed" << endl;
    }


// Converting Eigen::vector to std::vector 
    __solution[0] = vector<double> (sol.size(),0.0);
    for (int i=0; i<sol.size(); i++){
        __solution[0][i] = sol(i);
    }
    cout << "Solution finished" << endl;

    return;    
};

/* Default constructor for the Gauss_Seidel class. This class uses the Eigen library methods and data
* structures to solve the PDE.
*/
Seidel::Seidel(bool test_case):
Solver(SolverType::Gauss_Seidel, test_case){
    cout << "Gauss-Seidel solver initiated correctly\n";
};

/* Method to set up the right-hand side of the Finite Differences linear system Ax = b. This methods 
* sets the b. A 5-point stencil is used for the system. Current implementation accounts only for
* squared or rectangular shapes in cartesian coordinates. Therefore, only 4 boundaries are considered.
*/
void Seidel::setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh){

// Implementation for the cartesian coordinate case. 
    if(pde->getCoordinateSystem() == CoordinateSystem::Cartesian){

        int row;
        int t;
        std::list<array<int,2>> bc_order_idx {};

        double heat_ss = pde->getHeatSs();
        array<int,2> nodes = mesh->getNumNodes();
        size_t no_of_elems = nodes[0]*nodes[1];
        __b = vector<double> (no_of_elems, 0.0);

        vector<double> bc (4,0.0);
        vector<BoundaryTypes> bc_types = pde->getBcTypes();

// Set the order of implementation for boundary conditions. When two connecting lines of the
// domain have one Neumann and one Dirichlet conditions respectively, the Dirichlet has 
// predominance in the shared node. When there are two continuous Dirichlet conditions, the 
// order is not relevant, because the corner node is not considered for the 5-point stencil.
        if(pde->isTestCase()){
            bc[0] = pde->getBcs()[0];
            bc[1] = bc[0];
            bc[2] = bc[0];
            bc[3] = pde->getBcs()[1];
            bc_order_idx = {{3,nodes[1]-1},{0,0},{1,0},{2,nodes[0]-1}};
        }
        else {
            bc = pde->getBcs();
            for (int i=0;i<bc.size(); i++) {

                if (i==0 || i==1) t = 0;
                else if (i==2) t = nodes[0]-1;
                else if (i==3) t = nodes[1]-1;

                if (bc_types[i] == BoundaryTypes::Dirichlet){
                    bc_order_idx.push_back({i,t});
                }
                else if (bc_types[i] == BoundaryTypes::Neumann){
                    bc_order_idx.push_front({i,t});
                }
            }
        }

// Set the values for the vector b according to the respective boundary condition in 
// the predefined order.
        for (auto &elem: bc_order_idx) {
            if(elem[0] == 0 || elem[0] == 3){
                int j =  elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int i=0; i< nodes[0];i++) {
                        row = i + j* nodes[0];
                        __b[row] = bc[elem[0]];
                    }
                }
                else {
                    for (int i=0; i< nodes[0];i++) {
                        row = i + j* nodes[0];
                        __b[row] = 2*bc[elem[0]]*mesh->getStepSize()[1]*pow(mesh->getStepSize()[0],2);
                    }            
                }                
            }
            else if(elem[0] == 1 || elem[0] == 2) {
                int i = elem[1];
                if (bc_types[elem[0]] == BoundaryTypes::Dirichlet) {
                    for (int j=0; j< nodes[1];j++) {
                        row = i + j* nodes[0];
                        __b[row] = bc[elem[0]];
                    }  
                }
                else {
                    for (int j=0; j< nodes[1];j++) {
                        row = i + j* nodes[0];
                        __b[row] = 2*bc[elem[0]]*mesh->getStepSize()[0]*pow(mesh->getStepSize()[0],2);
                    }              
                }                
            }
        }

// Accounting only for point heat generation/sink. Internal nodes with no heat generation/sink or
// Neumann boundary have a zero entry in the b vector. A positive value input by the user accounts for 
// heat extraction of the system. 
        if (!pde->isHomogeneous()){ 
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

/* Sets up the parameters required for the solver, i.e.:
* Maximum number of iterations, Tolerance parameter to stop the solver and relaxation parameter.
* When the relaxation parameter is greater than 1, the SOR method is used. When it is 1, traditional
* Gauss-Seidel method is used.
*/
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

/* Solve the Finite Differences linear system by first setting the system itself, i.e.,
* A and b in Ax=b. Then performs the sparse LU factorization method and solves the system of equations. 
* According to Eigen references, the system is solved in O(N^3/2) operations.
*/
void Seidel::solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh) {

    Seidel::setupRhs(pde, mesh);
    Seidel::setupSolver();
    
// Initialize required variables
    bool flag = true;

    double res_norm;
    double tmp;
    double tmp2;
    int row;

    size_t iter = 0;

    array<int,2> nodes = mesh->getNumNodes();
    array<double,2> step_sq;

    std::transform(mesh->getStepSize().begin(), mesh->getStepSize().end(),
                step_sq.begin(),
                [](double &x){return pow(x,2);});

    double den = -2*(step_sq[0] + step_sq[1]);

    vector<double> tmp_sol(mesh->getTotalNodes() + 2*nodes[0], 0.0);

// Iterate until maximum iterations are reached or residual norm is lower than the user input tolerance.
    while (flag){

// Initialize the residual norm.
        res_norm = 0.0;
        
// Loop over every node of the mesh. The loop is performed in order of the nodes.
        for(int j=0; j < nodes[1];j++){
            for(int i=0; i < nodes[0];i++){
                row = i + j*nodes[0] + nodes[0];

// Skip the boundary nodes that have a dirichlet condition because their temperature is known.
                if((j==0 && pde->getBcTypes()[0]==BoundaryTypes::Dirichlet) ||
                (i==nodes[0]-1 && pde->getBcTypes()[2]==BoundaryTypes::Dirichlet)) {
                    tmp_sol[row] = __b[row - nodes[0]];
                    continue;
                }
                else if ((i==0 && pde->getBcTypes()[1]==BoundaryTypes::Dirichlet) ||
                (j==0 && pde->getBcTypes()[0]==BoundaryTypes::Dirichlet)) {
                    tmp_sol[row] = __b[row - nodes[0]];
                    continue;
                }
                else if ((i==nodes[0]-1 && pde->getBcTypes()[2]==BoundaryTypes::Dirichlet) ||
                (j==nodes[1]-1 && pde->getBcTypes()[3]==BoundaryTypes::Dirichlet)) {
                    tmp_sol[row] = __b[row - nodes[0]];
                    continue;
                }                
                else if ((j==nodes[1]-1 && pde->getBcTypes()[3]==BoundaryTypes::Dirichlet) ||
                (i==0 && pde->getBcTypes()[1]==BoundaryTypes::Dirichlet)) {
                    tmp_sol[row] = __b[row - nodes[0]];
                    continue;
                }

// Calculate the stencil for the current nodes. The stencil changes if the node belongs
// to a Neumann boundary and if it is an internal node.
                if (j==0 && pde->getBcTypes()[0]==BoundaryTypes::Neumann){
                    tmp2 = step_sq[1]*(tmp_sol[row-1] + tmp_sol[row+1])
                    + step_sq[0]*(2*tmp_sol[row+nodes[0]]);
                }
                else if (i==0 && pde->getBcTypes()[1]==BoundaryTypes::Neumann){
                    tmp2 = step_sq[1]*(2*tmp_sol[row+1])
                    + step_sq[0]*(tmp_sol[row-nodes[0]] + tmp_sol[row+nodes[0]]);
                }
                else if (i==nodes[0]-1 && pde->getBcTypes()[2]==BoundaryTypes::Neumann){
                    tmp2 = step_sq[1]*(2*tmp_sol[row-1])
                    + step_sq[0]*(tmp_sol[row-nodes[0]] + tmp_sol[row+nodes[0]]);
                }
                else if (j==nodes[1]-1 && pde->getBcTypes()[3]==BoundaryTypes::Neumann){
                    tmp2 = step_sq[1]*(tmp_sol[row-1] + tmp_sol[row+1])
                    + step_sq[0]*(2*tmp_sol[row-nodes[0]]);
                }
                else {
                    tmp2 = step_sq[1]*(tmp_sol[row-1] + tmp_sol[row+1])
                    + step_sq[0]*(tmp_sol[row-nodes[0]] + tmp_sol[row+nodes[0]]);
                }

// Use the stencil to update the solution.
                tmp_sol[row] = (1.0 - __omega)*tmp_sol[row] + __omega*(__b[row - nodes[0]] - tmp2)/den;
            
            }
        }

// Loop again over the nodes to calculate the residual norm of one loop of the method.
        for(int j=0; j < nodes[1];j++){
            for(int i=0; i < nodes[0];i++){
                row = i + j*nodes[0] + nodes[0];

// Skip the boundary nodes that have a dirichlet condition because their temperature is known.
                if((j==0 && pde->getBcTypes()[0]==BoundaryTypes::Dirichlet) ||
                (i==nodes[0]-1 && pde->getBcTypes()[2]==BoundaryTypes::Dirichlet)) {
                    continue;
                }
                else if ((i==0 && pde->getBcTypes()[1]==BoundaryTypes::Dirichlet) ||
                (j==0 && pde->getBcTypes()[0]==BoundaryTypes::Dirichlet)) {
                    continue;
                }
                else if ((i==nodes[0]-1 && pde->getBcTypes()[2]==BoundaryTypes::Dirichlet) ||
                (j==nodes[1]-1 && pde->getBcTypes()[3]==BoundaryTypes::Dirichlet)) {
                    continue;
                }                
                else if ((j==nodes[1]-1 && pde->getBcTypes()[3]==BoundaryTypes::Dirichlet) ||
                (i==0 && pde->getBcTypes()[1]==BoundaryTypes::Dirichlet)) {
                    continue;
                }

// Calculate the stencil for the current nodes. The stencil changes if the node belongs
// to a Neumann boundary and if it is an internal node.
                if (j==0 && pde->getBcTypes()[0]==BoundaryTypes::Neumann){
                    tmp2 = step_sq[1]*(tmp_sol[row-1] + tmp_sol[row+1])
                    + step_sq[0]*(2*tmp_sol[row+nodes[0]]);
                }
                else if (i==0 && pde->getBcTypes()[1]==BoundaryTypes::Neumann){
                    tmp2 = step_sq[1]*(2*tmp_sol[row+1])
                    + step_sq[0]*(tmp_sol[row-nodes[0]] + tmp_sol[row+nodes[0]]);
                }
                else if (i==nodes[0]-1 && pde->getBcTypes()[2]==BoundaryTypes::Neumann){
                    tmp2 = step_sq[1]*(2*tmp_sol[row-1])
                    + step_sq[0]*(tmp_sol[row-nodes[0]] + tmp_sol[row+nodes[0]]);
                }
                else if (j==nodes[1]-1 && pde->getBcTypes()[3]==BoundaryTypes::Neumann){
                    tmp2 = step_sq[1]*(tmp_sol[row-1] + tmp_sol[row+1])
                    + step_sq[0]*(2*tmp_sol[row-nodes[0]]);
                }
                else {
                    tmp2 = step_sq[1]*(tmp_sol[row-1] + tmp_sol[row+1])
                    + step_sq[0]*(tmp_sol[row-nodes[0]] + tmp_sol[row+nodes[0]]);
                }

// Update the residual norm using the stencil, the current value of the node and the right hand side of the system.
                tmp = __b[row - nodes[0]] - (tmp2 + tmp_sol[row]*den);
                res_norm += pow(tmp,2);
            }
        }

// Calculate residual norm, update iterations and validate if method must stop.
        res_norm = sqrt(res_norm);
        iter += 1;
        if (res_norm < __tol || iter == __max_iter) {
            flag = false;
            cout << "res_norm = " << res_norm;
            cout << "\niter = " << iter << endl;
        }

// Print current state of the residual norm and iterations.
        if (iter%100 == 0) {
            cout << "iteration " << iter << " done! \n";
            cout << "res_norm " << res_norm << "\n";
        }

    }

// Assign the finished solution to the solution parameter of the class.
    __solution[0] = vector<double>(tmp_sol.begin() + nodes[0],
                    tmp_sol.begin()+(mesh->getTotalNodes() + nodes[0]));

    return;
}