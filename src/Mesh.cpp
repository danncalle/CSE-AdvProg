#include "./headers/Mesh.h"
#include <algorithm>


Mesh::Mesh(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Domain>& domain) {
    // use some functionalities from the Utilities class
    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();
    
    vector<string> m = {"* Enter no. of nodes for first direction (integer, >1): ",
                        "* Enter no. of nodes for second direction (integer, >1): "};

    // requesting the number of nodes in each direction
    for(int i =0; i < 2; i++) {
        __n_of_nodes[i] = utils->requestInput('i', 2, static_cast<int>(INFINITY), m[i]);
    }

    // specify the step size and the total number of nodes depending on the domain coordinate system and the type of the shape
    /* == Polar case == */
    if(pde->getCoordinateSystem() == CoordinateSystem::Polar) {
        /*-- Circle --*/
        if(domain->getShape() == Shape::Circle) {
            double maxDim = domain->getDimensions()[0];
            __step_size = {maxDim/(__n_of_nodes[0]-1), 360.0/(__n_of_nodes[1])};
            __total_nodes = ((__n_of_nodes[0]-1)*__n_of_nodes[1]) + 1;
        }
        /*-- Oval --*/
        else if(domain->getShape() == Shape::Oval) {
            double a = domain->getDimensions()[0], b = domain->getDimensions()[1];
            __step_size[1] = 360.0/__n_of_nodes[1];
            __total_nodes = ((__n_of_nodes[0]-1)*__n_of_nodes[1]) + 1;
        }
        
        // invoke the polar meshing function
        _generatePolarMesh(domain->getDimensions(), domain->getShape());
    }
    /* == Cartesian case == */
    else if (pde->getCoordinateSystem() == CoordinateSystem::Cartesian) {
        array<double,2> maxDims = domain->getDimensions();
        __step_size = {maxDims[0]/(__n_of_nodes[0]-1), maxDims[1]/(__n_of_nodes[1]-1)};
        __total_nodes = __n_of_nodes[0]*__n_of_nodes[1];
        
        // invoke the cartesian meshing function
        _generateCartesianMesh(maxDims);
    }
}

void Mesh::_generateCartesianMesh (const array<double,2>& dims) {
    // @ Fill-in the __mesh matrix using the dimensions of the domain and the number of nodes in each direction
    
    // initialize the values to be used throughout from class members (for better readability of the remaining code)
    double dx = __step_size[0], dy = __step_size[1];
    int n_x = __n_of_nodes[0], n_y = __n_of_nodes[1];
    double x = dims[0], y = dims[1];
    
    
    for (int j = 0; j < n_y; j++) {
        for (int i = 0; i < n_x; i++) {   
            // Specify if the node is at the boundary or not (if yes, then directly used the x or y values)
            // Specify the coordinates of the current node
            if (i == 0) __mesh.push_back({1, 0, j*dy});
            else if (i == n_x - 1) __mesh.push_back({1, x, j*dy});
            else if (j == 0) __mesh.push_back({1, i*dx, 0});
            else if (j == n_y - 1) __mesh.push_back({1, i*dx, y});
            else __mesh.push_back({0, i*dx, j*dy});
        }
    }
}

void Mesh::_generatePolarMesh (const array<double,2>& dims, Shape shape) {
    // initialize the values to be used throughout from class members (for better readability of the remaining code)
    double dr = 0, dt = __step_size[1];
    int n_r = __n_of_nodes[0], n_t = __n_of_nodes[1];
    double r = 0;

    if(shape == Shape::Circle) {
        dr = __step_size[0];
        r = dims[0];
    }

    // Initialize center point
    __mesh.push_back({0, 0, 0});

    double theta, theta_rads;
    for (int j = 0; j < n_t; j++) {
        theta = j*dt;

        // for the oval shape, the radius changes at each theta. update the radius and the radius step size for each theta
        if(shape == Shape::Oval) {
            theta_rads = (theta*M_PI)/180.0;
            r = dims[0]*dims[1]/sqrt(pow(dims[1]*cos(theta_rads), 2)+pow(dims[0]*sin(theta_rads), 2));
            dr = r/(__n_of_nodes[0]-1);
        }

        // Specify if the node is at the boundary or not (if yes, then directly used the x or y values)
        // Specify the coordinates of the current node
        for (int i = 1; i < n_r; i++) {
            if (i == n_r - 1) __mesh.push_back({1,r,theta});
            else __mesh.push_back({0, i*dr, theta});

        }
    }
}

// Getters
vector<vector<double>> Mesh::getMesh() const {return __mesh;}
array<int,2> Mesh::getNumNodes() const {return __n_of_nodes;}
array<double,2> Mesh::getStepSize() const {return __step_size;}
size_t Mesh::getTotalNodes() const {return __total_nodes;}
