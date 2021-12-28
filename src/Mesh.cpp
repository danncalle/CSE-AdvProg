#include "headers/Mesh.h"


Mesh::Mesh(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Domain>& domain) {
    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();
    vector<string> m = {"* Enter no. of nodes for first direction (integer, >1): ",
                        "* Enter no. of nodes for second direction (integer, >1): "};

    for(int i =0; i < 2; i++) {
        _n_of_nodes[i] = utils->requestInput('i', 2, static_cast<int>(INFINITY), m[i]);
    }

    if(pde->getCoordinateSystem() == CoordinateSystem::Polar) {
        if(domain->getShape() == Shape::Circle) {
            double maxDim = domain->getDimensions()[0];
            _step_size = {maxDim/(_n_of_nodes[0]-1), 360.0/(_n_of_nodes[1])};
            _total_nodes = ((_n_of_nodes[0]-1)*_n_of_nodes[1]) + 1;
        }
        else if(domain->getShape() == Shape::Oval) {
            double a = domain->getDimensions()[0], b = domain->getDimensions()[1];
            _step_size[1] = 360.0/_n_of_nodes[1];
            _total_nodes = ((_n_of_nodes[0]-1)*_n_of_nodes[1]) + 1;
        }
        
        _generatePolarMesh(domain->getDimensions(), domain->getShape());
    }
    else if (pde->getCoordinateSystem() == CoordinateSystem::Cartesian) {
        array<double,2> maxDims = domain->getDimensions();
        _step_size = {maxDims[0]/(_n_of_nodes[0]-1), maxDims[1]/(_n_of_nodes[1]-1)};
        _total_nodes = _n_of_nodes[0]*_n_of_nodes[1];
        
        _generateCartesianMesh(maxDims);
    }
}

void Mesh::_generateCartesianMesh (const array<double,2>& dims) {
    double dx = _step_size[0], dy = _step_size[1];
    int n_x = _n_of_nodes[0], n_y = _n_of_nodes[1];
    double x = dims[0], y = dims[1];
    
    
    for (int j = 0; j < n_y; j++) {
        for (int i = 0; i < n_x; i++) {
            
            // Specify if the node is at the boundary or not (if yes, then directly used the x or y values)
            // Specify the coordinates of the current node
            if (i == 0) _mesh.push_back({1, 0, j*dy});
            else if (i == n_x - 1) _mesh.push_back({1, x, j*dy});
            else if (j == 0) _mesh.push_back({1, i*dx, 0});
            else if (j == n_y - 1) _mesh.push_back({1, i*dx, y});
            else _mesh.push_back({0, i*dx, j*dy});
        }
    }
}

void Mesh::_generatePolarMesh (const array<double,2>& dims, Shape shape) {
    double dr = 0, dt = _step_size[1];
    int n_r = _n_of_nodes[0], n_t = _n_of_nodes[1];
    double r = 0;

    if(shape == Shape::Circle) {
        dr = _step_size[0];
        r = dims[0];
    }

    // Initialize center point
    _mesh.push_back({0, 0, 0});

    double theta, theta_rads;
    for (int j = 0; j < n_t; j++) {
        theta = j*dt;

        if(shape == Shape::Oval) {
            theta_rads = (theta*M_PI)/180.0;
            r = dims[0]*dims[1]/sqrt(pow(dims[1]*cos(theta_rads), 2)+pow(dims[0]*sin(theta_rads), 2));
            dr = r/(_n_of_nodes[0]-1);
        }

        for (int i = 1; i < n_r; i++) {
            if (i == n_r - 1) _mesh.push_back({1,r,theta});
            else _mesh.push_back({0, i*dr, theta});

        }
    }
}

vector<vector<double>> Mesh::getMesh() const {return _mesh;}