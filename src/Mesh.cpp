#include "headers/Mesh.h"


Mesh::Mesh(const std::unique_ptr<Domain>& domain) {
    // TODO 1: request input for number of nodes

    if(domain->getType() == Type::Polar) {
        double maxDim = domain->getDimensions()[0];
        _step_size = {maxDim/(_n_of_nodes[0]-1), maxDim/(_n_of_nodes[1]-1)};
        _total_nodes = ((_n_of_nodes[0]-1)*_n_of_nodes[1]) + 1;
    }
    else if (domain->getType() == Type::Cartesian) {
        vector<double> maxDims = domain->getDimensions();
        _step_size = {maxDims[0]/(_n_of_nodes[0]-1), maxDims[1]/(_n_of_nodes[1]-1)};
        _total_nodes = _n_of_nodes[0]*_n_of_nodes[1];
    }

    generateMesh(std::move(domain));
}

void Mesh::generateMesh(const std::unique_ptr<Domain>& domain) {
    if(domain->getType() == Type::Cartesian) {
        _generateCartesianMesh(domain->getDimensions());
    }
    else if (domain->getType() == Type::Polar) {
        _generatePolarMesh(domain->getDimensions());
    }
}

void Mesh::_generateCartesianMesh (const vector<double>& dims) {
    double dx = _step_size[0], dy = _step_size[1];
    int n_x = _n_of_nodes[0], n_y = _n_of_nodes[1];
    double x = dims[0], y = dims[1];

    int current_node = 0;
    for (int j = 0; j < n_y; j++) {
        for (int i = 0; i < n_x; i++) {
            
            // Specify if the node is at the boundary or not (if yes, then directly used the x or y values)
            // Specify the coordinates of the current node
            if (i == 0) {
                _mesh[current_node][0] = 1;
                _mesh[current_node][1] = 0;
                _mesh[current_node][2] = j*dy;
            }
            else if (i == n_x - 1) {
                _mesh[current_node][0] = 1;
                _mesh[current_node][1] = x;
                _mesh[current_node][2] = j*dy;
            }
            else if (j == 0) {
                _mesh[current_node][0] = 1;
                _mesh[current_node][1] = i*dx;
                _mesh[current_node][2] = 0;
            }
            else if (j == n_y - 1) {
                _mesh[current_node][0] = 1;
                _mesh[current_node][1] = i*dx;
                _mesh[current_node][2] = y;
            }
            else {
                _mesh[current_node][0] = 0;
                _mesh[current_node][1] = i*dx;
                _mesh[current_node][2] = j*dy;
            }

            current_node++;
        }
    }
}

void Mesh::_generatePolarMesh (const vector<double>& dims) {
    double dr = _step_size[0], dt = _step_size[1];
    int n_r = _n_of_nodes[0], n_t = _n_of_nodes[1];
    double r = dims[0];

    // Initialize center point
    _mesh[0][0] = 0;
    _mesh[0][1] = 0;
    _mesh[0][2] = 0;
    
    int current_node = 1;
    
    for (int j = 0; j < n_t; j++) {
        for (int i = 1; i < n_r; i++) {
            
            // Specify if the node is at the boundary or not (if yes, then directly used the x or y values)
            // Specify the coordinates of the current node
            if (i == n_r - 1) {
                _mesh[current_node][0] = 1;
                _mesh[current_node][1] = r;
                _mesh[current_node][2] = j*dt;
            }
            else {
                _mesh[current_node][0] = 0;
                _mesh[current_node][1] = i*dr;
                _mesh[current_node][2] = j*dt;
            }

            current_node++;
        }
    }
}

vector<vector<double>> Mesh::getMesh() const {return _mesh;}