#ifndef MESH.H
#define MESH.H

#include <array>
#include <vector>
#include <memory>
#include <array>

#include "Domain.h"

using std::array;
using std::vector;
using std::array;

class Mesh {
    private:
        array<int,2> _n_of_nodes;
        array<double, 2> _step_size;
        size_t _total_nodes;

        vector<vector<double>> _mesh;

        void _generateCartesianMesh(const vector<double>&);
        void _generatePolarMesh(const vector<double>&);

    public:
        Mesh(const std::unique_ptr<Domain>& domain);
        void generateMesh(const std::unique_ptr<Domain>& domain);
        vector<vector<double>> getMesh() const;
};

#endif