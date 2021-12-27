#ifndef MESH.H
#define MESH.H

#include <array>
#include <vector>
#include <memory>
#include <array>

#include "Initiation.h"
#include "Domain.h"

using std::array;
using std::vector;
using std::array;

class Mesh {
    private:
        array<int,2> _n_of_nodes ={0,0};
        array<double, 2> _step_size = {0,0};
        size_t _total_nodes;

        vector<vector<double>> _mesh;

        void _generateCartesianMesh(const array<double,2>&);
        void _generatePolarMesh(const array<double,2>&, Shape);

    public:
        Mesh(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Domain>& domain);
        vector<vector<double>> getMesh() const;
};

#endif