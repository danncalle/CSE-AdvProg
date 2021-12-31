#ifndef __MESH_H_
#define __MESH_H_

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
    protected:
        array<int,2> __n_of_nodes ={0,0};
        array<double, 2> __step_size = {0,0};
        size_t __total_nodes;

        vector<vector<double>> __mesh;

        void _generateCartesianMesh(const array<double,2>&);
        void _generatePolarMesh(const array<double,2>&, Shape);

    public:
        Mesh(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Domain>& domain);
        vector<vector<double>> getMesh() const;
        array<int,2> getNumNodes() const;
        array<double,2> getStepSize() const;
        size_t getTotalNodes() const;
};

#endif