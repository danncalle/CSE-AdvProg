#ifndef __POSTPROCESSING_H_
#define __POSTPROCESSING_H_

#include <iostream>
#include <vector>
#include <memory>

#include "Utilities.h"
#include "Initiation.h"
#include "Mesh.h"
#include "Solver.h"

#include "../external/matplotlibcpp.h"

using std::vector;
using std::array;

class PostProcessing {
    private:
        // vector<vector<double>> _mesh;
        vector<vector<double>> _sol;

        Mesh* _mesh;    // raw pointer to the mesh domain
        Initiation* _pde_type; // raw pointer to the chosen pde type

        vector<double> _error; //error vector to be calculated for test case

    public:
        PostProcessing(const std::unique_ptr<Initiation> &,
        const std::unique_ptr<Mesh> &,
        const std::unique_ptr<Solver> &);

        void printError();
        void plotResult();
        void exportResult();
};

#endif