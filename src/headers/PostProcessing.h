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
        vector<vector<double>> _sol;  //col 1: t_num, col 2: t_anal, col 3: error

        Mesh* _mesh;
        Initiation* _pde_type;

        vector<double> _error;

    public:
        PostProcessing(const std::unique_ptr<Initiation> &,
        const std::unique_ptr<Mesh> &,
        const std::unique_ptr<Solver> &);

        void printError();
        void plotResult();
        void exportResult();
};

#endif