#ifndef ELLIPTICPDE.H
#define ELLIPTICPDE.H

#include <iostream>
#include <vector>
#include <array>

#include "Initiation.h"

using std::vector;
using std::array;

enum class BoundaryTypes {
    Dirichlet = 1,
    Neumann = 2
};

class EllipticPDE : public Initiation {
    private:
        vector<BoundaryTypes> _boundary_types; // Types of BCs for all possible boundries
        vector<double> _boundary_values; // Values of BCs for all possible boundries

        const double _heat_ss; // Heat source/sink
        const array <double, 2> _heat_ss_location = {}; // location of Heat source/sink
        const double _k; // Coefficient of thermal conductivity

    public:
        EllipticPDE(CoordinateSystem, bool test_case, bool is_homogeneous, PDEType PDE_type, double heat_ss, array <double, 2> location, double k);
        void setBCs ();

};

#endif