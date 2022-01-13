#ifndef __ELLIPTICPDE_H_
#define __ELLIPTICPDE_H_

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <limits>

#include "Initiation.h"
#include "Utilities.h"


using std::vector;
using std::array;
using std::cout;
using std::endl;

class Domain;

class EllipticPDE : public Initiation {
    private:
        vector<BoundaryTypes> _boundary_types = {}; // Types of BCs for all possible boundaries
        vector<double> _boundary_values = {}; // Values of BCs for all possible boundaries

        double _heat_ss; // Heat source/sink
        array <double, 2> _heat_ss_location = {}; // location of Heat source/sink
        double _k; // Coefficient of thermal conductivity

    public:
        EllipticPDE(CoordinateSystem, bool test_case, bool is_homogeneous);
        
        // Setting BCs and Inhomogeneous values
        void setBCs ();
        void setInHomogeneous(const std::unique_ptr<Domain>& domain);
        
        // Getters
        vector<BoundaryTypes> getBcTypes() const;
        vector<double> getBcs() const;
        double getHeatSs() const;
        double getK() const;
        array<double,2> getHeatSsLoc() const;
};

#endif