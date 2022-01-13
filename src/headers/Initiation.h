#ifndef __INITIATION_H_
#define __INITIATION_H_

#include <iostream>
#include <memory>
#include <vector>
#include <array>

using std::array;
using std::vector;

// enum for the type of the chosen Coordinate system
enum class CoordinateSystem {
    Cartesian = 1,
    Polar = 2
};

// enum for the type of the PDE to solve (default is 1 for this 2D heat transfer solver)
enum class PDEType {
    Elliptic = 1,
    Hyperbolic = 2,
    Parabolic = 3
};

// enum for the two possibe types of boundary conditions
enum class BoundaryTypes {
    Dirichlet = 1,
    Neumann = 2
};


// Declare the Domain class since both classes depend on each other (same is done in the Domain class header file)
class Domain;


class Initiation {
    protected:
        const CoordinateSystem __coordinate_system;
        const bool __test_case;     //specifies whether the problem should be restricted to a simple test case
        const bool __is_homogeneous;  // i.e. Laplace or Poisson's Equation
        const PDEType __PDE_type;   //Default value for this heat solver is Elliptic (implemented for future solver expansion)

    public:
        Initiation(CoordinateSystem coordinate_system, bool test_case, bool is_homogeneous, PDEType PDE_type);
        
        // Setters
        virtual void setInHomogeneous(const std::unique_ptr<Domain>& domain) = 0;
        virtual void setBCs() = 0;
        
        // Getters
        CoordinateSystem getCoordinateSystem() const;
        bool isTestCase() const;
        bool isHomogeneous() const;
            // pure virtual functions
        virtual vector<BoundaryTypes> getBcTypes() const = 0;
        virtual vector<double> getBcs() const = 0;
        virtual double getHeatSs() const = 0;
        virtual double getK() const = 0;
        virtual array<double,2> getHeatSsLoc() const = 0;


        virtual ~Initiation() = default;
};

#endif