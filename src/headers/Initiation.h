#ifndef INITIATION.H
#define INITIATION.H

#include <iostream>
#include <memory>

enum class CoordinateSystem {
    Cartesian = 1,
    Polar = 2
};

enum class PDEType {
    Elliptic = 1,
    Hyperbolic = 2,
    Parabolic = 3
};

class Domain;

class Initiation {
    protected:
        const CoordinateSystem __coordinate_system;
        const bool __test_case;
        const bool __is_homogeneous;
        const PDEType __PDE_type;   //Default value for this heat solver (implemented for future solver expansion)

    public:
        Initiation(CoordinateSystem coordinate_system, bool test_case, bool is_homogeneous, PDEType PDE_type);
        CoordinateSystem getCoordinateSystem() const;
        bool isTestCase() const;
        bool isHomogeneous() const;

        virtual void setInHomogeneous(const std::unique_ptr<Domain>& domain) = 0;
        virtual void setBCs() = 0;

        virtual ~Initiation() = default;
};

#endif