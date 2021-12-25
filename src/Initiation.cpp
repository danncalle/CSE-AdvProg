#include <iostream>

#include "headers/Initiation.h"

Initiation::Initiation(CoordinateSystem coordinate_system, bool test_case, bool is_homogeneous, PDEType PDE_type) 
: __coordinate_system(coordinate_system), __test_case(test_case), __is_homogeneous(is_homogeneous), __PDE_type(PDEType::Elliptic)
{}

CoordinateSystem Initiation::getCoordinateSystem () const {
    return __coordinate_system;
}

bool Initiation::isTestCase() const { return __test_case; }
bool Initiation::isHomogeneous() const { return __is_homogeneous; }