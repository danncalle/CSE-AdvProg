#include "headers/EllipticPDE.h"

EllipticPDE::EllipticPDE(CoordinateSystem coordinate_system, bool test_case, bool is_homogeneous, PDEType PDE_type, double heat_ss, array <double, 2> location, double k)
: Initiation(coordinate_system,test_case,is_homogeneous,PDE_type), _heat_ss(heat_ss), _heat_ss_location(location), _k(k)
{}

void EllipticPDE::setBCs () {
    // TODO: REQUEST TYPES  (using utilities class)
    _boundary_types.push_back({});
    // TODO: REQUEST VALUES (using utilities class)
    _boundary_values.push_back({});
}