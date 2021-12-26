#include "headers/EllipticPDE.h"
#include "headers/Utilities.h"

EllipticPDE::EllipticPDE(CoordinateSystem coordinate_system, bool test_case, bool is_homogeneous, PDEType PDE_type, double heat_ss, array <double, 2> location, double k)
: Initiation(coordinate_system,test_case,is_homogeneous,PDE_type), _heat_ss(heat_ss), _heat_ss_location(location), _k(k)
{}

void EllipticPDE::setBCs () {
    std::unique_ptr<Utilities> utilities = std::make_unique<Utilities>();
    vector<string> m = {};

    int max_n_BCs = 0; 
    if(__coordinate_system == CoordinateSystem::Cartesian) {
        max_n_BCs = 4;
        m = {"-Lower side BC Type (must be integer, 1 or 2): ",
            "-Left side BC Type (must be integer, 1 or 2): ",
            "-Right side BC Type (must be integer, 1 or 2): ",
            "-Top side BC Type (must be integer, 1 or 2): "};
    }
    else if (__coordinate_system == CoordinateSystem::Polar) {
        max_n_BCs = 1;
        m = {"-Circumference BC Type (must be integer, 1 or 2): "};
    }
    
    // REQUESTING TYPES OF BCs
    cout << "* Enter the type of BC (Dirichlet = 1, Neumann = 2) for each boundary." << endl;
    for (int i = 0; i < max_n_BCs; i++) {
        _boundary_types.push_back(static_cast<BoundaryTypes>(utilities->requestInput('i', 1, 2, m[i])));
    }
    
    // REQUESTING VALUES OF BCs
    cout << "* For each BC enter the value of Temperature (Dirichlet) or Heat flux (Neumann)." << endl;
    double infinity = std::numeric_limits<double>::infinity();
    for (int i = 0; i < max_n_BCs; i++) {
        _boundary_values.push_back(utilities->requestInput('d', -infinity, infinity, "m[i]"));
    }
}