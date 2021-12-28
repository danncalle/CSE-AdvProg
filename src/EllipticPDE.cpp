#include "headers/Initiation.h"
#include "headers/EllipticPDE.h"
#include "headers/Utilities.h"
#include "headers/Domain.h"

EllipticPDE::EllipticPDE(CoordinateSystem coordinate_system, bool test_case, bool is_homogeneous)
: Initiation(coordinate_system,test_case,is_homogeneous,PDEType::Elliptic)
{}

void EllipticPDE::setInHomogeneous(const std::unique_ptr<Domain>& domain) {
    if(__is_homogeneous || __test_case)
    return;
    
    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();

    vector<string> m = {"* Enter value of Heat source/sink: ",
                        "* Enter first coordinate for the location of the source/sink (must be inside the specified domain): ",
                        "* Enter second coordinate for the location of the source/sink (must be inside the specified domain): ",
                        "* Enter the thermal conductivity cofficient for this material: "};

    _heat_ss = utils->requestInput('d', static_cast<double>(-INFINITY), static_cast<double>(INFINITY), m[0]);
    if(__coordinate_system == CoordinateSystem::Cartesian) {
        _heat_ss_location[0] = utils->requestInput('d', 0.0, domain->getDimensions()[0], m[1]);
        _heat_ss_location[1] = utils->requestInput('d', 0.0, domain->getDimensions()[1], m[2]);
    
    }
    else if (__coordinate_system == CoordinateSystem::Polar) {
        //TODO: consider the case for circle and oval separaately
        _heat_ss_location[0] = utils->requestInput('d', -domain->getDimensions()[0], domain->getDimensions()[0], m[1]);
        _heat_ss_location[0] = utils->requestInput('d', 0.0, 360.0, m[2]);
    }   
    _k = utils->requestInput('d', static_cast<double>(-INFINITY), static_cast<double>(INFINITY), m[3]);
     
}

void EllipticPDE::setBCs () {
    std::unique_ptr<Utilities> utilities = std::make_unique<Utilities>();
    vector<string> m = {};

    int max_n_BCs = 0; 
    if(__coordinate_system == CoordinateSystem::Cartesian) {
        if(__test_case) {
            max_n_BCs = 2;
            _boundary_types.push_back(BoundaryTypes::Dirichlet);
            _boundary_types[1] = _boundary_types[0];
            _boundary_types[2] = _boundary_types[0];
            _boundary_types[3] = _boundary_types[0];
            m = {"-Bottom, Left, and Right side value of Temp: ",
                "-Top side value of Temp: "};
        }
        else {
            max_n_BCs = 4;
            m = {"-Bottom side BC Type (must be integer, 1 or 2): ",
                "-Left side BC Type (must be integer, 1 or 2): ",
                "-Right side BC Type (must be integer, 1 or 2): ",
                "-Top side BC Type (must be integer, 1 or 2): "};
        }
        
    }
    else if (__coordinate_system == CoordinateSystem::Polar) {
        max_n_BCs = 1;
        m = {"-Circumference BC Type (must be integer, 1 or 2): "};
    }
    
    // REQUESTING TYPES OF BCs
    if(!__test_case) {
        cout << "* Enter the type of BC (Dirichlet = 1, Neumann = 2) for each boundary." << endl;
        for (int i = 0; i < max_n_BCs; i++) {
            _boundary_types.push_back(static_cast<BoundaryTypes>(utilities->requestInput('i', 1, 2, m[i])));
        }
    }

    // REQUESTING VALUES OF BCs
    if(__test_case) cout << "* For each BC enter the value of Temperature." << endl;
    else cout << "* For each BC enter the value of Temperature (Dirichlet) or Heat flux (Neumann)." << endl;
    
    if (__coordinate_system == CoordinateSystem::Polar)  m = {"-Circumference BC value: "};

    for (int i = 0; i < max_n_BCs; i++) {
        _boundary_values.push_back(utilities->requestInput('d', static_cast<double>(-INFINITY), static_cast<double>(INFINITY), m[i]));
    }
}