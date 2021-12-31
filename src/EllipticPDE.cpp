#include "./headers/Initiation.h"
#include "./headers/EllipticPDE.h"
#include "./headers/Utilities.h"
#include "./headers/Domain.h"

EllipticPDE::EllipticPDE(CoordinateSystem coordinate_system, bool test_case, bool is_homogeneous)
: Initiation(coordinate_system,test_case,is_homogeneous,PDEType::Elliptic)
{}

void EllipticPDE::setInHomogeneous(const std::unique_ptr<Domain>& domain) {
    if(__is_homogeneous || __test_case)
    return;
    
    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();

    vector<string> message = {"* Enter value of Heat source/sink: ",
                        "* Enter first coordinate for the location of the source/sink (must be inside the specified domain): ",
                        "* Enter second coordinate for the location of the source/sink (must be inside the specified domain): ",
                        "* Enter the thermal conductivity cofficient for this material: "};

    _heat_ss = utils->requestInput('d', static_cast<double>(-INFINITY), static_cast<double>(INFINITY), message[0]);
    if(__coordinate_system == CoordinateSystem::Cartesian) {
        _heat_ss_location[0] = utils->requestInput('d', 0.0, domain->getDimensions()[0], message[1]);
        _heat_ss_location[1] = utils->requestInput('d', 0.0, domain->getDimensions()[1], message[2]);
    
    }
    else if (__coordinate_system == CoordinateSystem::Polar) {
        //TODO: consider the case for circle and oval separaately
        _heat_ss_location[0] = utils->requestInput('d', -domain->getDimensions()[0], domain->getDimensions()[0], message[1]);
        _heat_ss_location[0] = utils->requestInput('d', 0.0, 360.0, message[2]);
    }   
    _k = utils->requestInput('d', static_cast<double>(-INFINITY), static_cast<double>(INFINITY), message[3]);
     
}

void EllipticPDE::setBCs () {
    std::unique_ptr<Utilities> utilities = std::make_unique<Utilities>();
    vector<string> message = {};

    int max_n_BCs = 0; 
    if(__coordinate_system == CoordinateSystem::Cartesian) {
        if(__test_case) {
            max_n_BCs = 2;
            _boundary_types.push_back(BoundaryTypes::Dirichlet);
            _boundary_types[1] = _boundary_types[0];
            _boundary_types[2] = _boundary_types[0];
            _boundary_types[3] = _boundary_types[0];
            message = {"-Bottom, Left, and Right side value of Temp: ",
                "-Top side value of Temp: "};
        }
        else {
            max_n_BCs = 4;
            message = {"-Bottom side BC Type (must be integer, 1 or 2): ",
                "-Left side BC Type (must be integer, 1 or 2): ",
                "-Right side BC Type (must be integer, 1 or 2): ",
                "-Top side BC Type (must be integer, 1 or 2): "};
        }
        
    }
    else if (__coordinate_system == CoordinateSystem::Polar) {
        max_n_BCs = 1;
        message = {"-Circumference BC Type (must be integer, 1 or 2): "};
    }
    
    // REQUESTING TYPES OF BCs
    if(!__test_case) {
        cout << "* Enter the type of BC (Dirichlet = 1, Neumann = 2) for each boundary." << endl;
        for (int i = 0; i < max_n_BCs; i++) {
            _boundary_types.push_back(static_cast<BoundaryTypes>(utilities->requestInput('i', 1, 2, message[i])));
        }
    }

    // REQUESTING VALUES OF BCs
    if(__test_case) cout << "* For each BC enter the value of Temperature." << endl;
    else cout << "* For each BC enter the value of Temperature (Dirichlet) or Heat flux (Neumann)." << endl;
    
    if (__coordinate_system == CoordinateSystem::Polar) {
        message = {"-Circumference BC value: "};
    }
    else {
        message = {"-Bottom side BC value (double): ",
            "-Left side BC value (double): ",
            "-Right side BC value (double): ",
            "-Top side BC value (double): "};
    };
// Fix message so that it indicates which boundary and which type of BC it is
    for (int i = 0; i < max_n_BCs; i++) {
        _boundary_values.push_back(utilities->requestInput('d', static_cast<double>(-INFINITY), static_cast<double>(INFINITY), message[i]));
    }
}

vector<BoundaryTypes> EllipticPDE::getBcTypes(){return _boundary_types;}
vector<double> EllipticPDE::getBcs(){return _boundary_values;}
double EllipticPDE::getHeatSs(){return _heat_ss;}
double EllipticPDE::getK(){return _k;}
array<double,2> EllipticPDE::getHeatSsLoc(){return _heat_ss_location;}