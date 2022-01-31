#include "./headers/Initiation.h"
#include "./headers/EllipticPDE.h"
#include "./headers/Utilities.h"
#include "./headers/Domain.h"

EllipticPDE::EllipticPDE(CoordinateSystem coordinate_system, bool test_case, bool is_homogeneous)
: Initiation(coordinate_system,test_case,is_homogeneous,PDEType::Elliptic)
{}

void EllipticPDE::setInHomogeneous(const std::unique_ptr<Domain>& domain) {
    // @ set the location and the value of the heat source/sink for an inhomoegeneous problem

    /* if the homogeous or a simple test case is selected then directly exit this function */
    if(__is_homogeneous || __test_case)
    return;
    
    // use some functionalities from the Utilities class
    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();

    // storing the messages to be used for the next few input requests
    vector<string> message = {"* Enter value of Heat source/sink: ",
                        "* Enter first coordinate for the location of the source/sink (must be inside the specified domain): ",
                        "* Enter second coordinate for the location of the source/sink (must be inside the specified domain): ",
                        "* Enter the thermal conductivity cofficient for this material: "};

    // request value for heat source/sink
    _heat_ss = utils->requestInput('d', static_cast<double>(-INFINITY), static_cast<double>(INFINITY), message[0]);
    
    // request the location for the heat source/sink. this provided location must be inside the defined domain.
    // this applies separately for the Cartesian or the Polar system
    if(__coordinate_system == CoordinateSystem::Cartesian) {
        // Note: The bottom left corner of the rectangle/square is fixed at the origin.
        _heat_ss_location[0] = utils->requestInput('d', 0.0, domain->getDimensions()[0], message[1]);
        _heat_ss_location[1] = utils->requestInput('d', 0.0, domain->getDimensions()[1], message[2]);
    
    }
    else if (__coordinate_system == CoordinateSystem::Polar) {
        
        // request the angle first!
        _heat_ss_location[1] = utils->requestInput('d', 0.0, 360.0, message[2]);
        if(domain->getShape() == Shape::Circle) {
            // should be between 0 and r
            _heat_ss_location[0] = utils->requestInput('d', 0.0, domain->getDimensions()[0], message[1]);
        }
        else if (domain->getShape() == Shape::Oval) {
            // should be between 0 and r (Note: r is changing for each theta!)
            _heat_ss_location[0] = utils->requestInput('d', 0.0, domain->getRatTheta(_heat_ss_location[1]), message[1]);
        }
    }

    // request Coefficient of thermal conductivity   
    _k = utils->requestInput('d', static_cast<double>(-INFINITY), static_cast<double>(INFINITY), message[3]);
     
}

void EllipticPDE::setBCs () {
    // @ Populate the _boundary_types and _boundary_values based on user input

    // use some functionalities from the Utilities class
    std::unique_ptr<Utilities> utilities = std::make_unique<Utilities>();

    // storing messages that will be used in the input requets below
    vector<string> message = {};
    vector<string> message2 = {};

    // no. of BCs based on the shape
    int max_n_BCs = 0;

    /* == Cartesian case == */ 
    if(__coordinate_system == CoordinateSystem::Cartesian) {
        /*-- simple test case -- */
        if(__test_case) {
            max_n_BCs = 2;
            for (int i=0;i<4;i++) _boundary_types.push_back(BoundaryTypes::Dirichlet);
            message2 = {"-Bottom, Left, and Right side value of Temp: ",
                "-Top side value of Temp: "};
        }
        /*-- square/rectangle -- */
        else {
            max_n_BCs = 4;
            message = {"-Bottom side BC Type (must be integer, 1 or 2): ",
                "-Left side BC Type (must be integer, 1 or 2): ",
                "-Right side BC Type (must be integer, 1 or 2): ",
                "-Top side BC Type (must be integer, 1 or 2): "};

            message2 = {"-Bottom side BC value (double): ",
                "-Left side BC value (double): ",
                "-Right side BC value (double): ",
                "-Top side BC value (double): "};    
        }
    }
    /* == Polar case == */ 
    else if (__coordinate_system == CoordinateSystem::Polar) {
        max_n_BCs = 2;
        for (int i=0;i<max_n_BCs;i++) _boundary_types.push_back(BoundaryTypes::Dirichlet);
        // message = {"- Inner-most circle BC Type (must be integer, 1 or 2): ",
                    // "- Circumference BC Type (must be integer, 1 or 2):"};
    }
    
    // REQUESTING TYPES OF BCs
    if(!(__test_case || __coordinate_system == CoordinateSystem::Polar)) {
        cout << "* Enter the type of BC (Dirichlet = 1, Neumann = 2) for each boundary." << endl;
        for (int i = 0; i < max_n_BCs; i++) {
            _boundary_types.push_back(static_cast<BoundaryTypes>(utilities->requestInput('i', 1, 2, message[i])));
        }
    }

    // REQUESTING VALUES OF BCs
    if(__test_case || __coordinate_system == CoordinateSystem::Polar) cout << "* For each BC enter the value of Temperature." << endl;
    else cout << "* For each BC enter the value of Temperature (Dirichlet) or Heat flux (Neumann)." << endl;
    
    if (__coordinate_system == CoordinateSystem::Polar) {
        message2 = {"-Inner-most circle BC value: ","-Circumference BC value: "};
    }
    for (int i = 0; i < max_n_BCs; i++) {
        _boundary_values.push_back(utilities->requestInput('d', static_cast<double>(-INFINITY), static_cast<double>(INFINITY), message2[i]));
    }
}


// Getters
vector<BoundaryTypes> EllipticPDE::getBcTypes() const {return _boundary_types;}
vector<double> EllipticPDE::getBcs() const {return _boundary_values;}
double EllipticPDE::getHeatSs() const {return _heat_ss;}
double EllipticPDE::getK() const {return _k;}
array<double,2> EllipticPDE::getHeatSsLoc() const {return _heat_ss_location;}