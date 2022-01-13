#ifndef __DOMAIN_H_
#define __DOMAIN_H_

#include <vector>
#include <memory>
#include <array>

#include "Initiation.h"
#include "EllipticPDE.h"
#include "Utilities.h"

using std::vector;
using std::array;

// Enum to define the shape to be selected for the domain (will be restricted based on the chosen coordinate)
enum class Shape {
    Square = 1,
    Rectangle = 2,
    Circle = 3,
    Oval = 4
};

// Declare the Initiation class since both classes depend on each other (same is done in the Initiation class header file)
class Initiation;

class Domain {
    private:
        int _n_of_boundries;  // no. of boundaries of the domain depends on the chosen shape and is defined 
                              // automatically after the domain type is chosen
        array<double, 2> _major_dimensions = {0,0}; // `W` and `H` for cartesian. only `r` for circle. `a` and `b` for oval.

        Shape _shape;
        const CoordinateSystem _type;

    public:
        Domain(const std::unique_ptr<Initiation>& pde);    // Initialize the domain based on a certain PDE type
        
        // Getters
        Shape getShape() const;
        array<double,2> getDimensions() const;
        double getRatTheta(double theta) const;
};      

#endif