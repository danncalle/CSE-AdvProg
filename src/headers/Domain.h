#ifndef DOMAIN.H
#define DOMAIN.H

#include <vector>
#include <memory>
#include <array>

#include "Initiation.h"
#include "EllipticPDE.h"
#include "Utilities.h"

using std::vector;
using std::array;

enum class Shape {
    Square = 1,
    Rectangle = 2,
    Circle = 3,
    Oval = 4
};

class Initiation;

class Domain {
    private:
        int _n_of_boundries;
        array<double, 2> _major_dimensions = {0,0};

        Shape _shape;
        const CoordinateSystem _type;

    public:
        Domain(const std::unique_ptr<Initiation>& pde);
        Shape getShape() const;
        array<double,2> getDimensions() const;
};      

#endif