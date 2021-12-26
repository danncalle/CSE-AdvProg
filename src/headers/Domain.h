#ifndef DOMAIN.H
#define DOMAIN.H

#include <vector>

#include "EllipticPDE.h"

using std::vector;

enum class Shape {
    Square = 1,
    Rectangle = 2,
    Circle = 3,
    Oval = 4
};

class Domain {
    private:
        int _n_of_boundries;
        vector<double> _major_dimensions;

        Shape _shape;
        const CoordinateSystem _type;

    public:
        Domain(vector<double> m_dim, int boundaries, Shape shape, const EllipticPDE* pde);
        Shape getShape() const;
        vector<double> getDimensions() const;
};      

#endif