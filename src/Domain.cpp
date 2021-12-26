#include "headers/Domain.h"

Domain::Domain(vector<double> m_dim, int boundaries, Shape shape, const EllipticPDE* pde) 
: _type(pde->getCoordinateSystem())
{   
    // == Cartesian Case ==
    if(_type == CoordinateSystem::Cartesian) {
        _n_of_boundries = 4;
        
        // TODO1: request input for "Shape" 
        // type int
        // min: 1, max 2

        if (_shape == Shape::Square) {
            // TODO2: request ONE input for the side length
            // BUT INPUT AS 2 VAALUES IN ARRAY
        }
        else if (_shape == Shape::Rectangle) {
            for (size_t i = 0; i < 2; i++) {
                // TODO3: request input for each side (only 2)
            }
        }
    }

    // == Polar Case ==
    else if (_type == CoordinateSystem::Polar) {
        _n_of_boundries = 1;

        // TODO 4: request input for "Shape" 
        // type int
        // min: 1 (+2), max 2 (+2)

        if (_shape == Shape::Circle) {
            // TODO 5: request ONE input for the radius length
        }
        else if (_shape == Shape::Oval) {
            for (size_t i = 0; i < 3; i++) {
                // TODO 6: request input for radius, a, b
            }
        }
    }
}

Shape Domain::getShape() const { return _shape; }
vector<double> Domain::getDimensions() const { return _major_dimensions; }