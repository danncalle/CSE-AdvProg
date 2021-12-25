#include "headers/Domain.h"

Domain::Domain(vector<double> m_dim, int boundaries, Shape shape, Type type) 
: _type(type)
{   
    // == Cartesian Case ==
    if(_type == Type::Cartesian) {
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
    else if (_type == Type::Polar) {
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

Type Domain::getType() const { return _type; }
Shape Domain::getShape() const { return _shape; }
vector<double> Domain::getDimensions() const { return _major_dimensions; }