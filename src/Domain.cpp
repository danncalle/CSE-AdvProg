#include "headers/Domain.h"

Domain::Domain(const std::unique_ptr<Initiation>& pde) 
: _type(pde->getCoordinateSystem())
{       
    std::unique_ptr<Utilities> utils = std::make_unique<Utilities>();
    // == Cartesian Case ==
    if(_type == CoordinateSystem::Cartesian) {
        _n_of_boundries = 4;
        
        _shape = static_cast<Shape>(utils->requestInput('i',1, 2, "* Enter shape type (Square = 1, Rectangle = 2): "));

        if (_shape == Shape::Square) {
            _major_dimensions[0] = utils->requestInput('d', 0.0, static_cast<double>(INFINITY), "* Enter side length: ");
            _major_dimensions[1] = _major_dimensions[0];
        }
        else if (_shape == Shape::Rectangle) {
            vector<string> m = {"* Enter width (>0): ", "* Enter height (>0): "};
            for (size_t i = 0; i < 2; i++) {
                _major_dimensions[i] = utils->requestInput('d', 0.0, static_cast<double>(INFINITY), m[i]);        
            }
        }
    }

    // == Polar Case ==
    else if (_type == CoordinateSystem::Polar) {
        _n_of_boundries = 1;

        _shape = static_cast<Shape>(utils->requestInput('i',1, 2, "* Enter shape type (Circle = 1, Oval = 2): ")+2);

        if (_shape == Shape::Circle) {
            _major_dimensions[0] = utils->requestInput('d', 0.0, static_cast<double>(INFINITY), "* Enter radius length (>0): ");
        }
        else if (_shape == Shape::Oval) {
            for (size_t i = 0; i < 2; i++) {
                vector<string> m = {"* Enter radius in horizontal direction (>0): ", "* Enter radius in vertical direction (>0): "};
                _major_dimensions[i] = utils->requestInput('d', 0.0, static_cast<double>(INFINITY), m[i]);
            }
        }
    }
}

Shape Domain::getShape() const { return _shape; }
array<double,2> Domain::getDimensions() const { return _major_dimensions; }