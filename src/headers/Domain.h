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

enum class Shape {
    Square = 1,
    Rectangle = 2,
    Circle = 3,
    Oval = 4
};

class Initiation;
// Both the domain and the ellipticPDE have number of boundaries, i think it is better if 
// only domain has _n_of_boundaries and it is passed as an argument to ellipticPDE where it is set.
// This would recquire to have a get_num_bcs from the domain. Otherwise, just leave num_bcs to the EllipticPDE class.
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