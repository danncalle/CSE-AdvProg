#ifndef DOMAIN.H
#define DOMAIN.H

#include <vector>

using std::vector;

enum class Type {
    Cartesian = 1,
    Polar = 2
};

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
        const Type _type;

    public:
        Domain(vector<double> m_dim, int boundaries, Shape shape, Type type);
        Type getType() const;
        Shape getShape() const;
        vector<double> getDimensions() const;
};      

#endif