#ifndef UTILITIES.H
#define UTILITIES.H

#include <iostream>
#include <vector>

class Utilities {
    public:
        template <typename T> void print_vector(const T& v);
        template <typename T> void print_matrix(const T& m);
};

#endif