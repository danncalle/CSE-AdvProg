#include "headers/Utilities.h"

using std::cout;
using std::endl;

template <typename T> void Utilities::print_vector(const T& v) {
    // Print a vector of any type :)
    for (const auto& elem: v) {
        cout << elem << " ";
    }
    cout << endl;
}


 template <typename T> void Utilities::print_matrix(const T& m) {
    // Print a matrix of any type :)
    for (const auto& vec : m ) {
        print_vector(vec);
    }
}