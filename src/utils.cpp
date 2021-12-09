#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <tuple>

#include "./headers/utils.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cin;


template <typename T> void print_vector(const T& v) {
    // Print a vector of any type :)
    for (const auto& elem: v) {
        cout << elem << " ";
    }
    cout << endl;
}

template <typename T> void print_matrix(const T& m) {
    // Print a matrix of any type :)
    for (const auto& vec : m ) {
        print_vector(vec);
    }
}



vector<vector<double>> vec_2_mat(size_t n_x, size_t n_y, const vector<double>& x) {
// This function converts a vector to a 2D matrix. TODO
    
    vector<double> v (n_x,0);
    vector<vector<double>> M (n_y, v);
    return M;

}

