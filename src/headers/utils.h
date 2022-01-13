#ifndef __UTILS_H_
#define __UTILS_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <algorithm>

using std::vector;

template <typename T> void print_vector(const T& v);

template <typename T> void print_matrix(const T& m);

vector<vector<double>> vec_2_mat(size_t n_x, size_t n_y, const vector<double>& x);

// std::ostream &operator << (std::ostream &out, const vector<double> &v);

// std::ostream &operator << (std::ostream &out, const vector<vector<double>> &M);

#endif