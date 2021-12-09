#ifndef __NUMERICAL_H_
#define __NUMERICAL_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>

using std::vector;

vector<vector <double>> generate_mesh (double x, double y, size_t n_x, size_t n_y, int coordinate_sys);

vector<vector<double>> fdm_mesh(size_t n_x, size_t n_y, int coordinate_sys);

vector<vector<double>> test_LU(const vector<vector<double>>& L, const vector<vector<double>>& U);

vector<vector<double>> test_fdm_error(const vector<double>& x_analytical, const vector<double>& x_fdm);

vector <double> analytical_sol (vector<vector<double>>& mesh, double W, double H, double T1, double T2);

vector<double> fdm_rhs(size_t n_x, size_t n_y, int coordinate_sys,const vector<double>& bc);

vector<double> backward_sustitution(const vector<vector<double>>& M, const vector<double>& b);

vector<double> forward_sustitution(const vector<vector<double>>& M, const vector<double>& b);

vector<double> finite_diff_4_pts(double x, double y, size_t n_x, size_t n_y, int coordinate_sys, const vector<double>& bc);

void LU_factorization(const vector<vector<double>>& M, vector<vector<double>>& L, vector<vector<double>>& U);


#endif