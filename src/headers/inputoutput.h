#ifndef __INPUTOUTPUT_H_
#define __INPUTOUTPUT_H_

#include <vector>

using std::vector;


void initiateProgram ();

std::tuple<vector<double>, vector<int>, bool, vector<double>> requestInputs ();

void outputResults (const vector<vector<double>>& mesh, const vector<double>& t_sim, const vector<double>& t_analytic = {}, const vector<double>& err ={});

void plotResults(const vector<double>& dimensions, const vector<int>& nodes, const vector<double>& t_sim);

#endif