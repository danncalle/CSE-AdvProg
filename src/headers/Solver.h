#ifndef __SOLVER_H_
#define __SOLVER_H_

#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <memory>

#include "Mesh.h"

using std::vector;
using std::array;

enum class SolverType {
    Numerical = 1,
    Analytical = 2
};

class Solver {
    protected:
        const SolverType __solver_type;
        const bool __test_case;
        vector<vector<double>> __solution {{},{}};
        vector<double> __error {};

    public:
        Solver(SolverType solver_type, bool test_case);
        
        vector<double> getError() const;
        vector<vector<double>> getSolution() const;

        virtual void setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh) = 0;
        virtual void setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh) = 0;
        virtual void solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh) = 0;
        virtual void setError(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain) = 0;

        virtual ~Solver() = default;
};

class LU : public Solver{
    protected:
        vector<double> __b;
        vector<vector<double>> __M;
        vector<vector<double>> __L;
        vector<vector<double>> __U;

        vector<double> backwardSubstitution(const vector<double>& rhs);
        vector<double> forwardSubstitution(const vector<double>& rhs);
        void LUFactorization();
    
    public:
        LU(const bool test_case);

        void setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
        void setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
        void solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
        void setError(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain);
};

// class Seidel : public Solver{
//     protected:
//         vector<double> b;
//         vector<vector<double>> M;
//         vector<vector<double>> L;
//         vector<vector<double>> U;
    
//     public:
//         Seidel(const bool test_case);

//         void setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
//         void setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
//         void solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
//         void setError(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain);
// };

class Analytic : public Solver{
    protected:
        vector<double> __analytic_sol;
    public:
        Analytic();

        void solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain);
        vector<double> getAnalyticSol() const;

        void setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
        void setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
        void setError(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain);
        void solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);

};

#endif