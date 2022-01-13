#ifndef __SOLVER_H_
#define __SOLVER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <memory>
#include <eigen3/Eigen/Sparse>

#include "Utilities.h"
#include "Mesh.h"

using std::vector;
using std::array;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

enum class SolverType {
    LU_dir = 1,
    LU_sp = 2,
    Gauss_Seidel = 3
};
// Only have the setupRHS and setupMatrix in the LU solver class
// Give the user to load the LU matrix for a problem and we change the RHS and solve it.
class Solver {
    protected:
        const SolverType __solver_type;
        const bool __test_case;
        vector<vector<double>> __solution {{},{}};
        vector<double> __error {};

        void solve_analytical(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain);
    public:
        Solver(SolverType solver_type, const bool test_case);
        
        vector<double> getError() const;
        vector<vector<double>> getSolution() const;

        void setError(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh, const std::unique_ptr<Domain>& domain);
        virtual void solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh) = 0;
        virtual ~Solver() = default;
};

class LU_direct : public Solver{
    protected:
        vector<double> __b;
        vector<vector<double>> __M;
        vector<vector<double>> __L;
        vector<vector<double>> __U;

        vector<double> backwardSubstitution(const vector<double>& rhs);
        vector<double> forwardSubstitution(const vector<double>& rhs);
        void LUFactorization();
        void setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
        void setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);

    public:
        LU_direct(const bool test_case);

        void writeLU();
        void loadLU();
        void solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
};

class LU_sparse : public Solver{
    protected:
        VectorXd __b;
        SparseMatrix<double> __M;
        // SparseMatrix<double> __L;
        // SparseMatrix<double> __U;
        // SparseMatrix<double> __D;

        void setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
        void setupMatrix(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);

    public:
        LU_sparse(const bool test_case);

        // void writeLU();
        // void loadLU();

        void solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
};

class Seidel : public Solver{
    protected:
        vector<double> __b;
        size_t __max_iter;
        double __tol;
        double __omega;

        void setupRhs(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
    
    public:
        Seidel(const bool test_case);

        void setupSolver();
        void solve(const std::unique_ptr<Initiation>& pde, const std::unique_ptr<Mesh>& mesh);
};

#endif