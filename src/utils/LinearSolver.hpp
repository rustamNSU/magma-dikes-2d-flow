#pragma once
#include <Eigen/Sparse>


namespace LinearSolver{
    #if defined(USE_PARDISO)
    void solvePardiso(const Eigen::SparseMatrix<double>& mat, const Eigen::VectorXd& rhs, Eigen::VectorXd& sol);
    #endif
    #if defined(USE_UMFPACK)
    void solveUmfpack(const Eigen::SparseMatrix<double>& mat, const Eigen::VectorXd& rhs, Eigen::VectorXd& sol);
    #endif
}; // LinearSolver