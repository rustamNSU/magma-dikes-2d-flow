#pragma once
#include <Eigen/Sparse>


namespace LinearSolver{
    struct solverInfo{
        bool is_converge = false;
    };
    solverInfo solveUmfpack(const Eigen::SparseMatrix<double>& mat, const Eigen::VectorXd& rhs, Eigen::VectorXd& sol);
}; // LinearSolver