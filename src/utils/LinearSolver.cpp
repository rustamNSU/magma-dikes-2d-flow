#include "LinearSolver.hpp"
#include <Eigen/UmfPackSupport>


LinearSolver::solverInfo LinearSolver::solveUmfpack(const Eigen::SparseMatrix<double>& mat, const Eigen::VectorXd& rhs, Eigen::VectorXd& sol){
    Eigen::UmfPackLU<Eigen::SparseMatrix<double>> eigen_solver;
    eigen_solver.analyzePattern(mat);
	eigen_solver.factorize(mat);
	sol = eigen_solver.solve(rhs);
	LinearSolver::solverInfo info;
    if (eigen_solver.info() != Eigen::Success) info.is_converge = false;
	else info.is_converge = true;
    return info;
}