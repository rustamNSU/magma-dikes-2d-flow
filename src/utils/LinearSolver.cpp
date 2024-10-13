#include "LinearSolver.hpp"

#if defined(USE_PARDISO)
#define EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif // defined(USE_PARDISO)

#if defined(USE_UMFPACK)
#include <Eigen/UmfPackSupport>
#endif


#if defined(USE_PARDISO)
void LinearSolver::solvePardiso(const Eigen::SparseMatrix<double>& mat, const Eigen::VectorXd& rhs, Eigen::VectorXd& sol){
    Eigen::PardisoLU<Eigen::SparseMatrix<double>> eigen_solver;
    auto pardiso_parameters = eigen_solver.pardisoParameterArray();
	pardiso_parameters[0] = 1; // Use default values if 0
	pardiso_parameters[1] = 2; // Fill-in reducing ordering for the input matrix.
	pardiso_parameters[2] = 0; // Reserved. Set to zero (num_procs)
	pardiso_parameters[3] = 0;
	pardiso_parameters[7] = 2; // Max numbers of iterative refinement steps
	pardiso_parameters[10] = 1; // Scaling vectors
	pardiso_parameters[12] = 1; // mproved accuracy using (non-) symmetric weighted matching
	pardiso_parameters[17] = -1; // Report the number of non-zero elements in the factors.
	pardiso_parameters[18] = -1; // Report number of floating point operations (in 106 floating point operations) that are necessary to factor the matrix A. Enable report if iparm[18] < 0 on entry. This increases the reordering time.
	pardiso_parameters[23] = 10; // Parallel factorization control.
	pardiso_parameters[24] = 2; // Parallel Forward/Backward Solve
	pardiso_parameters[26] = 1; // Matrix checker,  PARDISO checks whether column indices are sorted in increasing order within each row
	pardiso_parameters[27] = 0; // Input arrays (a, x and b) and all internal arrays must be presented in double precision
	pardiso_parameters[33] = 0; // Optimal number of OpenMP threads for conditional numerical reproducibility (CNR) mode
	pardiso_parameters[59] = 0;
    eigen_solver.analyzePattern(mat);
	eigen_solver.factorize(mat);
	sol = eigen_solver.solve(rhs);
    if (eigen_solver.info() != Eigen::Success){
		throw std::invalid_argument("pardiso linear solver doesn't converge\n");
	}
    return;
}
#endif

#if defined(USE_UMFPACK)
void LinearSolver::solveUmfpack(const Eigen::SparseMatrix<double>& mat, const Eigen::VectorXd& rhs, Eigen::VectorXd& sol){
    Eigen::UmfPackLU<Eigen::SparseMatrix<double>> eigen_solver;
    eigen_solver.analyzePattern(mat);
	eigen_solver.factorize(mat);
	sol = eigen_solver.solve(rhs);
    if (eigen_solver.info() != Eigen::Success){
		throw std::invalid_argument("UMFPACK linear solver doesn't converge\n");
	}
    return;
}
#endif