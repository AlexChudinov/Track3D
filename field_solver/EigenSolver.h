#pragma once
#ifndef _EIGEN_SOLVER_H_
#define _EIGEN_SOLVER_H_

#include <Eigen/IterativeLinearSolvers>

//Solver uses Eigenlib for solving laplacian problems
class EigenSolver : public COperator
{
public:
	EigenSolver(const CFieldOperator::Matrix& matrix, const CMeshAdapter& mesh);
	~EigenSolver();

	virtual void applyToField(Field& f, double * tol = nullptr) const;

	virtual Hist applyToFieldNTimes(Field& f0, size_t N, ThreadPool::Progress * p) const;

private:

	Eigen::VectorXd mB;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> mSolver;
};

#endif