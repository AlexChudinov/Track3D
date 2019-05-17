#include "stdafx.h"
#include "MeshData.h"
#include "EigenSolver.h"

EigenSolver::EigenSolver(const CFieldOperator::Matrix & matrix, const CMeshAdapter& mesh)
	:
	mB(matrix.size())
{
	const size_t n = matrix.size();
	Eigen::SparseMatrix<double> local(n, n);
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(n);
	for (size_t i = 0; i < n; ++i)
	{
		const Eigen::Index nRow = static_cast<Eigen::Index>(i);
		for (const CFieldOperator::MatrixCoef& c : matrix[i])
		{
			const Eigen::Index nCol = static_cast<Eigen::Index>(c.first);
			triplets.push_back(Eigen::Triplet<double>(nRow, nCol, c.second));
		}
	}
	local.setFromTriplets(triplets.begin(), triplets.end());
	Field f(matrix.size());
	mesh.boundaryMesh()->applyBoundaryVals(f);
	for (size_t i = 0; i < f.size(); ++i)
	{
		mB(static_cast<Eigen::Index>(i)) = f[i];
	}
	mSolver.compute(local);
}

EigenSolver::~EigenSolver()
{
}

void EigenSolver::applyToField(Field & f, double * tol) const
{
	Eigen::VectorXd x0(f.size());
	for (size_t i = 0; i < f.size(); ++i)
	{
		x0(static_cast<Eigen::Index>(i)) = f[i];
	}
	x0 = mSolver.solveWithGuess(mB, x0);
	double xMax = DBL_MIN;
	for (size_t i = 0; i < f.size(); ++i)
	{
		f[i] = x0(static_cast<Eigen::Index>(i));
		xMax = max(f[i], xMax);
	}
	if (tol)
	{
		*tol = mSolver.error() / xMax;
	}

}

EigenSolver::Hist EigenSolver::applyToFieldNTimes(Field & f0, size_t N, ThreadPool::Progress * p) const
{
	return Hist();
}
