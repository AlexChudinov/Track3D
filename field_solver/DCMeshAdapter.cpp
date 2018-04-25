#include "stdafx.h"
#include "AnsysMesh.h"
#include "DCMeshAdapter.h"

DCMeshAdapter::DCMeshAdapter(const Elements & es, const Nodes & ns, const DirTess & tess)
	:
	CMeshAdapter(es, ns),
	m_tess(tess)
{
}

CMeshAdapter::PScalFieldOp DCMeshAdapter::createOperator(ScalarOperatorType type, const BaseOperatorParams*)
{
#define CREATE_GRAD_OP(X) PScalFieldOp(new ScalarFieldOperator(grad##X()))
	switch (type)
	{
	case LaplacianSolver: return PScalFieldOp(new ScalarFieldOperator(laplacian()));
	case LaplacianSolver1: return PScalFieldOp(new ScalarFieldOperator(laplacian1()));
	case GradX: return CREATE_GRAD_OP(X);
	case GradY: return CREATE_GRAD_OP(Y);
	case GradZ: return CREATE_GRAD_OP(Z);
	}
	return PScalFieldOp();
}


#define SET_NODE_GRAD_FUN(X, x)\
CMeshAdapter::InterpCoefs DCMeshAdapter::grad##X(Label idx) const\
{\
	InterpCoefs res;\
	const Node* pNode = m_nodes[idx];\
	if (pNode->vNbrNodes.size() == m_tess.get_cell(idx)->nFaceCount)\
	{\
		double fNorm = 0.0;\
\
		for (uint32_t nFaceIdx = 0; nFaceIdx < pNode->vNbrNodes.size(); ++nFaceIdx)\
		{\
			size_t l = pNode->vNbrNodes[nFaceIdx];\
			Vector3D n = (m_nodes[l]->pos - m_nodes[idx]->pos); n.normalize();\
			double fWeight = n.x * m_tess.get_cell(idx)->pFaceSquare[nFaceIdx];\
			fNorm += fWeight;\
			res[l] = fWeight;\
		}\
		res[idx] = fNorm;\
		return mul(.5/m_tess.get_cell(idx)->fVolume, res);\
	}\
	else return CMeshAdapter::grad##X(idx);\
}

SET_NODE_GRAD_FUN(X,x)
SET_NODE_GRAD_FUN(Y,y)
SET_NODE_GRAD_FUN(Z,z)

#define SET_GRAD_FUN(X) \
CMeshAdapter::ScalarFieldOperator DCMeshAdapter::grad##X() const\
{\
	progressBar()->set_job_name("Creating gradient " #X "-component operator...");\
	progressBar()->set_progress(0);\
\
	ScalarFieldOperator result;\
	result.m_matrix.resize(m_nodes.size());\
\
	ThreadPool::splitInPar(m_nodes.size(),\
		[&](size_t nNodeIdx)->void { result.m_matrix[nNodeIdx] = grad##X(nNodeIdx); },\
		progressBar());\
\
	return result;\
}

SET_GRAD_FUN(X)
SET_GRAD_FUN(Y)
SET_GRAD_FUN(Z)


CMeshAdapter::ScalarFieldOperator DCMeshAdapter::laplacian() const
{
	progressBar()->set_job_name("Creating laplacian #1 operator...");
	progressBar()->set_progress(0);

	ScalarFieldOperator result;
	result.m_matrix.resize(m_nodes.size());
	std::vector<NodeType> vNodeTypes(nodeTypes());

	ThreadPool::splitInPar(m_nodes.size(),
		[&](size_t nNodeIdx)->void 
	{
		switch (vNodeTypes[nNodeIdx])
		{
		case FirstTypeBoundaryNode: //Fixed value BC
			result.m_matrix[nNodeIdx][nNodeIdx] = 1.0;
			break;
		case SecondTypeBoundaryNode: //Zero gradient
		{
			const Vector3D& n = boundaryMesh()->normal(nNodeIdx);
			InterpCoefs coefs = std::move(add(
				mul(n.x, gradX(nNodeIdx)),
				add(mul(n.y, gradY(nNodeIdx)),
					mul(n.z, gradZ(nNodeIdx)))));
			double fSum = coefs[nNodeIdx];
			coefs.erase(nNodeIdx);
			mul(-1. / fSum, coefs);
			result.m_matrix[nNodeIdx] = std::move(coefs);
			break;
		}
		default: //It is inner point
		{
			InterpCoefs coefs;
			double fNorm = 0.0;

			for (uint32_t nFaceIdx = 0; nFaceIdx < m_nodes[nNodeIdx]->vNbrNodes.size(); ++nFaceIdx)
			{
				size_t l = m_nodes[nNodeIdx]->vNbrNodes[nFaceIdx];
				double fWeight = m_tess.get_cell(nNodeIdx)->pFaceSquare[nFaceIdx]
					/ (m_nodes[l]->pos - m_nodes[nNodeIdx]->pos).length();
				coefs[l] = fWeight;
				fNorm += fWeight;
			}
			mul(1. / fNorm, coefs);
			result.m_matrix[nNodeIdx] = std::move(coefs);
			break;
		}
		}
	},
		progressBar());

	return result;
}

CMeshAdapter::ScalarFieldOperator DCMeshAdapter::laplacian1() const
{
	ScalarFieldOperator
		opGradX = gradX(), opGradY = gradY(), opGradZ = gradZ(), opLaplace;
	progressBar()->set_job_name("Create laplacian #2 operator");
	opLaplace.m_matrix.resize(m_nodes.size());

	std::vector<NodeType> vNodeTypes(nodeTypes());

	ThreadPool::splitInPar(m_nodes.size(),
		[&](size_t nCurNodeIdx)
	{
		switch (vNodeTypes[nCurNodeIdx])
		{
		case FirstTypeBoundaryNode: //Fixed value BC
			opLaplace.m_matrix[nCurNodeIdx][nCurNodeIdx] = 1.0;
			break;
		case SecondTypeBoundaryNode: //Zero gradient
		{
			const Vector3D& n = boundaryMesh()->normal(nCurNodeIdx);
			InterpCoefs coefs = std::move(add(
				mul(n.x, opGradX.m_matrix[nCurNodeIdx]),
				add(mul(n.y, opGradY.m_matrix[nCurNodeIdx]),
					mul(n.z, opGradZ.m_matrix[nCurNodeIdx]))));
			double fSum = coefs[nCurNodeIdx];
			coefs.erase(nCurNodeIdx);
			mul(-1. / fSum, coefs);
			opLaplace.m_matrix[nCurNodeIdx] = std::move(coefs);
			break;
		}
		default:
		{
			InterpCoefs coefsX, coefsY, coefsZ;
			for (const auto& c : opGradX.m_matrix[nCurNodeIdx])
			{
				add(coefsX, mul(c.second, opGradX.m_matrix[c.first]));
			}
			for (const auto& c : opGradY.m_matrix[nCurNodeIdx])
			{
				add(coefsY, mul(c.second, opGradY.m_matrix[c.first]));
			}
			for (const auto& c : opGradZ.m_matrix[nCurNodeIdx])
			{
				add(coefsZ, mul(c.second, opGradZ.m_matrix[c.first]));
			}
			opLaplace.m_matrix[nCurNodeIdx] = add(coefsX, add(coefsY, coefsZ));
			double fSum = opLaplace.m_matrix[nCurNodeIdx][nCurNodeIdx];
			opLaplace.m_matrix[nCurNodeIdx].erase(nCurNodeIdx);
			mul(-1. / fSum, opLaplace.m_matrix[nCurNodeIdx]);
		}
		}

	},
		progressBar());


	return opLaplace;
}
