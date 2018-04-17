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
	double fNorm = 0.0;\
\
	for (Label l : m_nodes[idx]->vNbrNodes)\
	{\
		Vector3D n = (m_nodes[l]->pos - m_nodes[idx]->pos); n.normalize();\
		double fWeight = n.x * m_tess.get_cell(idx)->pFaceSquare[l];\
		fNorm += fWeight;\
		res[l] = fWeight;\
	}\
	res[idx] = fNorm;\
	return mul(.5/m_tess.get_cell(idx)->fVolume, res);\
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
	progressBar()->set_job_name("Creating laplacian operator...");
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
			for (uint32_t l : m_nodes[nNodeIdx]->vNbrNodes)
			{
				double fWeight = m_tess.get_cell(nNodeIdx)->pFaceSquare[l]
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