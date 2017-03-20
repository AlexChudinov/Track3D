#include "stdafx.h"
#include <queue>
#include "MeshData.h"

bool BoundaryMesh::empty() const
{
	return m_mapBoundariesList.empty(); 
}

uint32_t BoundaryMesh::maxLabel() const
{
	return m_mapReversedBoundariesList.rbegin()->first;
}

void BoundaryMesh::addBoundary(
	const std::string& strName, 
	const std::vector<uint32_t>& vLabels,
	const std::vector<Vector3D>& vNormals,
	BoundaryType type)
{
	if (vLabels.size() != vNormals.size())
		throw std::runtime_error("BoundaryMesh::addBoundary:"
			" Sizes of normals and labels vectors are different.");
	if (isBoundary(strName)) removeBoundary(strName);
	m_mapBoundariesList[strName] = std::make_pair(type, SetLabels(vLabels.begin(), vLabels.end()));
	BoundariesMap::const_iterator it = m_mapBoundariesList.find(strName);
	for (size_t i = 0; i < vLabels.size(); ++i)
	{
		std::pair<Vector3D, NamesList>& entry = m_mapReversedBoundariesList[vLabels[i]];
		entry.first = vNormals[i];
		entry.second.insert(std::cref(it->first));
	}
}

void BoundaryMesh::removeBoundary(const std::string& strName)
{
	if(m_pBoundaryVals) m_pBoundaryVals->removeBoundary(strName);
	for (uint32_t l : m_mapBoundariesList.at(strName).second)
	{
		m_mapReversedBoundariesList[l].second.erase(strName);
		if (m_mapReversedBoundariesList[l].second.empty()) m_mapReversedBoundariesList.erase(l);
	}
	m_mapBoundariesList.erase(strName);
}

void BoundaryMesh::boundaryType(const std::string& strName, BoundaryType type)
{
	m_mapBoundariesList.at(strName).first = type;
}

BoundaryMesh::BoundaryType BoundaryMesh::boundaryType(const std::string& strName) const
{
	return m_mapBoundariesList.at(strName).first;
}

const BoundaryMesh::NamesList& BoundaryMesh::boundaryNames(uint32_t l) const
{
	return m_mapReversedBoundariesList.at(l).second;
}

const BoundaryMesh::SetLabels& BoundaryMesh::boundaryLabels(const std::string& strName) const
{
	return m_mapBoundariesList.at(strName).second;
}

BoundaryMesh::const_iterator BoundaryMesh::begin()const{ return m_mapReversedBoundariesList.begin(); }
BoundaryMesh::iterator BoundaryMesh::begin() { return m_mapReversedBoundariesList.begin(); }
BoundaryMesh::const_iterator BoundaryMesh::end() const { return m_mapReversedBoundariesList.end(); }
BoundaryMesh::iterator BoundaryMesh::end() { return m_mapReversedBoundariesList.end(); }

bool BoundaryMesh::isBoundary(const std::string& sName) const
{
	return m_mapBoundariesList.find(sName) != m_mapBoundariesList.end();
}

bool BoundaryMesh::isBoundary(uint32_t l) const
{
	return m_mapReversedBoundariesList.find(l) != m_mapReversedBoundariesList.end();
}


bool BoundaryMesh::isFirstType(const std::string& strName) const
{
	return m_mapBoundariesList.at(strName).first == FIXED_VAL;
}

bool BoundaryMesh::isFirstType(uint32_t l) const
{
	return std::accumulate(m_mapReversedBoundariesList.at(l).second.begin(),
		m_mapReversedBoundariesList.at(l).second.end(), false,
		[=](bool val, const std::string& strName)->bool
	{
		return val |= isFirstType(strName);
	});
}

BoundaryMesh::Vector3D BoundaryMesh::normal(uint32_t l) const
{
	return m_mapReversedBoundariesList.at(l).first;
}

size_t BoundaryMesh::patchSize(const std::string & strName) const
{
	return m_mapBoundariesList.at(strName).second.size();
}

const CMeshAdapter::Element * CMeshAdapter::element(
	const CMeshAdapter::Vector3D & v,
	CMeshAdapter::Label & nCurNode,
	const CMeshAdapter::Label & nPrevNode) const
{
	std::vector<bool> visitedNodes(m_nodes.size(), false);
	std::set<Label> visitedElemets;
	std::queue<Label> nodesQueue;
	nodesQueue.push(nPrevNode);
	visitedNodes[nPrevNode] = true;

	while (!nodesQueue.empty())
	{
		nCurNode = nodesQueue.front(); nodesQueue.pop();
		const Labels& vNbrElemIdxs = m_nodes[nCurNode]->nbr_elems();
		for (UINT nElemIdx : vNbrElemIdxs)
			if (visitedElemets.find(nElemIdx) == visitedElemets.end())
				if (m_elems[nElemIdx]->inside(v)) return m_elems[nElemIdx];
				else
				{
					visitedElemets.insert(nElemIdx);
					const UINT * pNodeIdxStart = m_elems[nElemIdx]->nodes();
					const UINT * pNodeIdxEnd = pNodeIdxStart + m_elems[nElemIdx]->get_node_count();
					for (; pNodeIdxStart != pNodeIdxEnd; ++pNodeIdxStart)
						if (!visitedNodes[*pNodeIdxStart])
						{
							nodesQueue.push(*pNodeIdxStart);
							visitedNodes[*pNodeIdxStart] = true;
						}
				}
	}
	return nullptr;
}

CMeshAdapter::InterpCoefs & CMeshAdapter::add(InterpCoefs & ic1, const InterpCoefs & ic2)
{
	for (const InterpCoef& c : ic2)
	{
		InterpCoefs::iterator it = ic1.find(c.first);
		if (it == ic1.end())
			ic1[c.first] = c.second;
		else
			it->second += c.second;
	}
	return ic1;
}

CMeshAdapter::InterpCoefs & CMeshAdapter::mul(double h, InterpCoefs & ic)
{
	for (auto& c : ic) c.second *= h;
	return ic;
}

CMeshAdapter::Matrix3D CMeshAdapter::covarianceOfDirections(Label l) const
{
	Matrix3D result; result.fill(0.0);
	for (Label i : m_lazyGraph->neighbor(l))
	{
		Vector3D v = m_nodes[i]->pos - m_nodes[l]->pos;
		double fSqLength = v.sqlength(); fSqLength *= fSqLength;
		result(0, 0) += (v.x*v.x) / fSqLength;
		result(1, 1) += (v.y*v.y) / fSqLength;
		result(2, 2) += (v.z*v.z) / fSqLength;
		result(0, 1) = result(1, 0) += (v.x*v.y) / fSqLength;
		result(0, 2) = result(2, 0) += (v.x*v.z) / fSqLength;
		result(1, 2) = result(2, 1) += (v.y*v.z) / fSqLength;
	}
	return result;
}

CMeshAdapter::Vector3DOp CMeshAdapter::finDiffDirCov(Label l) const
{
	Vector3DOp result;
	for (Label i : m_lazyGraph->neighbor(l))
	{
		Vector3D v = m_nodes[i]->pos - m_nodes[l]->pos;
		double fSqLength = v.sqlength(); fSqLength *= fSqLength;
		add(result[0], mul(v.x / fSqLength, InterpCoefs{ { uint32_t(i), 1.0 },{ uint32_t(l), -1.0 } }));
		add(result[1], mul(v.y / fSqLength, InterpCoefs{ { uint32_t(i), 1.0 },{ uint32_t(l), -1.0 } }));
		add(result[2], mul(v.z / fSqLength, InterpCoefs{ { uint32_t(i), 1.0 },{ uint32_t(l), -1.0 } }));
	}
	return result;
}

CMeshAdapter::InterpCoefs CMeshAdapter::gradX(Label l) const
{
	InterpCoefs result;
	Matrix3D m = covarianceOfDirections(l);
	double fDet = math::det(m);
	if (fDet*fDet <= eps())
		throw std::runtime_error("CMeshAdapter::gradX: Determinant is zero!");
	Vector3DOp vDiff = finDiffDirCov(l);
	mul(1. / fDet,
		add(result,
			add(mul(math::det(Matrix2D{ { m(1,1), m(1,2) },{ m(2,1), m(2,2) } }), vDiff[0]),
				add(mul(-math::det(Matrix2D{ { m(0,1), m(0,2) },{ m(2,1), m(2,2) } }), vDiff[1]),
					mul(math::det(Matrix2D{ { m(0,1), m(0,2) },{ m(1,1), m(1,2) } }), vDiff[2])))));
	return result;
}

CMeshAdapter::InterpCoefs CMeshAdapter::gradY(Label l) const
{
	InterpCoefs result;
	Matrix3D m = covarianceOfDirections(l);
	double fDet = math::det(m);
	if (fDet*fDet <= eps())
		throw std::runtime_error("CMeshAdapter::gradY: Determinant is zero!");
	Vector3DOp vDiff = finDiffDirCov(l);
	mul(1. / fDet,
		add(result,
			add(mul(-math::det(Matrix2D{ { m(1,0), m(1,2) },{ m(2,0), m(2,2) } }), vDiff[0]),
				add(mul(math::det(Matrix2D{ { m(0,0), m(0,2) },{ m(2,0), m(2,2) } }), vDiff[1]),
					mul(-math::det(Matrix2D{ { m(0,0), m(0,2) },{ m(1,0), m(1,2) } }), vDiff[2])))));
	return result;
}

CMeshAdapter::InterpCoefs CMeshAdapter::gradZ(Label l) const
{
	InterpCoefs result;
	Matrix3D m = covarianceOfDirections(l);
	double fDet = math::det(m);
	if (fDet*fDet <= eps())
		throw std::runtime_error("CMeshAdapter::gradZ: Determinant is zero!");
	Vector3DOp vDiff = finDiffDirCov(l);
	mul(1. / fDet,
		add(result,
			add(mul(math::det(Matrix2D{ { m(1,0), m(1,1) },{ m(2,0), m(2,1) } }), vDiff[0]),
				add(mul(-math::det(Matrix2D{ { m(0,0), m(0,1) },{ m(2,0), m(2,1) } }), vDiff[1]),
					mul(math::det(Matrix2D{ { m(0,0), m(0,1) },{ m(1,0), m(1,1) } }), vDiff[2])))));
	return result;
}

CMeshAdapter::ScalarFieldOperator CMeshAdapter::gradX()
{
	lazyGraphCreation();
	m_pProgressBar->set_job_name("Creating gradient x-component operator...");
	m_pProgressBar->set_progress(0);

	ScalarFieldOperator result; 
	result.m_matrix.resize(m_nodes.size());

	ThreadPool::getInstance().splitInPar(m_nodes.size(),
		[&](size_t nNodeIdx)->void { result.m_matrix[nNodeIdx] = gradX(nNodeIdx); },
		m_pProgressBar.get());

	return result;
}

CMeshAdapter::ScalarFieldOperator CMeshAdapter::gradY()
{
	lazyGraphCreation();
	m_pProgressBar->set_job_name("Creating gradient y-component operator...");
	m_pProgressBar->set_progress(0);

	ScalarFieldOperator result;
	result.m_matrix.resize(m_nodes.size());

	ThreadPool::getInstance().splitInPar(m_nodes.size(),
		[&](size_t nNodeIdx)->void { result.m_matrix[nNodeIdx] = gradY(nNodeIdx); },
		m_pProgressBar.get());

	return result;
}

CMeshAdapter::ScalarFieldOperator CMeshAdapter::gradZ()
{
	lazyGraphCreation();
	m_pProgressBar->set_job_name("Creating gradient z-component operator...");
	m_pProgressBar->set_progress(0);

	ScalarFieldOperator result;
	result.m_matrix.resize(m_nodes.size());

	ThreadPool::getInstance().splitInPar(m_nodes.size(),
		[&](size_t nNodeIdx)->void { result.m_matrix[nNodeIdx] = gradZ(nNodeIdx); },
		m_pProgressBar.get());

	return result;
}

void CMeshAdapter::lazyGraphCreation()
{
	if (!m_lazyGraph)
	{
		m_pProgressBar->set_job_name("Creating connectivity graph...");
		m_pProgressBar->set_progress(0);
		m_lazyGraph.reset(new Graph);
		for (const Element* e : m_elems)
		{
			if (e->nInd % 1000 == 0) m_pProgressBar->set_progress(e->nInd * 100 / m_elems.size());
			if (m_pProgressBar->get_terminate_flag()) return m_lazyGraph.reset();
			m_lazyGraph->addElem(e);
		}
	}
}

std::vector<CMeshAdapter::NodeType> CMeshAdapter::nodeTypes() const
{
	std::vector<NodeType> res(m_nodes.size());
	ThreadPool::getInstance().splitInPar(m_nodes.size(),
		[&](size_t nNodeIdx) 
	{
		res[nNodeIdx] = boundaryMesh()->isBoundary(nNodeIdx) ?
			(boundaryMesh()->isFirstType(nNodeIdx) ? FirstTypeBoundaryNode : SecondTypeBoundaryNode)
			: InnerNode;
	});
	return res;
}

CMeshAdapter::ScalarFieldOperator CMeshAdapter::laplacianSolver0()
{
	ScalarFieldOperator result;
	lazyGraphCreation();

	m_pProgressBar->set_job_name("Laplacian solver #0 creation...");
	m_pProgressBar->set_progress(0);

	result.m_matrix.resize(m_nodes.size());
	for (const Node* n : m_nodes)
	{
		Label nNodeIdx = n->nInd;
		if (nNodeIdx % 1000 == 0) m_pProgressBar->set_progress(nNodeIdx * 100 / m_nodes.size());
		if (m_pProgressBar->get_terminate_flag()) break;
		if (m_pBoundary->isBoundary(nNodeIdx))
		{
			if (m_pBoundary->isFirstType(nNodeIdx)) //Fixed value BC
				result.m_matrix[nNodeIdx][nNodeIdx] = 1.0;
			else
			{
				double fTotDist = 0.0, fDist;
				InterpCoefs coefs;
				for (Label l : m_lazyGraph->neighbor(nNodeIdx))
				{
					fDist = m_pBoundary->isBoundary(l) ?
						1. / (m_nodes[l]->pos - n->pos).length()
						: 2. / (m_nodes[l]->pos - n->pos).length();
					fTotDist += fDist;
					add(coefs, { InterpCoef(uint32_t(l), fDist) });
				}
				mul(1. / fTotDist, coefs);
				result.m_matrix[nNodeIdx] = std::move(coefs);
			}
		}
		else
		{
			double fTotDist = 0.0, fDist;
			InterpCoefs coefs;
			for (Label l : m_lazyGraph->neighbor(nNodeIdx))
			{
				fDist = 1. / (m_nodes[l]->pos - n->pos).sqlength();
				fTotDist += fDist;
				add(coefs, { InterpCoef(uint32_t(l), fDist) });
			}
			mul(1. / fTotDist, coefs);
			result.m_matrix[nNodeIdx] = std::move(coefs);
		}
	}
	return result;
}

CMeshAdapter::ScalarFieldOperator CMeshAdapter::laplacianSolver1()
{
	ScalarFieldOperator result;
	lazyGraphCreation();

	m_pProgressBar->set_job_name("Laplacian solver #1 creation...");
	m_pProgressBar->set_progress(0);

	result.m_matrix.resize(m_nodes.size());

	std::vector<NodeType> vNodeTypes(nodeTypes());

	ThreadPool::getInstance().splitInPar(m_nodes.size(),
		[&](size_t nNodeIdx) 
	{
		Label nCurNodeIdx = static_cast<Label>(nNodeIdx), nNextNodeIdx;
		const Node* n = m_nodes[nCurNodeIdx];
		switch (vNodeTypes[nCurNodeIdx])
		{
		case FirstTypeBoundaryNode: //Fixed value BC
			result.m_matrix[nCurNodeIdx][nCurNodeIdx] = 1.0;
			break;
		case SecondTypeBoundaryNode: //Zero gradient
		{
			Vector3D norm = boundaryMesh()->normal(nCurNodeIdx);
			InterpCoefs normalDerivative = std::move(add(
				mul(norm.x, gradX(nCurNodeIdx)),
				add(mul(norm.y, gradY(nCurNodeIdx)),
					mul(norm.z, gradZ(nCurNodeIdx)))));
			double fSum = normalDerivative[nCurNodeIdx];
			normalDerivative.erase(nCurNodeIdx);
			mul(-1. / fSum, normalDerivative);
			result.m_matrix[nCurNodeIdx] = std::move(normalDerivative);
			break;
		}
		default: //It is inner point
		{
			double
				hx1 = optimalStep({ -1.0, 0.0, 0.0 }, nCurNodeIdx),
				hx2 = optimalStep({ 1.0, 0.0, 0.0 }, nCurNodeIdx),
				hy1 = optimalStep({ 0.0, -1.0, 0.0 }, nCurNodeIdx),
				hy2 = optimalStep({ 0.0, 1.0, 0.0 }, nCurNodeIdx),
				hz1 = optimalStep({ 0.0, 0.0, -1.0 }, nCurNodeIdx),
				hz2 = optimalStep({ 0.0, 0.0, 1.0 }, nCurNodeIdx);

			Vector3D
				x0{ n->pos.x - hx1, n->pos.y, n->pos.z },
				x1{ n->pos.x + hx2, n->pos.y, n->pos.z },
				y0{ n->pos.x, n->pos.y - hy1, n->pos.z },
				y1{ n->pos.x, n->pos.y + hy2, n->pos.z },
				z0{ n->pos.x, n->pos.y, n->pos.z - hz1 },
				z1{ n->pos.x, n->pos.y, n->pos.z + hz2 };

			InterpCoefs coefsX = mul(1. / hx1, interpCoefs(x0, nNextNodeIdx, nCurNodeIdx));
			add(coefsX, mul(1. / hx2, interpCoefs(x1, nNextNodeIdx, nCurNodeIdx)));
			mul(1. / (hx2 + hx1), coefsX);

			InterpCoefs coefsY = mul(1. / hy1, interpCoefs(y0, nNextNodeIdx, nCurNodeIdx));
			add(coefsY, mul(1. / hy2, interpCoefs(y1, nNextNodeIdx, nCurNodeIdx)));
			mul(1. / (hy2 + hy1), coefsY);

			InterpCoefs coefsZ = mul(1. / hz1, interpCoefs(z0, nNextNodeIdx, nCurNodeIdx));
			add(coefsZ, mul(1. / hz2, interpCoefs(z1, nNextNodeIdx, nCurNodeIdx)));
			mul(1. / (hz2 + hz1), coefsZ);

			mul(1. / (
				(1. / hx1 + 1. / hx2) / (hx1 + hx2)
				+ (1. / hy1 + 1. / hy2) / (hy1 + hy2)
				+ (1. / hz1 + 1. / hz2) / (hz1 + hz2)),
				add(coefsX, add(coefsY, coefsZ)));
			result.m_matrix[nCurNodeIdx] = std::move(coefsX);
		}
		}
	}, m_pProgressBar.get());

	return result;
}

CMeshAdapter::ScalarFieldOperator CMeshAdapter::directedDerivative(const Vector3D & dir)
{
	ScalarFieldOperator result;
	lazyGraphCreation();

	std::string strJobName = std::string("Derivative in direction {"
		+ std::to_string(dir.x) + ", " + std::to_string(dir.y) + ", " + std::to_string(dir.z)
		+ "} creation...");
	m_pProgressBar->set_job_name(strJobName.c_str());
	m_pProgressBar->set_progress(0);

	result.m_matrix.resize(m_nodes.size());

	for (const Node* n : m_nodes)
	{
		size_t nNodeIdxNext;
		if (n->nInd % 1000 == 0) m_pProgressBar->set_progress(n->nInd * 100 / m_nodes.size());
		if (m_pProgressBar->get_terminate_flag()) break;
		
		add(result.m_matrix[n->nInd],
			add(mul(dir.x, gradX(n->nInd)),
				add(mul(dir.y, gradY(n->nInd)),
					mul(dir.z, gradZ(n->nInd)))));
	}

	return result;
}

CMeshAdapter::InterpCoefs CMeshAdapter::interpCoefs(const Vector3D & pos, const Element * e) const
{
	InterpCoefs res;
	double s, t, u;
	if (!e->param(pos, s, t, u))
		throw std::runtime_error("CMeshAdapter::interpCoefs"
			": CElem3D::param returned false.");
	for (size_t i = 0; i < e->get_node_count(); ++i)
	{
		double w = e->shape_func(s, t, u, i);
		res.insert(InterpCoef(e->get_node_index(i), w));
	}
	return res;
}

const CMeshAdapter::Element * CMeshAdapter::lookInClosestElements(const Vector3D & pos, Label l) const
{
	const Labels& vElemIdxs = m_nodes[l]->nbr_elems();
	for (Label nElemIdx : vElemIdxs)
		if (m_elems[nElemIdx]->inside(pos)) return m_elems[nElemIdx];
	return nullptr;
}

CMeshAdapter::CMeshAdapter(const Elements & es, const Nodes & ns, double fSmallStepFactor)
	:
	m_elems(es),
	m_nodes(ns),
	m_pBoundary(new BoundaryMesh),
	m_pProgressBar(new ProgressBar),
	m_fSmallStepFactor(fSmallStepFactor),
	m_fEpsilon(100*std::numeric_limits<double>::epsilon())
{
}

CMeshAdapter::ProgressBar * CMeshAdapter::progressBar() const
{
	return m_pProgressBar.get();
}

void CMeshAdapter::releaseGraph()
{
	m_lazyGraph.reset();
}

double CMeshAdapter::smallStepFactor() const
{
	return m_fSmallStepFactor;
}

void CMeshAdapter::smallStepFactor(double fVal)
{
	m_fSmallStepFactor = fVal;
}

double CMeshAdapter::eps() const
{
	return m_fEpsilon;
}

void CMeshAdapter::eps(size_t nFactor)
{
	m_fEpsilon = std::numeric_limits<double>::epsilon()*nFactor;
}

double CMeshAdapter::minElemSize(Label l) const
{
	const Labels& vElemIdxs = m_nodes[l]->nbr_elems();
	Labels::const_iterator pMinElem = std::min_element(vElemIdxs.begin(), vElemIdxs.end(),
		[=](Label l1, Label l2)->bool
	{
		Vector3D
			vBox1 = m_elems[l1]->box.vMax - m_elems[l1]->box.vMin,
			vBox2 = m_elems[l2]->box.vMax - m_elems[l2]->box.vMin;
		double
			fBox1 = min(vBox1.x, min(vBox1.y, vBox1.z)),
			fBox2 = min(vBox2.x, min(vBox2.y, vBox2.z));
		return fBox1 < fBox2;
	});
	Vector3D vMinBox = m_elems[*pMinElem]->box.vMax - m_elems[*pMinElem]->box.vMin;
	return min(vMinBox.x, min(vMinBox.y, vMinBox.z));
}

double CMeshAdapter::minEdgeLength(Label l) const
{
	const CMeshConnectivity::NodeConnections& vNeighborLabels = m_lazyGraph->neighbor(l);
	return (m_nodes[*std::min_element(vNeighborLabels.begin(), vNeighborLabels.end(),
		[=](Label l1, Label l2)->bool
	{
		Vector3D 
			df1 = m_nodes[l1]->pos - m_nodes[l]->pos,
			df2 = m_nodes[l2]->pos - m_nodes[l]->pos;
		return df1.sqlength() < df2.sqlength();
	})]->pos - m_nodes[l]->pos).length();
}

double CMeshAdapter::maxEdgeLength(Label l) const
{
	const CMeshConnectivity::NodeConnections& vNeighborLabels = m_lazyGraph->neighbor(l);
	return (m_nodes[*std::max_element(vNeighborLabels.begin(), vNeighborLabels.end(),
		[=](Label l1, Label l2)->bool
	{
		Vector3D
			df1 = m_nodes[l1]->pos - m_nodes[l]->pos,
			df2 = m_nodes[l2]->pos - m_nodes[l]->pos;
		return df1.sqlength() < df2.sqlength();
	})]->pos - m_nodes[l]->pos).length();
}

double CMeshAdapter::optimalStep(const Vector3D& dir, Label l, Label deep) const
{
	double a = 0.0, b = 2.0*maxEdgeLength(l);
	while (lookInNeighbor(m_nodes[l]->pos + dir*b, l, deep)) b *= 2.;
	while (b - a > eps())
	{
		double m = (b + a) / 2.;
		if (lookInNeighbor(m_nodes[l]->pos + dir*m, l, deep))
			a = m;
		else
			b = m;
	}
	return a;
}

const CMeshAdapter::Element * CMeshAdapter::lookInNeighbor(const Vector3D & pos, Label l, Label deep) const
{
	const Element* result = lookInClosestElements(pos, l);
	if (result || deep == 0) return result;

	std::set<Label> 
		visitedNodes, 
		visitedElemets(m_nodes[l]->nbr_elems().begin(), m_nodes[l]->nbr_elems().end());
	visitedNodes.insert(l);
	
	for(Label nCurDeep = 0; nCurDeep < deep; ++nCurDeep)
	{
		std::set<Label> nextElemsToVisit;
		for (Label nElemIdx : visitedElemets)
		{
			if (visitedNodes.size() == m_nodes.size()) return nullptr;
			const Label
				*pFirstNodeIdx = m_elems[nElemIdx]->nodes(),
				*pLastNodeIdx = pFirstNodeIdx + m_elems[nElemIdx]->get_node_count();
			for (; pFirstNodeIdx < pLastNodeIdx; ++pFirstNodeIdx)
				if (visitedNodes.find(*pFirstNodeIdx) == visitedNodes.end())
				{
					visitedNodes.insert(*pFirstNodeIdx);
					for (Label nElemIdx : m_nodes[*pFirstNodeIdx]->nbr_elems())
						if (nextElemsToVisit.find(nElemIdx) == nextElemsToVisit.end()
							&& visitedElemets.find(nElemIdx) == visitedElemets.end())
							if (m_elems[nElemIdx]->inside(pos)) return m_elems[nElemIdx];
							else nextElemsToVisit.insert(nElemIdx);
				}
		}
		visitedElemets = std::move(nextElemsToVisit);
	}
	return nullptr;
}

BoundaryMesh * CMeshAdapter::boundaryMesh() const
{
	return m_pBoundary.get();
}

CMeshAdapter::InterpCoefs CMeshAdapter::interpCoefs(
	const Vector3D& v,
	Label& nCurNode,
	const Label& nPrevNode
) const
{
	const Element * e = element(v, nCurNode, nPrevNode);
	if (e) return interpCoefs(v, e);
	return InterpCoefs();
}

CMeshAdapter::PScalFieldOp CMeshAdapter::createOperator(ScalarOperatorType type)
{
	switch (type)
	{
	case CMeshAdapter::LaplacianSolver0:
		return PScalFieldOp(new ScalarFieldOperator(laplacianSolver0()));
	case CMeshAdapter::LaplacianSolver1:
		return PScalFieldOp(new ScalarFieldOperator(laplacianSolver1()));
	case CMeshAdapter::GradX:
		return PScalFieldOp(new ScalarFieldOperator(gradX()));
	case CMeshAdapter::GradY:
		return PScalFieldOp(new ScalarFieldOperator(gradY()));
	case CMeshAdapter::GradZ:
		return PScalFieldOp(new ScalarFieldOperator(gradZ()));
	default:
		throw std::runtime_error("CMeshAdapter::createOperator: Unsupported operator type.");
	}
}

void CMeshConnectivity::addTet(Label n0, Label n1, Label n2, Label n3)
{
	addEdge(n0, n1);
	addEdge(n0, n2);
	addEdge(n0, n3);
	addEdge(n1, n2);
	addEdge(n1, n3);
	addEdge(n2, n3);
}

void CMeshConnectivity::addPyr(Label n0, Label n1, Label n2, Label n3, Label n4)
{
	addEdge(n0, n1);
	addEdge(n0, n3);
	addEdge(n0, n4);
	addEdge(n1, n2);
	addEdge(n1, n4);
	addEdge(n2, n3);
	addEdge(n2, n4);
	addEdge(n3, n4);
}

void CMeshConnectivity::addWedge(Label n0, Label n1, Label n2, Label n3, Label n4, Label n5)
{
	addEdge(n0, n1);
	addEdge(n0, n2);
	addEdge(n0, n3);
	addEdge(n1, n2);
	addEdge(n1, n4);
	addEdge(n2, n5);
	addEdge(n3, n4);
	addEdge(n3, n5);
	addEdge(n4, n5);
}

void CMeshConnectivity::addHexa(Label n0, Label n1, Label n2, Label n3, Label n4, Label n5, Label n6, Label n7)
{
	addEdge(n0, n1);
	addEdge(n0, n3);
	addEdge(n0, n4);
	addEdge(n1, n2);
	addEdge(n1, n5);
	addEdge(n2, n3);
	addEdge(n2, n6);
	addEdge(n3, n7);
	addEdge(n4, n5);
	addEdge(n4, n7);
	addEdge(n5, n6);
	addEdge(n6, n7);
}

CMeshConnectivity::Label CMeshConnectivity::size() const
{
	return m_graph.size();
}

void CMeshConnectivity::addEdge(Label i, Label j)
{
	Label max = max(i, j);
	Label min = min(i, j);
	if (max >= m_graph.size()) m_graph.resize(max + 1);
	m_graph[max].insert(min);
	m_graph[min].insert(max);
}

void CMeshConnectivity::addElem(const Elem * e)
{
	const UINT * pNodeIdx = e->nodes();
	switch (e->get_node_count())
	{
	case 4: //Tetrahedron
		return addTet(pNodeIdx[0], pNodeIdx[1], pNodeIdx[2], pNodeIdx[3]);
	case 5: //Pyramide
		return addPyr(pNodeIdx[0], pNodeIdx[1], pNodeIdx[2], pNodeIdx[3], pNodeIdx[4]);
	case 6: //Wedge
		return addWedge(pNodeIdx[0], pNodeIdx[1], pNodeIdx[2], pNodeIdx[3], pNodeIdx[4], pNodeIdx[5]);
	case 8: //Hexahedron
		return addHexa(pNodeIdx[0], pNodeIdx[1], pNodeIdx[2], pNodeIdx[3], 
			pNodeIdx[4], pNodeIdx[5], pNodeIdx[6], pNodeIdx[7]);
	default:
		throw std::runtime_error(std::string("CMeshConnectivity::addElem: ")
			+ "Unexpected number of vertices - " + std::to_string(e->get_node_count())
			+ " in element #" + std::to_string(e->nInd) + ".");
	}
}

const CMeshConnectivity::NodeConnections& CMeshConnectivity::neighbor(Label i) const
{
	return m_graph.at(i);
}
