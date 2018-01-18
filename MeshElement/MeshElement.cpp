#include "stdafx.h"
#include "MeshElement.h"

const MeshTet::Nodes* MeshTet::s_globalNodes;

void MeshTet::setGlobalNodes(const Nodes* nodes)
{
	s_globalNodes = nodes;
}

MeshTet::MeshTet(const IdxList & nodeIdxs)
	:
	m_nodeIdxs(nodeIdxs)
{
	if (!s_globalNodes) throw std::runtime_error(
		"MeshTet::MeshTet :"
		"Global nodes pointer was not initialised.");
	for (Idx i : m_nodeIdxs)
	{
		if (i >= s_globalNodes->size()) throw std::runtime_error(
			"MeshTet::MeshTet :"
			"There is no node with number: " + std::to_string(i) + ".");
	}
}

MeshTet::DblList MeshTet::barycentricCoords(const Vector3D p0) const
{
	const Nodes& nodes = *s_globalNodes;
	Matrix3D T = barycentricTensor();
	Vector3D r = p0 - nodes[m_nodeIdxs[0]]->pos;
	
	double TDet = T.det();
	Matrix3D rT = T; rT.set_col(r, 0);
	double l1 = rT.det() / TDet;

	rT.set_col(T.get_col(0), 0); rT.set_col(r, 1);
	double l2 = rT.det() / TDet;

	rT.set_col(T.get_col(1), 1); rT.set_col(r, 2);
	double l3 = rT.det() / TDet;

	return { 1.0 - l1 - l2 - l3, l1, l2, l3 };
}

const MeshTet::IdxList & MeshTet::nodeIdxs() const
{
	return m_nodeIdxs;
}

MeshTet::Tets MeshTet::splitAnsysElement(const Elem & e)
{
	const Idx* pNodeIdxs;

	if (!(pNodeIdxs = e.nodes()))
		throw std::runtime_error("MeshTet::splitAnsysElement : Null pointer for node indices.");

	switch (e.get_node_count())
	{
	case 4: return splitAnsysTet(pNodeIdxs[0], pNodeIdxs[1], pNodeIdxs[2], pNodeIdxs[3]);
	case 5: return splitAnsysPyr(pNodeIdxs[0], pNodeIdxs[1], pNodeIdxs[2], pNodeIdxs[3], pNodeIdxs[4]);
	case 6: return splitAnsysWedge(pNodeIdxs[0], pNodeIdxs[1], pNodeIdxs[2], pNodeIdxs[3], pNodeIdxs[4], pNodeIdxs[5]);
	case 8: return splitAnsysHex(pNodeIdxs[0], pNodeIdxs[1], pNodeIdxs[2], pNodeIdxs[3], pNodeIdxs[4], pNodeIdxs[5], pNodeIdxs[6], pNodeIdxs[7]);
	default: throw std::runtime_error("MeshTet::splitAnsysElement : "
		"Unexpected number of nodes: " + std::to_string(e.get_node_count()) + ".");
	}
}

MeshTet::Tets MeshTet::splitAnsysTet(Idx n0, Idx n1, Idx n2, Idx n3)
{
	return { MeshTet({n0, n1, n2, n3}) };
}

MeshTet::Tets MeshTet::splitAnsysPyr(Idx n0, Idx n1, Idx n2, Idx n3, Idx n4)
{
	Tets
		tet1 = splitAnsysTet(n0, n1, n3, n4),
		tet2 = splitAnsysTet(n1, n2, n3, n4);
	tet1.insert(tet1.end(), tet2.begin(), tet2.end());
	return tet1;
}

MeshTet::Tets MeshTet::splitAnsysWedge(Idx n0, Idx n1, Idx n2, Idx n3, Idx n4, Idx n5)
{
	Tets
		pyr = splitAnsysPyr(n0, n3, n4, n1, n5),
		tet = splitAnsysTet(n0, n1, n2, n5);
	pyr.insert(pyr.end(), tet.begin(), tet.end());
	return pyr;
}

MeshTet::Tets MeshTet::splitAnsysHex(Idx n0, Idx n1, Idx n2, Idx n3, Idx n4, Idx n5, Idx n6, Idx n7)
{
	Tets
		Wedge1 = splitAnsysWedge(n0, n1, n3, n4, n5, n7),
		Wedge2 = splitAnsysWedge(n1, n2, n3, n5, n6, n7);
	Wedge1.insert(Wedge1.end(), Wedge2.begin(), Wedge2.end());
	return Wedge1;
}

MeshTet::Matrix3D MeshTet::barycentricTensor() const
{
	const Nodes& nodes = *s_globalNodes;
	return Matrix3D(
		nodes[m_nodeIdxs[1]]->pos - nodes[m_nodeIdxs[0]]->pos,
		nodes[m_nodeIdxs[2]]->pos - nodes[m_nodeIdxs[0]]->pos,
		nodes[m_nodeIdxs[3]]->pos - nodes[m_nodeIdxs[0]]->pos);
}
