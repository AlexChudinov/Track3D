#pragma once
#ifndef _MESH_ELEMENT_
#define _MESH_ELEMENT_

#include <list>
#include <array>
#include <Elements.h>

/**
 * Possible new mesh element implementation [AC 16.01.2018]
 */

class MeshTet
{
public:
	typedef UINT Idx;
	typedef std::array<Idx, 4> IdxList;
	typedef std::list<MeshTet> Tets;
	typedef EvaporatingParticle::CNode3D Node;
	typedef std::vector<Node*> Nodes;
	typedef EvaporatingParticle::Vector3D Vector3D;
	typedef EvaporatingParticle::Matrix3D Matrix3D;
	typedef std::array<double, 4> DblList;
	typedef EvaporatingParticle::CElem3D Elem;

private:
	static const Nodes * s_globalNodes;
public:
	static void setGlobalNodes(const Nodes* nodes);

	MeshTet(const IdxList& nodeIdxs);

	/**
	 * Returns baricentric coordinates of a point p0
	 */
	DblList barycentricCoords(const Vector3D p0) const;

	/**
	 * Returns node indices
	 */
	const IdxList& nodeIdxs() const;

	/**
	* ANSYS elemet splitting into tets
	*/
	static Tets splitAnsysElement(const Elem& e);

private:
	IdxList m_nodeIdxs;

	/**
	 * Creates baricentric tensor
	 */
	Matrix3D barycentricTensor() const;

	/**
	 * Transformations from ANSYS elements
	 */
	static Tets splitAnsysTet(Idx n0, Idx n1, Idx n2, Idx n3);
	static Tets splitAnsysPyr(Idx n0, Idx n1, Idx n2, Idx n3, Idx n4);
	static Tets splitAnsysWedge(Idx n0, Idx n1, Idx n2, Idx n3, Idx n4, Idx n5);
	static Tets splitAnsysHex(Idx n0, Idx n1, Idx n2, Idx n3, Idx n4, Idx n5, Idx n6, Idx n7);
};

#endif