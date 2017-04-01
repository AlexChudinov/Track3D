#pragma once
#ifndef _BOUNDARY_MESH_
#define _BOUNDARY_MESH_

#include <map>
#include <set>
#include <memory>
#include <vector>
#include <numeric>
#include <algorithm>

#include "../track3d/vector3d.hpp"
#include "../utilities/MemoryPool.h"

//Boundary conditions for the mesh
class BoundaryMesh
{
	using StrRef = std::reference_wrapper<const std::string>;

	//Declare reference wrappers less operator
	class StrRefLess : public std::binary_function<const StrRef&, const StrRef&, bool>
	{
	public: bool operator()(const StrRef& s1, const StrRef& s2) const 
	{ 
		return s1.get() < s2.get(); 
	}
	};

	using Label = uint32_t;
	using BoundaryValsEntry = std::map<Label, double, std::less<Label>,
		Allocator<std::pair<Label, double>>>;
	using BoundaryValsTab = std::map<std::string, BoundaryValsEntry, std::less<std::string>,
		Allocator<std::pair<std::string, BoundaryValsEntry>>>;

public:
	enum BoundaryType { ZERO_GRAD, FIXED_VAL };

	using Vector3D = EvaporatingParticle::Vector3<double>;
	using SetLabels = std::set<Label, std::less<Label>, Allocator<Label>>;
	using BoundaryDescription = std::pair<BoundaryType, SetLabels>;
	using BoundariesMap = std::map<std::string, BoundaryDescription, std::less<std::string>,
		Allocator<std::pair<std::string, BoundaryDescription>>>;
	using NamesList = std::set<StrRef, StrRefLess, Allocator<StrRef>>;
	using ReversedBoundariesMap = std::map<Label, std::pair<Vector3D, NamesList>, std::less<Label>,
		Allocator<std::pair<Label, std::pair<Vector3D, NamesList>>>>;
	using BoundaryNormals = std::map<Label, Vector3D, std::less<Label>,
		Allocator<std::pair<Label, Vector3D>>>;

	using iterator = ReversedBoundariesMap::iterator;
	using const_iterator = ReversedBoundariesMap::const_iterator;

private:
	BoundariesMap m_mapBoundariesList;
	ReversedBoundariesMap m_mapReversedBoundariesList;
	BoundaryValsTab m_BoundaryVals;

public:
	//Checks if boundary is empty
	bool empty() const;
	//Returns biggest number of boundary
	Label maxLabel() const;

	//Adds new boundary patch
	void addBoundary(
		const std::string& strName,
		const std::vector<uint32_t>& vLabels,
		const std::vector<Vector3D>& vNormals,
		BoundaryType type = FIXED_VAL);

	//Removes existing boundary patch
	void removeBoundary(const std::string& strName);

	//Sets type of a boundary with a name strName
	void boundaryType(const std::string& strName, BoundaryType type);

	//Returns type of a boundary with a name strName
	BoundaryType boundaryType(const std::string& strName) const;

	//Returns a set of boundaries connected to a given label
	const NamesList& boundaryNames(Label l) const;

	//Returns numbers of nodes which are belong to a given boundary
	const SetLabels& boundaryLabels(const std::string& strName) const;

	//Get iterators for boundary
	const_iterator begin() const;
	iterator begin();
	const_iterator end() const;
	iterator end();

	//Checks if the name is in the boundaries list
	bool isBoundary(const std::string& sName) const;
	//Checks if the label belongs to the boundary
	bool isBoundary(uint32_t l) const;

	//Checks if it is a first-type (Dirichlet) boundary condition
	bool isFirstType(const std::string& sName) const;
	bool isFirstType(Label l) const;

	//Gets the normal for given node label
	const Vector3D& normal(Label l) const;

	//Returns size of a patch strName
	size_t patchSize(const std::string& strName) const;

	//Sets boundary vals
	void boundaryVals(const std::string& strName, const std::vector<double>& vVals);

	void applyBoundaryVals(std::vector<double>& f) const;
};

#endif // !_BOUNDARY_MESH_
