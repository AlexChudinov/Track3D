#include <stdafx.h>
#include "BoundaryMesh.h"

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
	if (m_BoundaryVals.find(strName) != m_BoundaryVals.end()) 
		m_BoundaryVals.erase(strName);
	for (Label l : m_mapBoundariesList.at(strName).second)
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

const BoundaryMesh::NamesList& BoundaryMesh::boundaryNames(Label l) const
{
	return m_mapReversedBoundariesList.at(l).second;
}

const BoundaryMesh::SetLabels& BoundaryMesh::boundaryLabels(const std::string& strName) const
{
	return m_mapBoundariesList.at(strName).second;
}

BoundaryMesh::const_iterator BoundaryMesh::begin()const { return m_mapReversedBoundariesList.begin(); }
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

bool BoundaryMesh::isFirstType(Label l) const
{
	return std::accumulate(m_mapReversedBoundariesList.at(l).second.begin(),
		m_mapReversedBoundariesList.at(l).second.end(), false,
		[=](bool val, const std::string& strName)->bool
	{
		return val |= isFirstType(strName);
	});
}

const BoundaryMesh::Vector3D& BoundaryMesh::normal(Label l) const
{
	return m_mapReversedBoundariesList.at(l).first;
}

size_t BoundaryMesh::patchSize(const std::string & strName) const
{
	return m_mapBoundariesList.at(strName).second.size();
}

void BoundaryMesh::boundaryVals(const std::string & strName, const std::vector<double>& vVals)
{
	if (!isBoundary(strName))
		throw std::runtime_error("BoundaryMesh::addBoundaryVals:"
			" Boundary " + strName + " was not found.");
	const SetLabels& vBoundaryLabels = boundaryLabels(strName);
	if (vVals.size() != vBoundaryLabels.size())
		throw std::runtime_error("BoundaryMesh::addBoundaryVals:"
			" Input vector of values and vector of boundary labels have different sizes.");
	SetLabels::const_iterator itLabels = vBoundaryLabels.begin();
	for (size_t i = 0; i < vVals.size(); ++i)
		m_BoundaryVals[strName][*(itLabels++)] = vVals[i];
}

void BoundaryMesh::applyBoundaryVals(std::vector<double>& f) const
{
	if (maxLabel() >= f.size())
		throw std::runtime_error("BoundaryMesh::applyBoundaryVals: "
			"Field size is too small.");

	for (const auto& boundaryPatch : *this)
	{
		int nPrimCondCount = 0;
		double fPrimCondAcc = 0.0;
		const NamesList& listNames = boundaryPatch.second.second;

		for (const auto& name : listNames)
			switch (boundaryType(name))
			{
			case ZERO_GRAD: break;
			case FIXED_VAL:
				fPrimCondAcc += m_BoundaryVals.at(name).at(boundaryPatch.first);
				nPrimCondCount++;
				break;
			default:
				throw std::runtime_error("BoundaryMesh::applyBoundaryVals:"
					" Unexpected boundary condition type.");
			}

		if (nPrimCondCount != 0)
			f[boundaryPatch.first] = fPrimCondAcc / nPrimCondCount;
		else
			f[boundaryPatch.first] = 0;
	}
}