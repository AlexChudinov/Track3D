#pragma once
#ifndef _FIELD_H_
#define _FIELD_H_

#ifdef _DEBUG
	#include <iostream>
#endif // _DEBUG

#include <map>
#include <memory>
#include <limits>
#include <linearAlgebra\matrixTemplate.h>

#include "mesh_geometry.h"

/**
* Field manipulation class
*/
template<typename field_type>
//Field type could be scalar or vector
class field
{
	template <typename field_type> friend class FieldLinearOp;
public:
	using mesh_geom = mesh_geometry<double, uint32_t>;
	using data_vector = std::vector<field_type>;
	using node_types_list = std::vector<bool>; //true if it is inner point and false if it is a boundary
	using node_labels_list = std::set<uint32_t>;
	using vector3f = math::vector_c<double, 3>;

	using BoundaryMesh = typename mesh_geom::BoundaryMesh;
	using BoundaryMeshSharedPtr = std::shared_ptr<BoundaryMesh>;
	using MeshSharedPtr = std::shared_ptr<mesh_geom>;
	using BoundaryValues = std::map<std::string, std::map<uint32_t, field_type>>;
private:
	//Keep reference to a space mesh
	MeshSharedPtr m_pMeshGeometry;
	BoundaryMeshSharedPtr m_pBoundaryMesh;

	data_vector _data; //Field data itself
	node_types_list _node_types; //Types of a field nodes, true if it is inner point and false if it is a boundary
	BoundaryValues m_boundaryFieldVals;
public:
	/**
	 * Creates zero filled field
	 */
	field(const MeshSharedPtr& meshGeometry)
		: 
		m_pMeshGeometry(meshGeometry),
		m_pBoundaryMesh(new BoundaryMesh(meshGeometry->createBoundary())),
		_data(m_pMeshGeometry->size(), field_type(0.0)),
		_node_types(m_pMeshGeometry->size(), true)
	{}

	const data_vector& data() const { return _data; }
	data_vector& data() { return _data; }

	//Returns field data size
	size_t size() const { return _data.size(); }

	/**
	 * Adds new boundary to a field
	 */
	void add_boundary(
		const std::string& sName, 
		const std::vector<uint32_t>& vLabels,
		const std::vector<vector3f>& vNormals
	)
	{
		m_pBoundaryMesh->addBoundary(sName, vLabels, vNormals);
		std::map<uint32_t, field_type>& boundaryPatch = m_boundaryFieldVals[sName];
		for (uint32_t l : vLabels)
		{
			boundaryPatch[l] = field_type(0.0);
			_node_types[l] = false;
		}
	}

	/**
	 * Sets values on a boundary
	 */
	void set_boundary_uniform_val(const std::string& sName, const field_type& val)
	{
		auto& boundaryPatch = m_boundaryFieldVals.at(sName);
		for (auto& boundaryNode : boundaryPatch) boundaryNode.second = val;
	}

	void set_boundary_vals(const std::string& sName, const data_vector& vals)
	{
		auto& boundaryPatch = m_boundaryFieldVals.at(sName);
		if (vals.size() != boundaryPatch.size())
			throw std::runtime_error("Boundary and input vector sizes mismatch.\n");
		size_t i = 0;
		for (auto& boundaryNode : boundaryPatch) boundaryNode.second = vals[i++];
	}

	/**
	 * Sets a boundary type
	 */
	void set_boundary_type(const std::string& sName, BoundaryMesh::BoundaryType type)
	{
		m_pBoundaryMesh->boundaryType(sName, type);
	}

	//Applies boundary conditions to a mesh
	//Puts averaged fixed values at FIXED_VAL boundary conditions and initializes ZERO_GRAD with zeros
	void applyBoundaryConditions()
	{
		for (const auto& boundaryLabel : *m_pBoundaryMesh)
		{
			int primaryCondition = 0;
			field_type primaryCondAcc = 0.0;
			const typename BoundaryMesh::NamesList& listNames = boundaryLabel.second.second;
			for (const auto& name : listNames)
			{
				switch (m_pBoundaryMesh->boundaryType(name))
				{
				case BoundaryMesh::ZERO_GRAD:
					break;
				case BoundaryMesh::FIXED_VAL:
					primaryCondAcc += m_boundaryFieldVals[name][boundaryLabel.first];
					primaryCondition++;
					break;
				default:
					throw std::runtime_error("Field::applyBoundaryConditions : Unexpected boundary condition type.");
				}
			}
			if (primaryCondition != 0)
				_data[boundaryLabel.first] = primaryCondAcc / primaryCondition;
			else
				_data[boundaryLabel.first] = 0;
		}
	}

	/**
	 * Diffusion of a field using squared distances to a neighbour points
	 * returns new point value
	 */
	field_type diffuse_one_point(uint32_t l) const
	{
		double totalSqrDist = 0.0;
		field_type result = 0.0;
		if (_node_types[l])
		{
			m_pMeshGeometry->visit_neigbour(l, [&](uint32_t l1)
			{
				double w = 1 / 
					math::sqr(m_pMeshGeometry->spacePositionOf(l1) - m_pMeshGeometry->spacePositionOf(l));
				result += _data[l1] * w;
				totalSqrDist += w;
			});

		}
		else
		{
			for (const auto& name : m_pBoundaryMesh->boundaryNames(l))
			{
				if (m_pBoundaryMesh->boundaryType(name) == BoundaryMesh::FIXED_VAL) return _data[l];
			}
			m_pMeshGeometry->visit_neigbour(l, [&](uint32_t l1)
			{
				double w = 1 /
					math::sqr(m_pMeshGeometry->spacePositionOf(l1) - m_pMeshGeometry->spacePositionOf(l));
				result += _node_types[l1] ? _data[l1] * 2.0 * w : _data[l1] * w;
				totalSqrDist += _node_types[l1] ? 2.0 * w : w;
			});
		}
		return result / totalSqrDist;
	}

	/**
	 * Returns diffused field
	 */
	field diffuse() const
	{
		field result(*this);
		for (UINT i = 0; i < _data.size(); ++i)
			result._data[i] = diffuse_one_point(i);
		return result;
	}

	/**
	 * Interpolate field value into a given point
	 * It is better when track_label is a clossest point to a {x,y,z}
	 */
	field_type interpolate(double x, double y, double z, uint32_t * track_label = nullptr) const
	{
		uint32_t start_label = track_label ? *track_label : 0;
		//track_info result = closest_plane_interpolation(x, y, z, start_label);
		mesh_geom::InterpCoefs coefs = m_pMeshGeometry->interpCoefs(x, y, z, start_label);

		if (track_label) *track_label = coefs.begin()->first;
		return std::accumulate(coefs.begin(),coefs.end(),0.0,
			[&](field_type val, mesh_geom::InterpCoef c)->field_type
		{
			return val += _data[c.first] * c.second;
		});
	}

};

#endif // !_FIELD_H