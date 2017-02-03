#pragma once
#ifndef _FIELD_OPERATOR_
#define _FIELD_OPERATOR_

#include "Field.h"

//Implementation of basic field operations in the shape of linear transforamtions

template<typename field_type>
class FieldLinearOp
{
public:
	using Field = field<field_type>;
	using mesh_geom = mesh_geometry<double, uint32_t>;
	using MatrixElem = typename mesh_geom::InterpCoef;
	using MatrixRow = typename mesh_geom::InterpCoefs;
	using Matrix = std::vector<MatrixRow>;
	using BoundaryMeshSharedPtr = std::shared_ptr<mesh_geometry<double, uint32_t>::BoundaryMesh>; 
	using MeshSharedPtr = std::shared_ptr<mesh_geometry<double, uint32_t>>;
	using InterpCoef = mesh_geometry<double, uint32_t>::InterpCoef;
	using InterpCoefs = mesh_geometry<double, uint32_t>::InterpCoefs;
	using NodeTypes = std::vector<bool>;
	using vector3f = mesh_geom::vector3f;

private:
	Matrix m_matrix;
	MeshSharedPtr m_pMeshGeometry;
	BoundaryMeshSharedPtr m_pBoundaryMesh;

	NodeTypes m_nodeTypes;

	//Adds two matrix rows
	static MatrixRow& add(MatrixRow& r1, const MatrixRow& r2)
	{

		for (const InterpCoef& c : r2)
		{
			MatrixRow::iterator it = r1.find(c.first);
			if (it == r1.end())
			{
				r1[c.first] = c.second;
			}
			else
			{
				it->second += c.second;
			}
		}
		return r1;
	}

	//Multiplicates a matrix row by a number
	static MatrixRow& mul(double h, MatrixRow& r)
	{
		for (std::pair<const uint32_t, double>& e : r) e.second *= h;
		return r;
	}

public:
	FieldLinearOp(const Field& field)
		: 
		m_matrix(field.m_pMeshGeometry->size()),
		m_pMeshGeometry(field.m_pMeshGeometry),
		m_pBoundaryMesh(field.m_pBoundaryMesh),
		m_nodeTypes(field._node_types)
	{}

	//Gets the size of a field
	size_t size() const { return m_matrix.size(); }

	//Sets inner matrix to identity
	FieldLinearOp& setToIdentity()
	{
		uint32_t i = 0;
		for (MatrixRow& row : m_matrix)
		{
			row.clear();
			row[i++] = 1.0;
		}
		return *this;
	}

	//Creates solver for equations system Ax=0, where A is laplacian
	FieldLinearOp& laplacianSolver()
	{
		for (uint32_t i = 0; i < size(); ++i)
		{
			double h = m_pMeshGeometry->shortestEdgeLength(i) / 2.0; //calculate small step
			if (m_pBoundaryMesh->isBoundary(i))
			{
				if (m_pBoundaryMesh->isFirstType(i))
				{
					m_matrix[i] = InterpCoefs{ InterpCoef{i, 1.0} };
				}
				else
				{//Zero gradient condition					
					vector3f r = m_pMeshGeometry->spacePositionOf(i) + h*m_pBoundaryMesh->normal(i);
					m_matrix[i] = m_pMeshGeometry->interpCoefs(r[0], r[1], r[2], i);
				}
			}
			else
			{
				vector3f r = m_pMeshGeometry->spacePositionOf(i);
				m_matrix[i] = m_pMeshGeometry->interpCoefs(r[0] + h, r[1], r[2], i);
				m_matrix[i] = add(m_matrix[i], m_pMeshGeometry->interpCoefs(r[0] - h, r[1], r[2], i));
				m_matrix[i] = add(m_matrix[i], m_pMeshGeometry->interpCoefs(r[0], r[1] + h, r[2], i));
				m_matrix[i] = add(m_matrix[i], m_pMeshGeometry->interpCoefs(r[0], r[1] - h, r[2], i));
				m_matrix[i] = add(m_matrix[i], m_pMeshGeometry->interpCoefs(r[0], r[1], r[2] + h, i));
				m_matrix[i] = add(m_matrix[i], m_pMeshGeometry->interpCoefs(r[0], r[1], r[2] - h, i));
				m_matrix[i] = mul(1. / 6., m_matrix[i]);
			}
		}
		return *this;
	}

	//Applies linear operator to a field
	void applyToField(Field& field) const
	{
		if (field.size() != size()) throw
			std::runtime_error("FieldLinearOp::applyToField:"
				"Field and operator sizes mismatch.");
		typename Field::data_vector data = field.data();
		size_t i = 0;
		for (const auto& row : m_matrix)
		{
			field_type& dataElem = data[i++];
			dataElem = 0.0;
			for (const auto& interpCoef : row)
			{
				dataElem += field._data[interpCoef.first] * interpCoef.second;
			}
		}
		field.data() = data;
	}
};

#endif //_FIELD_OPERATOR_