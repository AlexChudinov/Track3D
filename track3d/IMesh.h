#ifndef _IMESH_
#define _IMESH_

#include "vector3d.hpp"
#include <vector>
#include "Base.h"
#include "Symmetry.hpp"

/**
 * AC 02/07/2016 Mesh interface
 */

typedef EvaporatingParticle::Vector3D V3D;
typedef std::vector<V3D> VV3D;
typedef std::vector<double> Vd;

/**
 * Mesh interface class
 */
class IMesh 
{ 
public: 
	virtual ~IMesh(){} 

	//Sets the node for a given position. 
	virtual bool setNode(const V3D& pos) = 0;

	//Returns vector field values connected to a mesh node
	virtual const VV3D& getVectorVals() const = 0;

	//Returns scalar field values connected to a mesh node
	virtual const Vd& getScalarVals() const = 0;
};

/**
 * General description of a mesh class
 */
template
<
	class _M, //Class containig information about mesh values and mesh connectivity
	class _E, //Type of a mesh elementary volume element
	class _N  //Type of a mesh node element
>
class Mesh : public IMesh
{
protected:
	typedef std::vector<_E*> _VE;
	typedef bool(_E::*_IsInside)(const V3D&) const;
	typedef _VE (_E::*_GetNeighborElems)() const; //with coppying, TODO: Change CElem3D implementation to hold vector of neigbor elements
	typedef Vd(_E::*_GetScalars)(const V3D&) const;
	typedef VV3D(_E::*_GetVectors)(const V3D&) const;
	typedef std::size_t(_E::*_GetIdx)() const;
	typedef const _E* (_M::*_SearchFromBegin)(const V3D&) const;
	typedef const _E* (_M::*_SearchFromIdx)(std::size_t,const V3D&) const;
	typedef int (_M::*_Symmetry)() const;

	const _M* m_mesh;		//constant reference to a mesh
	const _E* m_currentElement;	//current mesh element
	V3D m_prevPosition; //previous ion position in the mesh
	Vd m_fieldScalarVals;   //holds scalar values of field in a current point
	VV3D m_fieldVectorVals; //holds vector values of field in a current point
	_IsInside m_inside;     //checks that a given point is inside in the m_currentElement
	_GetNeighborElems m_neighborElems; //get elements incident to a given element
	_GetScalars m_getScalars; //scalar field values at the point calculated by setNode procedure
	_GetVectors m_getVectors; //vector field values at the point calculated by setNode procedure
	_GetIdx m_getELemIdx;  //get index of a current element
	_SearchFromIdx m_searchFromIdx; //searchs the element in the mesh starting from the some index value
	_SearchFromBegin m_searchFromBegin; //searchs for the location through whole mesh
	
	V3D m_reflectionCoefs;
	math::Reflector<V3D, double, _M> m_reflector; //Providing mesh reflection transformations

	void calculateFieldValues(const V3D& reflPos)
	{
		m_fieldScalarVals = (m_currentElement->*m_getScalars)(reflPos);
		m_fieldVectorVals = (m_currentElement->*m_getVectors)(reflPos);
		for (VV3D::iterator it = m_fieldVectorVals.begin();
			it != m_fieldVectorVals.end(); ++it)
		{
			*it = (*it) && m_reflectionCoefs;
		}
	}

public:
	/**
	 * Interface constructor
	 */
	Mesh
	(
		const _M* mesh, 
		_IsInside inside, 
		_GetNeighborElems neighborElems,
		_GetScalars getScalars,
		_GetVectors getVectors,
		_GetIdx getElemIdx,
		_SearchFromIdx searchFromIdx,
		_SearchFromBegin searchFromBegin,
		_Symmetry meshSymmetry
	)
		:
		m_mesh(mesh), 
		m_currentElement(NULL),
		m_prevPosition(),
		m_inside(inside),
		m_neighborElems(neighborElems),
		m_getScalars(getScalars),
		m_getVectors(getVectors),
		m_getELemIdx(getElemIdx),
		m_searchFromIdx(searchFromIdx),
		m_searchFromBegin(searchFromBegin),
		m_reflector(m_mesh,(m_mesh->*meshSymmetry)())
	{
		setNode(m_prevPosition);
	}

	~Mesh(){}

	virtual bool setNode(const V3D& pos)
	{
		if (m_prevPosition == pos) return true;
		else m_prevPosition = pos;

		m_reflectionCoefs = m_reflector.reflectionCoefs(pos);
		V3D reflPos = m_reflectionCoefs && pos;
		if (m_currentElement)
		{
			if ((m_currentElement->*m_inside)(reflPos))
			{
				calculateFieldValues(reflPos);
				return true;
			}
			else
			{
				_VE elems = (m_currentElement->*m_neighborElems)();
				for (_VE::const_iterator itE = elems.begin();
					itE != elems.end(); ++itE)
				{
					if (((*itE)->*m_inside)(reflPos))
					{
						m_currentElement = *itE;
						calculateFieldValues(reflPos);
						return true;
					}
				}
			}
			m_currentElement = (m_mesh->*m_searchFromIdx)
			(
				(m_currentElement->*m_getELemIdx)(),
				reflPos
			);
			if (m_currentElement)
			{
				calculateFieldValues(reflPos);
				return true;
			}
			return false;
		}
		m_currentElement = (m_mesh->*m_searchFromBegin)(reflPos);
		if (m_currentElement)
		{
			calculateFieldValues(reflPos);
			return true;
		}
		return false;
	}
	virtual const VV3D& getVectorVals() const
	{ return m_fieldVectorVals; }
	virtual const Vd& getScalarVals() const
	{ return m_fieldScalarVals; }
};

/**
 * To obtain field information from the node someone needs 
 * to do some transformation with an initial field value
 */
class IField 
{
public:
	virtual ~IField(){}

	/**
	 * Sets node in the mesh
	 */
	virtual bool setNode(const V3D& pos) = 0;

	/**
	 * Returns vector fields at a location pos at a time. [MS] 20-07-2016 added the ion velocity and the current to the prototype
   * of this function to support the Coulomb field computation in the axially-symmetric case.
	 */
	virtual VV3D vectorFields(const V3D& pos, const V3D& vel, double time, double curr) const = 0;
	/**
	 * Returns scalar fields at a location pos at a time
	 */
	virtual Vd   scalarFields(const V3D& pos, double time) const = 0;

	virtual double fieldHighestFrequency() const = 0;
};

#endif