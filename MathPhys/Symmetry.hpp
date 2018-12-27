#pragma once
#ifndef _SYMMETRY_
#define _SYMMETRY_

/**
 * Field and acceleration symmetry transformation
 * AC 24/06/07
 */

#include "Exception.h"

namespace math
{
	/**
	 * Reflection on symmetry planes. It is supposed that field calculated only for positive
	 * domain
	 */
	template<class Vector3D, class NumType, class Model>
	class Reflector
	{
	public:
		typedef unsigned int uint32;
		typedef Vector3D(Model::*VectorField) (const Vector3D&, const Vector3D&, const NumType&);
		typedef NumType(Model::*ScalarField) (const Vector3D&, const Vector3D&, const NumType&);
	
	private:
		/**
		 * Field dependence on coordinate, velocity and time: U(r,v,t)
		 */
		const Model* m_model;
		VectorField m_fieldVector;
		ScalarField m_fieldScalar;

		uint32 m_symType;

		/**
		* An axis along which reflection is supposed. There are three axes X,Y,Z,
		* which is correspond to a reflection on YZ, ZX and XY planes
		*/
		enum ReflAxis
		{
      _None = 0,  // [MS] 18-07-2016, No symmetry at all.
			_Z    = 1, /// Plane orthogonal to an Z axis
			_Y    = 2, /// Plane orthogonal to an Y axis
			_X    = 4  /// Plane orthogonal to an X axis
		};

	public:

		/**
		 * Constructs reflector object, using field functor and symmetry type
		 */
		Reflector
		(
			const Model* model,
			uint32 symType,
			VectorField fieldVector = NULL,
			ScalarField fieldScalar = NULL
		)
			:
			m_model(model), 
			m_fieldVector(fieldVector),
			m_fieldScalar(fieldScalar),
			m_symType(symType) 
		{}

		/**
		 * Get/set operations
		 */
		inline ScalarField& scalarField()
		{ return m_fieldScalar; }

		inline const ScalarField& scalarField() const
		{ return m_fieldScalar; }

		inline VectorField& vectorField()
		{ return m_fieldVector; }

		inline const VectorField& vectorField() const
		{ return m_fieldVector; }

		/**
		 * Returns true if some of the types of symmetry is occured
		 */
		inline bool X() const { return _X & m_symType; }
		inline bool Y() const { return _Y & m_symType; }
		inline bool Z() const { return _Z & m_symType; }

		/**
		 * Calculates reflection coeffs
		 * @param r Radius vector
		 * @return Vector of reflection coefs: -1 - reflection occurs 1 - reflection does not occur
		 */
		inline Vector3D reflectionCoefs(const Vector3D& r) const
		{
			Vector3D rCoefs;
			rCoefs.x = X() && r.x < 0.0 ? -1.0 : 1.0;
			rCoefs.y = Y() && r.y < 0.0 ? -1.0 : 1.0;
			rCoefs.z = Z() && r.z < 0.0 ? -1.0 : 1.0;
			return rCoefs;
		}

		/**
		 * Calculates vector field in the given symmetry type
		 * Note, operator&& stands for a element wise product of two vectors
		 */
		inline Vector3D vectorField
		(
			const Vector3D& r, 
			const Vector3D& v, 
			const NumType& time
		) const
		{
			if (!m_fieldVector) 
				throw Exception("From reflector object: function for a vector field was not set.");
			Vector3D rCoefs = reflectionCoefs(r);
			return rCoefs && (m_model->*m_fieldVector)(rCoefs && r, rCoefs && v, time);
		}

		/**
		 * Same as above for scalar field
		 */
		inline NumType scalarField
		(
			const Vector3D& r,
			const Vector3D& v,
			const NumType& time
		) const
		{
			if(!m_fieldScalar)
				throw Exception("From reflector object: function for a vector field was not set.");
			Vector3D rCoefs = reflectionCoefs(r);
			return (m_model->*m_fieldScalar)(rCoefs && r, rCoefs && v, time);
		}
	};
}

#endif // !_SYMMETRY_