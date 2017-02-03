#pragma once
#ifndef _MONTECARLO_
#define _MONTECARLO_
/**
 * @brief Monte Carlo random generators
 */

#include <cstdlib>
#include <cmath>
#include "constant.hpp"
#include "vector3d.hpp"

namespace math
{
	using EvaporatingParticle::Vector3D;
	/**
	 * @brief Initializes random number generators
	 * @param seed Seed for random number generators
	 */
	inline void initGenerators(unsigned int seed)
	{
		std::srand(seed);
	}

	/**
	* Uniform [0,1)
	*/
	inline double rand()
	{
		return double(std::rand()) / double(RAND_MAX);
	}

	/**
	* Normal N [0,1) Box Muller algorithm
	*/
	inline double randn()
	{
		return std::sqrt(-2 * std::log(rand())) * std::cos(Const_2PI*rand());
	}

	/**
	 * @brief Generates molecule with normally distributed velocity components.
	 * @param uIon normalized ion velocity (in sqrt(kT/mGas))
	 * @returns Molecule velocity components
	 */
	inline Vector3D randMolecule(const Vector3D& uIon)
	{
		Vector3D testMol(randn(), randn(), randn());

		while (abs(uIon - testMol) < 6.*rand())
		{
			testMol = Vector3D(randn(), randn(), randn());
		}

		return testMol;
	}

	/**
	 * @brief Calculates vector of unit length with isotropically distributed direction
	 * @return Vector of unit length with isotropically distributed direction
	 */
	inline Vector3D randOnSphere()
	{
		double phi    = Const_2PI * rand();
		double thetta = std::acos(2.*rand() - 1.);

		return Vector3D(std::cos(thetta),
			std::sin(thetta) * std::cos(phi),
			std::sin(thetta) * std::sin(phi));
	}

	/**
	 * @brief Generates random point inside a circle
	 * @param r0 Circle radius
	 * @param pos Circle position
	 * @param norm Circle normal
	 */
	inline Vector3D randInCircle
	(
		double r0, 
		const Vector3D& pos, 
		const Vector3D& norm
	)
	{
		double phi = Const_2PI*rand();

		if (norm.length() < Const_Almost_Zero)
		{ //Generate circle with the z-normal this case
			return r0*std::sqrt(rand())*Vector3D(std::cos(phi), std::sin(phi), 0) + pos;
		}

		Vector3D en = norm / norm.length();

		Vector3D x, y; //Coordinate vectors of the circle
		do
		{
			x = Vector3D(rand(), rand(), rand());
			y = x*en;
		} while (y.length() < Const_Almost_Zero);
		x.normalize();
		y.normalize();

		return r0*std::sqrt(rand())*(x*std::cos(phi) + y*std::sin(phi)) + pos;
	}
}

#endif // !_MONTECARLO_