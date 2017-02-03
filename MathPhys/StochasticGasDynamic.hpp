#pragma once
#ifndef _STOCHASTICGASDYNAMIC_
#define _STOCHASTICGASDYNAMIC_

/**
* @brief Defines stochastical gas molecular interaction
*/

#include "constant.hpp"
#include "MonteCarlo.hpp"
#include "mathematics.h"

namespace math
{
	/**
	 * @brief Calculates collision between molecule and ion
	 * @param uIon1 Ion velocity measured in sqrt(kT/mGas) after collision
	 * @param uIon0 Ion velocity measured in sqrt(kT/mGas) before collision
	 * @param mIon  Ion mass
	 * @param mMol  Gas molecule mass
	 */
	template<class Vector3D>
	inline void collision
	(
		Vector3D& uIon1,
		const Vector3D& uIon0,
		const double& mIon,
		const double& mMol
	)
	{
		Vector3D uMol = randMolecule(uIon0);
		uIon1 = (mIon*uIon0 + mMol*uMol) / (mMol + mIon)
			+ mMol / (mMol + mIon) * abs(uIon0 - uMol) * randOnSphere();
	}

	
	/**
	 * @brief Mean relative ion-gas speed averaged over a Maxwell distribution of a gas [see Simion documentation]
	 * @param uIon Ion's absolute speed in sqrt(kT/m)
	 */
	inline double meanRelativeGasSpeed(double uIon)
	{
		if(uIon == 0.0)
      return sqrt(8 / Const_PI);

		return (uIon + 1 / uIon) * erf(uIon / sqrt(2.)) +
			sqrt(2 / Const_PI) * exp(-uIon * uIon / 2);
	}

	/*
	[AC] 20/07/2016 Two ways of ion mobility and diffusion coefficient estimation:
	1) Takes into account current ion temperature dependent on ion's drift velocity
	2) Takes into account just ambient gas conditions
	*/
	/**
	Estimates ion mobility taking into account ion temperature
	*/
	template<class T>
	class MobilityCalculatorUsingIonTemperature
	{
		const T m_k; //Ion mobility

	public:
		inline MobilityCalculatorUsingIonTemperature
		(
			const T& k0, //Ion mobility at STP
			const T& Tg, //gas temperature
			const T& P,  //gas pressure
			const T& Ti  //ion temperature
		)
			:
			m_k(k0 * Tg / std::sqrt(Const_T0*Ti) * Const_One_Atm_CGS / P)
		{}

		inline const T& ionMobility () const
		{
			return m_k;
		}

		inline T ionDiffusion
		(
			const T& Ti, //Ion temperature
			const T& Tg, //Gas temperature
			const T& z //Ion charge
		) const
		{
			return Const_Boltzmann*ionMobility()*Ti / z;
		}
	};
	/**
	Estimates ion mobility using just ambient gas temperature
	*/
	template <class T>
	struct MobilityCalculatorUsingGasTemperature
	{
		const T m_k; //Ion mobility

	public:
		inline MobilityCalculatorUsingGasTemperature
		(
			const T& k0, //Ion mobility at STP
			const T& Tg, //gas temperature
			const T& P,  //gas pressure
			const T& Ti  //ion temperature
		)
			:
			m_k(k0 * std::sqrt(Tg / Const_T0) * Const_One_Atm_CGS / P)
		{}

		inline const T& ionMobility() const
		{
			return m_k;
		}

		inline T ionDiffusion
		(
			const T& Ti, //Ion temperature
			const T& Tg, //Gas temperature
			const T& z //Ion charge
		) const
		{
			return Const_Boltzmann*ionMobility()*Tg / z;
		}
	};
}

#endif

