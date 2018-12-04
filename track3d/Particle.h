#pragma once
#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "../utilities/Integrator.h"
#include "../utilities/MemoryPool.h"
#include "vector3d.hpp"

namespace EvaporatingParticle
{
	class CTracker;
}

class Particle
{
public:
	using Vector3D = EvaporatingParticle::Vector3D;
	using CTracker = EvaporatingParticle::CTracker;

	static void operator delete(void* ptr, size_t n);
	virtual void deleteObj() = 0;

protected:
	mutable int mIdx;		// mesh element which this item belongs to.

	Vector3D mPos,      // position 
		mVel;	   // velocity of the item.
};

class Ion : 
	public Particle, 
	public BlockAllocator<Ion>,
	boost::additive1< Ion,
	boost::multiplicative2<Ion, double> >
{
public:
	using Particle::operator delete;
	virtual void deleteObj();

	Ion& operator*=(double h);

	Ion& operator+=(const Ion& s);

	static void diff(const Ion& s, Ion& ds, double t);

private:
	double mTemp,     // ion temperature.
		mTempInf,     // steady-state (equilibrium) ion temperature.
		mUnfragm,     // part of unfragmented ions.
		mMob;         // ion mobility, random diffusion velocity jumps support.
};

extern template class Integrator<Ion>;

#endif // !_PARTICLE_H_
