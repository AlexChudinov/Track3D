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

	inline const Vector3D& pos() const 
	{
		return mPos;
	}

	inline const Vector3D& vel() const
	{
		return mVel;
	}

protected:
	mutable int mElemIdx;		// mesh element which this item belongs to.

	Vector3D mPos,      // position 
		mVel;			// velocity of the item.
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

	inline double phase() const
	{
		return mPhase;
	}

	inline void setPhase(double phase)
	{
		mPhase = phase;
	}

	inline double temp() const
	{
		return mTemp;
	}
private:
	double mPhase,	  //Start phase of the rf field
		mTemp,		  // ion temperature.
		mCurr;		  // current at given track.
};

class Droplet :
	public Particle,
	public BlockAllocator<Droplet>,
	boost::additive1< Droplet,
	boost::multiplicative2<Droplet, double> >
{
public:
	using Particle::operator delete;
	virtual void deleteObj();

	Droplet& operator*=(double h);

	Droplet& operator+=(const Droplet& s);

	static void diff(const Droplet& s, Droplet& ds, double t);

	inline double mass() const
	{
		return mMass;
	}

	inline double temp() const
	{
		return mTemp;
	}
private:
	double mMass, //droplet mass
		mTemp;    //droplet temperature
};

extern template class Integrator<Ion>;

extern template class Integrator<Droplet>;

#endif // !_PARTICLE_H_
