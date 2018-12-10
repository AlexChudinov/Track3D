#include "stdafx.h"
#include "Particle.h"
#include "ParticleTracking.h"

void Particle::operator delete(void * ptr, size_t n)
{
	((Particle*)ptr)->deleteObj();
}

void Ion::deleteObj()
{
	BlockPool<Ion>::freeBlock(this);
}

Ion & Ion::operator*=(double h)
{
	mPos *= h;
	mVel *= h;
	mTemp *= h;
	return *this;
}

Ion & Ion::operator+=(const Ion & s)
{
	mPos += s.mPos;
	mVel += s.mVel;
	mTemp += s.mTemp;
	return *this;
}

void Ion::diff(const Ion & s, Ion & ds, double t)
{
	using namespace EvaporatingParticle;

	const CElem3D * elem = CAnsysMesh::find_global_elem
	(
		CAnsysMesh::get_global_elements()[s.mElemIdx],
		s.mPos
	);

	const CTracker * tr = CParticleTrackingApp::Get()->GetTracker();
	const double h = tr->get_time_step();

	if (elem)
	{
		s.mElemIdx = elem->nInd;
		ds.mPos = s.mVel;
		CNode3D node;
		elem->interpolate(s.mPos, node);
		double fExpCoeff, fMob;
		ds.mVel = tr->get_ion_accel(node, s.mVel, t, h, s.mPhase, s.mCurr, fExpCoeff, fMob);
	}
	else
	{
		s.mElemIdx = -1;
	}
}

void Droplet::deleteObj()
{
	BlockPool<Droplet>::freeBlock(this);
}

Droplet & Droplet::operator*=(double h)
{
	mPos *= h;
	mVel *= h;
	return *this;
}

Droplet & Droplet::operator+=(const Droplet & s)
{
	mPos += s.mPos;
	mVel += s.mVel;
	return *this;
}

void Droplet::diff(const Droplet & s, Droplet & ds, double t)
{
}
