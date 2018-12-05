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
	EvaporatingParticle::CTracker * tr = CParticleTrackingApp::Get()->GetTracker();
	const EvaporatingParticle::CNodesCollection& nodes = tr->get_nodes();
	EvaporatingParticle::CNode3D node = *nodes[s.mIdx];

	double phase = 0.0;

	ds.mPos = s.mVel;

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
