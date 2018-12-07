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
	EvaporatingParticle::CNode3D node;
	EvaporatingParticle::CElem3D * elem = tr->get_elems()[s.mIdx];
	
	Vector3D a = tr->get_ion_accel(node, )

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
