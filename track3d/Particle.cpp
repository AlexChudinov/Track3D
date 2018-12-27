#include "stdafx.h"
#include "Symmetry.hpp"
#include "Particle.h"
#include "EvaporationModel.h"
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
	using Reflector = math::Reflector<Vector3D, double, CTracker>;

	const CTracker * tr = CParticleTrackingApp::Get()->GetTracker();
	Reflector reflector
	(
		tr, 
		static_cast<Reflector::uint32>(tr->get_sym_plane())
	);

	Vector3D
		reflCoefs = reflector.reflectionCoefs(s.mPos),
		reflPos = reflCoefs && s.mPos,
		reflVel = reflCoefs && s.mVel;

	const CElem3D * elem = CAnsysMesh::find_global_elem
	(
		CAnsysMesh::get_global_elements()[s.mElemIdx],
		reflPos
	);

	const double h = tr->get_time_step();

	if (elem)
	{
		s.mElemIdx = elem->nInd;
		ds.mPos = s.mVel;
		CNode3D node;
		elem->interpolate(reflPos, node);
		double fExpCoeff, fMob;
		ds.mVel = reflCoefs &&
			tr->get_ion_accel(node, reflVel, t, h, s.mPhase, s.mCurr, fExpCoeff, fMob);
		double tInf;
		ds.mTemp = tr->get_dTi(node, s.mVel, s.mTemp, fExpCoeff, tInf);
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
	mMass *= h;
	mTemp *= h;
	return *this;
}

Droplet & Droplet::operator+=(const Droplet & s)
{
	mPos += s.mPos;
	mVel += s.mVel;
	mMass += mMass;
	mTemp += mTemp;
	return *this;
}

void Droplet::diff(const Droplet & s, Droplet & ds, double t)
{
	using namespace EvaporatingParticle;

	using Reflector = math::Reflector<Vector3D, double, CTracker>;

	const CTracker * tr = CParticleTrackingApp::Get()->GetTracker();
	Reflector reflector
	(
		tr,
		static_cast<Reflector::uint32>(tr->get_sym_plane())
	);

	Vector3D
		reflCoefs = reflector.reflectionCoefs(s.mPos),
		reflPos = reflCoefs && s.mPos,
		reflVel = reflCoefs && s.mVel;

	const CElem3D * elem = CAnsysMesh::find_global_elem
	(
		CAnsysMesh::get_global_elements()[s.mElemIdx],
		reflPos
	);

	const double h = tr->get_time_step();

	if (elem)
	{
		s.mElemIdx = elem->nInd;
		ds.mPos = s.mVel;
		CNode3D node;
		elem->interpolate(reflPos, node);
		double fD = tr->get_particle_diameter(s.mMass), fRe;
		ds.mVel = reflCoefs && tr->get_accel(node, reflVel, s.mMass, fD, t, fRe);
		CEvaporationModel * evapor = tr->get_evapor_model();
		ds.mMass = -evapor->get_evaporation_rate(node, s.mTemp, fD, fRe);
		ds.mTemp = -evapor->get_cooling_rate(node, s.mTemp, fD, fRe);
	}
	else
	{
		s.mElemIdx = -1;
	}
}
