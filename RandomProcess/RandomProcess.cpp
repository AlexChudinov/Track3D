#include "stdafx.h"
#include "RandomProcess.h"
#include "constant.hpp"

RandomProcess::PProcess RandomProcess::createRandomProcess(const RandomProcessParams & params,
	RandomProcessType type)
{
	switch (type) {
	case DIFFUSION_VELOCITY_JUMP:
		return PProcess(new DiffusionVelocityJump((const DiffusionParams &)params));
	default:
		return PProcess();
	}
}

RandomProcess::RandomProcess(RndGenResult seed)
	:
	m_generator(seed),
	m_distribution(0.0, 1.0)
{
}

double RandomProcess::rand()
{
	return m_distribution(m_generator);
}

DiffusionVelocityJump::DiffusionVelocityJump(const DiffusionParams & params)
	:
	RandomProcess(params.seed),
	m_fIonMass(params.ionMass),
	m_fIonMobility(params.ionMobility),
	m_fIonCharge(params.ionCharge)
{
}

DiffusionVelocityJump::Item DiffusionVelocityJump::randomJump(const Item & i1, const Item & i2)
{
	//Average ion temperature
	double T = .5*(i2.temp + i1.temp);

	//Small time step
	double h = i2.time - i1.time;

	//Average ion velocity
	double v0 = sqrt(2. * m_fIonCharge * Const_Boltzmann * T * h / m_fIonMobility) / m_fIonMass;

	//Calculate vector randomly distributed on a sphere
	double phi = 2. * Const_PI * rand();
	double thetta = acos(1 - 2 * rand());
	EP::Vector3D v = v0 * EP::Vector3D{ sin(phi)*sin(thetta), cos(phi)*sin(thetta), cos(thetta) };

	Item ret = i2;
	ret.vel += v;
	return ret;
}

DiffusionCoordJump::DiffusionCoordJump(const DiffusionParams & params)
	:
	RandomProcess(params.seed),
	m_fIonMobility(params.ionMobility),
	m_fIonCharge(params.ionCharge)
{
}

DiffusionCoordJump::Item DiffusionCoordJump::randomJump(const Item & i1, const Item & i2)
{
	Item ret = i2;

	//Average ion temperature
	double T = .5*(i2.temp + i1.temp);

	//Small time step
	double h = i2.time - i1.time;

	//Diffusion jump calculation
	double dx = sqrt(2 * m_fIonMobility * Const_Boltzmann * T * h / m_fIonCharge);

	//Calculate vector randomly distributed on a sphere
	double phi = 2. * Const_PI * rand();
	double thetta = acos(1 - 2 * rand());
	EP::Vector3D dr = dx * EP::Vector3D{ sin(phi)*sin(thetta), cos(phi)*sin(thetta), cos(thetta) };

	ret.pos += dr;
	return ret;
}
