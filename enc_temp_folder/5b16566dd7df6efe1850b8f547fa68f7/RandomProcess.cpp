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

const char* RandomProcess::rndProcName(RandomProcessType nType)
{
  switch(nType)
  {
    case DIFFUSION_VELOCITY_JUMP: return _T("Velocity Jump");
    case DIFFUSION_COORD_JUMP: return _T("Coordinates Jump");
  }

  return _T("None");
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

  // [MS] 16-05-2017
  double b = .5*(i2.mob + i1.mob);  // ion mobility recalculated to current conditions, P and T.
  double D = Const_Boltzmann * T * b / m_fIonCharge;  // ion diffusion coefficient, the Einstein relation.

	//Average ion velocity
	double v0 = sqrt(6. * D / h);
  // [/MS]

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

  // [MS] 16-05-2017
  double b = .5*(i2.mob + i1.mob);  // ion mobility recalculated to current conditions, P and T.
  double D = Const_Boltzmann * T * b / m_fIonCharge;  // ion diffusion coefficient, the Einstein relation.

	//Diffusion jump calculation
	double dx = sqrt(6. * D * h);
  // [/MS]

	//Calculate vector randomly distributed on a sphere
	double phi = 2. * Const_PI * rand();
	double thetta = acos(1 - 2 * rand());
	EP::Vector3D dr = dx * EP::Vector3D{ sin(phi)*sin(thetta), cos(phi)*sin(thetta), cos(thetta) };

	ret.pos += dr;
	return ret;
}
