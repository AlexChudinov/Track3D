#include "stdafx.h"
#include "RandomProcess.h"
#include "constant.hpp"

RandomProcess::PProcess RandomProcess::createRandomProcess(const RandomProcessParams & params,
	RandomProcessType type)
{
	switch (type) {
	case DIFFUSION_VELOCITY_JUMP:
		return PProcess(new DiffusionVelocityJump((const DiffusionParams &)params));
	case DIFFUSION_COORD_JUMP:
		return PProcess(new DiffusionCoordJump((const DiffusionParams&)params));
	case COLLISION:
		return PProcess(new Collision((const Collision::CollisionParams&)params));
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
	m_uniDistrib(0.0, 1.0),
	m_normalDistrib(0.0, 1.0)
{
}

double RandomProcess::rand()
{
	return m_uniDistrib(m_generator);
}

double RandomProcess::randn()
{
	return m_normalDistrib(m_generator);
}

EP::Vector3D RandomProcess::randOnSphere()
{
	double phi = 2. * Const_PI * rand();
	double thetta = acos(1 - 2 * rand());
	return Vector3D{ sin(phi)*sin(thetta), cos(phi)*sin(thetta), cos(thetta) };
}

RandomProcess::Vector3D RandomProcess::rndMol(const Vector3D & vrel, double Tgas, double Mgas)
{
	Vector3D res;

	double sigma = sqrt(Const_Boltzmann * Tgas / Mgas);

	do {
		res = sigma * Vector3D{ randn(), randn(), randn() };
	} while (res.length() < 3.0 * sigma * rand());

	return res;
}


DiffusionVelocityJump::DiffusionVelocityJump(const DiffusionParams & params)
	:
	RandomProcess(params.seed),
	m_fIonMass(params.ionMass),
	m_fIonMobility(params.ionMobility),
	m_fIonCharge(params.ionCharge)
{
}

const char * DiffusionVelocityJump::rndProcName() const
{
	return _T("Velocity Jump");
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

	Vector3D v = v0 * randOnSphere();

	Item ret = i2;
	ret.vel += v;
	return ret;
}

DiffusionVelocityJump::Item DiffusionVelocityJump::gasDependedRndJmp(const ItemNode & i1, const ItemNode & i2)
{
	return randomJump(i1.first, i2.first);
}

DiffusionCoordJump::DiffusionCoordJump(const DiffusionParams & params)
	:
	RandomProcess(params.seed),
	m_fIonMobility(params.ionMobility),
	m_fIonCharge(params.ionCharge)
{
}

const char * DiffusionCoordJump::rndProcName() const
{
	return _T("Coordinates Jump");
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

	Vector3D dr = dx * randOnSphere();

	ret.pos += dr;
	return ret;
}

DiffusionCoordJump::Item DiffusionCoordJump::gasDependedRndJmp(const ItemNode & i1, const ItemNode & i2)
{
	return randomJump(i1.first, i2.first);
}

Collision::Collision(const CollisionParams & params)
	:
	RandomProcess(params.seed),
	m_fIonCrossSection(params.ionCrossSection),
	m_fIonMass(params.ionMass),
	m_fGasMass(params.gasMass)
{
}

const char * Collision::rndProcName() const
{
	return _T("Collisions");
}

Collision::Item Collision::gasDependedRndJmp(const ItemNode & in1, const ItemNode & in2)
{
	Item res = in2.first;
	double h = in2.first.time - in1.first.time;
	double Tgas = (in1.second.temp + in2.second.temp) / 2.;
	Vector3D velIonRel = ((in1.first.vel - in1.second.vel) + (in2.first.vel - in2.second.vel)) / 2.;
	
	//if(meanRelativeSpeed(velIonRel, Tgas) * m_fIonCrossSection)

	return res;
}

double Collision::meanRelativeSpeed(const Vector3D & velIonRel, double Tgas) const
{
	double
		A = m_fGasMass / 2. / Const_Boltzmann / Tgas,
		meanGasSpeed = 2 / sqrt(A * Const_PI),
		s = velIonRel.length() * sqrt(A);
	return meanGasSpeed*((s + 1. / (2.*s))*sqrt(Const_PI) / 2.*EP::CMath::erf(s) + 1. / 2. * exp(-s*s));
}
