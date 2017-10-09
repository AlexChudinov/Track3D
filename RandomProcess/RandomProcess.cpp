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
	switch (nType)
	{
	case DIFFUSION_VELOCITY_JUMP: return _T("Velocity Jump");
	case DIFFUSION_COORD_JUMP: return _T("Coordinates Jump");
	case COLLISION: return _T("Collision");
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
	double vbigSqr = (vrel & vrel) + 9.0 * sigma * sigma;

	do {
		res = sigma * Vector3D{ randn(), randn(), randn() };
	} while ((res - vrel).sqlength() < vbigSqr * rand());

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

DiffusionVelocityJump::RandomProcessType DiffusionVelocityJump::rndProcType() const
{
	return DIFFUSION_VELOCITY_JUMP;
}

DiffusionVelocityJump::Item DiffusionVelocityJump::randomJump(const Item & i1, const Item & i2)
{
	//Average ion temperature
	double T = .5*(i2.temp + i1.temp);

	//Small time step
	double h = i2.time - i1.time;

	// ion mobility recalculated to current conditions, P and T.
	double kappa = .5*(i2.mob + i1.mob);

	// ion diffusion coefficient, the Einstein relation.
	double D = Const_Boltzmann * T * kappa / m_fIonCharge;

	//Average ion velocity													
	double v0 = sqrt(6. * D / h);

	Vector3D v = v0 * randOnSphere();

	Item ret = i2;
	ret.vel += v;
	return ret;
}

DiffusionVelocityJump::Item DiffusionVelocityJump::gasDependedRndJmp(const ItemNode & in1, const ItemNode & in2)
{
	Item ret = in2.first;

	//Average ion temperature
	double T = .5*(in2.first.temp + in1.first.temp);

	//Average gas pressure
	double P = .5*(in2.second.press + in1.second.press);

	//Small time step
	double h = in2.first.time - in1.first.time;

	//Ion mobility
	double kappa = m_fIonMobility * sqrt(T / Const_T0) * Const_One_Atm_CGS / P;

	//Characteristic relaxation time
	double tau = kappa * m_fIonMass / m_fIonCharge;

	// Diffusion coefficient
	double D = Const_Boltzmann * T * kappa / m_fIonCharge;

	//Abdolute jump value
	double dvAbs = sqrt(6. * D * h) / (tau * (1 - exp(-h / tau)));

	Vector3D dv = dvAbs * randOnSphere();

	ret.vel += dv;

	return ret;
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

DiffusionCoordJump::RandomProcessType DiffusionCoordJump::rndProcType() const
{
	return DIFFUSION_COORD_JUMP;
}

DiffusionCoordJump::Item DiffusionCoordJump::randomJump(const Item & i1, const Item & i2)
{
	Item ret = i2;

	//Average ion temperature
	double T = .5*(i2.temp + i1.temp);

	//Small time step
	double h = i2.time - i1.time;

	// ion mobility recalculated to current conditions, P and T.
	double kappa = .5*(i2.mob + i1.mob);  

	// ion diffusion coefficient, the Einstein relation.
	double D = Const_Boltzmann * T * kappa / m_fIonCharge; 

	//Diffusion jump calculation													
	double dx = sqrt(6. * D * h);

	Vector3D dr = dx * randOnSphere();

	ret.pos += dr;
	return ret;
}

DiffusionCoordJump::Item DiffusionCoordJump::gasDependedRndJmp(const ItemNode & in1, const ItemNode & in2)
{
	Item ret = in2.first;

	//Average ion temperature
	double T = .5*(in2.first.temp + in1.first.temp);

	//Average gas pressure
	double P = .5*(in2.second.press + in1.second.press);

	//Small time step
	double h = in2.first.time - in1.first.time;

	//Ion mobility
	double kappa = m_fIonMobility * sqrt(T / Const_T0) * Const_One_Atm_CGS / P;

	// Diffusion coefficient
	double D = Const_Boltzmann * T * kappa / m_fIonCharge; 

	//Abdolute jump value
	double dx = sqrt(6. * D * h);

	Vector3D dr = dx * randOnSphere();

	ret.pos += dr;

	return ret;
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

Collision::RandomProcessType Collision::rndProcType() const
{
	return COLLISION;
}

Collision::Item Collision::randomJump(const Item & i1, const Item & i2)
{
	return i2;
}

Collision::Item Collision::gasDependedRndJmp(const ItemNode & in1, const ItemNode & in2)
{
	Item res = in2.first;
	double h = in2.first.time - in1.first.time;
	double Tgas = (in1.second.temp + in2.second.temp) / 2.;
	double Pgas = (in1.second.press + in2.second.press) / 2.;
	Vector3D Ugas = (in1.second.vel + in2.second.vel) / 2.;
	double nGas = Pgas / Tgas / Const_Boltzmann;
	Vector3D velIonRel = ((in1.first.vel - in1.second.vel) + (in2.first.vel - in2.second.vel)) / 2.;

	if (meanRelativeSpeed(velIonRel, Tgas) * nGas * m_fIonCrossSection * h > rand()) {
		velIonRel = in2.first.vel - Ugas;
		velIonRel = collide(velIonRel, Tgas);
		res.vel = velIonRel + Ugas;
	}

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

Collision::Vector3D Collision::collide(const Vector3D & v, double Tgas)
{
	//Molecular velocity
	Vector3D vmol = rndMol(v, Tgas, m_fGasMass);

	//Absolute ion-molecular velocity
	double dvAbs = (v - vmol).length();

	//Total mass of the system
	double M = m_fIonMass + m_fGasMass;

	//Center mass velocity
	Vector3D vC = (m_fIonMass * v + m_fGasMass * vmol) / M;

	return vC + m_fGasMass / M * dvAbs * randOnSphere();
}
