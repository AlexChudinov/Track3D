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
	case COLLISION_ANY_PRESS:
		return PProcess(new CollisionAnyPress((const Collision::CollisionParams&)params));
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
	case COLLISION: return _T("Random Collisions");
  case COLLISION_ANY_PRESS: return _T("Collisions for any Pressure"); // [MS] 16-10-2017 UI support.
	}

	return _T("None");
}

RandomProcess::RandomProcess(RndGenResult seed)
	:
	m_generator(seed),
	m_uniDistrib(0.0, 1.0),
	m_normalDistrib(0.0, 1.0),
	m_poisonDistrib()
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

unsigned long long RandomProcess::randPoisson(double lamdba)
{
	m_poisonDistrib.param(PoissonDistribution::param_type(lamdba));
	return m_poisonDistrib(m_generator);
}

EP::Vector3D RandomProcess::randOnSphere()
{
	double phi = 2. * Const_PI * rand();
// [MS] 29-05-2019, 2D support. 
  if(m_b2D)
    return Vector3D(sin(phi), cos(phi), 0);
// [/MS]
	double thetta = acos(1 - 2 * rand());
	return Vector3D{ sin(phi)*sin(thetta), cos(phi)*sin(thetta), cos(thetta) };
}

RandomProcess::Vector3D RandomProcess::rndMol(const Vector3D & vrel, double Tgas, double Mgas)
{
	Vector3D res;

	double sigma = sqrt(Const_Boltzmann * Tgas / Mgas);
	double vbig = sqrt((vrel & vrel) + 9.0 * sigma * sigma);

	do {
// [MS] 29-05-2019, 2D support. 
		res = m_b2D ? sigma * Vector3D{ randn(), randn(), 0 } : sigma * Vector3D{ randn(), randn(), randn() };
// [/MS]
	} while ((res - vrel).length() < vbig * rand());

	return res;
}


DiffusionVelocityJump::DiffusionVelocityJump(const DiffusionParams & params)
	:
	RandomProcess(params.seed),
	m_fIonMass(params.ionMass),
	m_fIonMobility(params.ionMobility),
	m_fIonCharge(params.ionCharge)
{
  set_2D(params.b2D); // [MS] 29-05-2019, 2D support.
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
  set_2D(params.b2D); // [MS] 29-05-2019, 2D support.
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
  set_2D(params.b2D); // [MS] 29-05-2019, 2D support.
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
	Vector3D velIonRel = (in1.first.vel + in2.first.vel) / 2. - Ugas;

	if (meanRelativeSpeed(velIonRel, Tgas) * nGas * m_fIonCrossSection * h > rand()) {
		Vector3D velMol = rndMol(velIonRel, Tgas, m_fGasMass);
		res.vel -= in2.second.vel;
		res.vel = collide(m_fIonMass, res.vel, m_fGasMass, velMol);
		res.vel += in2.second.vel;
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

Collision::Vector3D Collision::collide(double m1, const Vector3D & v1, double m2, const Vector3D& v2)
{
	//Absolute ion-molecular velocity
	double dvAbs = (v1 - v2).length();

	//Total mass of the system
	double M = m1 + m2;

	//Center mass velocity
	Vector3D vC = (m1 * v1 + m2 * v2) / M;

	return vC + m2 / M * dvAbs * randOnSphere();
}

CollisionAnyPress::CollisionAnyPress(const CollisionParams & params)
	:
	Collision(params)
{
  set_2D(params.b2D); // [MS] 29-05-2019, 2D support.
}

const char * CollisionAnyPress::rndProcName() const
{
	return _T("Collisions for any preassure.");
}

CollisionAnyPress::RandomProcessType CollisionAnyPress::rndProcType() const
{
	return COLLISION_ANY_PRESS;
}

CollisionAnyPress::Item CollisionAnyPress::randomJump(const Item & i1, const Item & i2)
{
	return i2;
}

CollisionAnyPress::Item CollisionAnyPress::gasDependedRndJmp(const ItemNode & in1, const ItemNode & in2)
{
	Item res = in2.first;
	double h = in2.first.time - in1.first.time;
	double Tgas = (in1.second.temp + in2.second.temp) / 2.;
	double Pgas = (in1.second.press + in2.second.press) / 2.;
	Vector3D Ugas = (in1.second.vel + in2.second.vel) / 2.;
	double nGas = Pgas / Tgas / Const_Boltzmann;
	Vector3D velIonRel = (in1.first.vel + in2.first.vel) / 2. - Ugas;

	double lambda = meanRelativeSpeed(velIonRel, Tgas) * nGas * m_fIonCrossSection * h;
	unsigned long long nCollisions = randPoisson(lambda);

	if (nCollisions != 0) {
		double fQuasiMolMass = nCollisions*m_fGasMass;
		Vector3D velMol = rndMol(velIonRel, Tgas, fQuasiMolMass);
		res.vel -= in2.second.vel;
		res.vel = collide(m_fIonMass, res.vel, fQuasiMolMass, velMol);
		res.vel += in2.second.vel;
	}

	return res;
}
