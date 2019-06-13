#pragma once
#ifndef _RANDOM_PROCESS_
#define _RANDOM_PROCESS_

#include <random>

#include "../track3d/mathematics.h"
#include "../track3d/Elements.h"
#include "TrackItem.h"

#define EP EvaporatingParticle

class RandomProcess {
public:
	using RndGen = std::mt19937_64;
	using RndGenResult = RndGen::result_type;
	using Item = EP::CIonTrackItem;
	using Node = EP::CNode3D;
	using Vector3D = EP::Vector3D;
	using ItemNode = std::pair<Item, Node>;
	using UniformDistribution = std::uniform_real_distribution<double>;
	using NormalDistribution = std::normal_distribution<double>;
	using PoissonDistribution = std::poisson_distribution<unsigned long long>;
	using PProcess = std::unique_ptr<RandomProcess>;

	//Base params for random process initialisation
	struct RandomProcessParams {
		RandomProcess::RndGenResult seed;
    bool b2D;   // [MS] 29-05-2019, 2D support.
	};

	//Random process type
	enum RandomProcessType {
		DIFFUSION_VELOCITY_JUMP,
		DIFFUSION_COORD_JUMP,
		COLLISION,
		COLLISION_ANY_PRESS
	};

	//Creates custom random process
	static PProcess createRandomProcess(const RandomProcessParams & params, RandomProcessType type);

	static const char* rndProcName(RandomProcessType nType);

	//Returns a name of a random process
	virtual const char * rndProcName() const = 0;
	//Returns a type of a random process
	virtual RandomProcessType rndProcType() const = 0;

	//Constructs with given seed
	RandomProcess(RndGenResult seed);

	//Destructor
	virtual ~RandomProcess() {}

	//Returns recalculated CTrackItem taking into account start conditions {i1, n1}
	//...and finish conditions {i2, n2}
	//Calculates random jump depended only on particle
	virtual Item randomJump(const Item& i1, const Item& i2) = 0;
	//Calculates random jump depended also on gas conditions
	virtual Item gasDependedRndJmp(const ItemNode& i1, const ItemNode& i2) = 0;

protected:
	//Returns uniformly distributed random value
	double rand();
	//Returns normally distributed random value
	double randn();
	//Returns number distributed with poisson distribution with mean lambda
	unsigned long long randPoisson(double lambda);

	//Returns vector isotropically distributed on unite sphere
	Vector3D randOnSphere();

	//Returns random molecule with respect to a given ion velocity
	Vector3D rndMol(const Vector3D& vrel, double Tgas, double Mgas);

  void set_2D(bool b2D)  { m_b2D = b2D; }   // [MS] 29-05-2019, 2D support.

private:
	RndGen m_generator;
	UniformDistribution m_uniDistrib;
	NormalDistribution m_normalDistrib;
	PoissonDistribution m_poisonDistrib;
  bool m_b2D;   // [MS] 29-05-2019, 2D support.
};

struct DiffusionParams : public RandomProcess::RandomProcessParams {
	double ionMass;
	double ionMobility;
	double ionCharge;
};

//Diffusion jump in velocity
class DiffusionVelocityJump : public RandomProcess {
public:

	DiffusionVelocityJump(const DiffusionParams & params);

	virtual const char * rndProcName() const;

	virtual RandomProcessType rndProcType() const;

	virtual Item randomJump(const Item& i1, const Item& i2);

	virtual Item gasDependedRndJmp(const ItemNode& i1, const ItemNode& i2);

private:
	double m_fIonMass;
	double m_fIonMobility;
	double m_fIonCharge;
};

//Diffusion jump in coordinate
class DiffusionCoordJump : public RandomProcess {
public:
	DiffusionCoordJump(const DiffusionParams & params);

	virtual const char * rndProcName() const;

	virtual RandomProcessType rndProcType() const;

	//Without potential correction
	virtual Item randomJump(const Item& i1, const Item& i2);

	virtual Item gasDependedRndJmp(const ItemNode& i1, const ItemNode& i2);

private:
	double m_fIonMobility;
	double m_fIonCharge;
};


//Collisional random process
class Collision : public RandomProcess {
public:
	struct CollisionParams : public RandomProcessParams {
		double ionMass;
		double ionCrossSection;
		double gasMass;
	};

	//Creates collision instance
	Collision(const CollisionParams& params);

	virtual const char * rndProcName() const;

	virtual RandomProcessType rndProcType() const;

	//Dummy function
	virtual Item randomJump(const Item& i1, const Item& i2);

	//Generates collision
	virtual Item gasDependedRndJmp(const ItemNode& in1, const ItemNode& in2);

protected:
	//Calculates average relative ion molecular speed
	double meanRelativeSpeed(const Vector3D& velIonRel, double Tgas) const;
	double m_fIonCrossSection;
	double m_fIonMass;
	double m_fGasMass;

	//Recalculate first particle velocity according to a hard sphere collision
	//with a random molecule in the gas with temperature Tgas
	Vector3D collide(double m1, const Vector3D& v1, double m2, const Vector3D& v2);
};

//Modelling collisions for every preassure
class CollisionAnyPress : public Collision
{
public:
	CollisionAnyPress(const CollisionParams& params);

	virtual const char* rndProcName() const;

	virtual RandomProcessType rndProcType() const;

	//Dummy function
	virtual Item randomJump(const Item& i1, const Item& i2);

	//Generates collision
	virtual Item gasDependedRndJmp(const ItemNode& in1, const ItemNode& in2);
};

#endif