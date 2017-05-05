#pragma once
#ifndef _RANDOM_PROCESS_
#define _RANDOM_PROCESS_

#include <random>
#include "TrackItem.h"

#define EP EvaporatingParticle

class RandomProcess {
public:
	using RndGen = std::mt19937_64;
	using RndGenResult = RndGen::result_type;
	using Item = EP::CIonTrackItem;
	using Distribution = std::uniform_real_distribution<double>;
	using PProcess = std::unique_ptr<RandomProcess>;

	//Base params for random process initialisation
	struct RandomProcessParams {
		RandomProcess::RndGenResult seed;
	};

	//Random process type
	enum RandomProcessType {
		DIFFUSION,
	};

	//Creates custom random process
	static PProcess createRandomProcess(const RandomProcessParams & params,
		RandomProcessType type);

	//Constructs with given seed
	RandomProcess(RndGenResult seed);

	//Returns recalculated CTrackItem taking into account start conditions {i1, n1}
	//...and finish conditions {i2, n2}
	virtual Item randomJump(const Item& i1, const Item& i2) = 0;

	//Returns uniformly distributed random value
	double rand();

	//
	virtual ~RandomProcess(){}

private:
	RndGen m_generator;
	Distribution m_distribution;
};

class Diffusion : public RandomProcess {
public:
	struct DiffusionParams : public RandomProcess::RandomProcessParams {
		double ionMass;
		double ionMobility;
		double ionCharge;
	};

	Diffusion(const DiffusionParams & params);

	virtual Item randomJump(const Item& i1, const Item& i2);

private:
	double m_fIonMass;
	double m_fIonMobility;
	double m_fIonCharge;
};

#endif