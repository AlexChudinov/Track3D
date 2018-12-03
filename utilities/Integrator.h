#pragma once
#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include <boost/numeric/odeint.hpp>

class IntegratorState :
	boost::additive1< IntegratorState,
	boost::additive2< IntegratorState, double,
	boost::multiplicative< IntegratorState, double > > >
{
public:
	virtual IntegratorState& operator*=(double h);
	virtual IntegratorState& operator+=(const IntegratorState&);

	//Differentiates state
	static void diff(const IntegratorState& s, IntegratorState& ds, const double t);

private:
	//Calculates state derivatives at time t and saves them into dS
	virtual void d(IntegratorState& dS, double t) const;
};

class Integrator
{
public:
	using Ptr = std::unique_ptr<Integrator>;

	enum Type
	{
		ExplicitEuler,
		ModifiedMidpoint,
		RungeKutta4,
		RungeKuttaCashKarp54,
		RungeKuttaDopri5,
		RungeKuttaFehlberg78
	};

	static Ptr create(Type type);

	//Change s0 with value of state after step dt was calculated
	virtual void doStep(IntegratorState& s0, double t0, double dt) = 0;
};

template<template<class...> class Stepper>  
class IntegratorImpl : public Integrator
{
	Stepper<
		IntegratorState,
		double,
		IntegratorState,
		double,
		boost::numeric::odeint::vector_space_algebra> mStepper;

public:
	void doStep(IntegratorState& s0, double t0, double dt)
	{
		mStepper.do_step(IntegratorState::diff, s0, t0, dt);
	}
};

#endif