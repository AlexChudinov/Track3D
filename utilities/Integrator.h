#pragma once
#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include <boost/numeric/odeint.hpp>

template<template<class...> class Stepper, class State>
class IntegratorImpl;

template<class State>
class Integrator
{
public:
	using Ptr = std::unique_ptr<Integrator>;

	enum Type
	{
		Euler,
		ModifiedMidpoint,
		RungeKutta4,
		RungeKuttaCashKarp54,
		RungeKuttaDopri5,
		RungeKuttaFehlberg78
	};

	static Ptr create(Type type)
	{
		using namespace boost::numeric::odeint;
		switch (type)
		{
		case Euler: return Ptr(new IntegratorImpl<euler>);
		case ModifiedMidpoint: return Ptr(new IntegratorImpl<modified_midpoint>);
		case RungeKutta4: return Ptr(new IntegratorImpl<runge_kutta4>);
		case RungeKuttaCashKarp54: return Ptr(new IntegratorImpl<runge_kutta_cash_karp54>);
		case RungeKuttaDopri5: return Ptr(new IntegratorImpl<runge_kutta_dopri5>);
		case RungeKuttaFehlberg78: return Ptr(new IntegratorImpl<runge_kutta_fehlberg78>);
		}
		return Ptr();
	}

	//Change s0 with value of state after step dt was calculated
	virtual void doStep(State& s0, double t0, double dt) = 0;
};

template<template<class...> class Stepper, class State>  
class IntegratorImpl : public Integrator<State>
{
	Stepper<
		State,
		double,
		State,
		double,
		boost::numeric::odeint::vector_space_algebra> mStepper;

public:
	void doStep(State& s0, double t0, double dt)
	{
		mStepper.do_step(State::diff, s0, t0, dt);
	}
};

#endif