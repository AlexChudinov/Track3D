#pragma once
#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

class IntegratorState
{
public:
	using Ptr = std::unique_ptr<IntegratorState>;

	virtual IntegratorState& operator*=(double h) = 0;
	virtual IntegratorState& operator+=(const IntegratorState&) = 0;

	static Ptr& operator*=(Ptr& s, double h);
	static Ptr& operator+=(Ptr& s0, const Ptr& s1);
	//Differentiates state
	static void diff(const IntegratorState& s, IntegratorState& ds, const double t);

private:
	//Calculates state derivatives at time t and saves them into dS
	virtual void d(IntegratorState& dS, double t) const = 0;
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
		RungeKuttaDopru5,
		RungeKuttaFehlberg78
	};

	static Ptr create(Type type);

	//Change s0 with value of state after step dt was calculated
	virtual void do_step(IntegratorState& s0, double t0, double dt) = 0;
};

class IntegratorExplicitEuler
{
	euler<IntegratorState*, double, IntegratorState*, double, vector_space_algebra>
		mStepper;
public:
	virtual void do_step(IntegratorState* s0, double t0, double dt);
};

#endif