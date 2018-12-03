#include "stdafx.h"
#include "Integrator.h"

IntegratorState & IntegratorState::operator*=(double /*h*/)
{
	assert(false); //Dummy function
	return *this;
}

IntegratorState & IntegratorState::operator+=(const IntegratorState &)
{
	assert(false); //Dummy function
	return *this;
}

void IntegratorState::diff(const IntegratorState & s, IntegratorState & ds, const double t)
{
	s.d(ds, t);
}

void IntegratorState::d(IntegratorState & /*dS*/, double /*t*/) const
{
	assert(false); //Dummy function
}

Integrator::Ptr Integrator::create(Type type)
{
	using namespace boost::numeric::odeint;
	switch (type)
	{
	case ExplicitEuler: return Ptr(new IntegratorImpl<euler>);
	case ModifiedMidpoint: return Ptr(new IntegratorImpl<modified_midpoint>);
	case RungeKutta4: return Ptr(new IntegratorImpl<runge_kutta4>);
	case RungeKuttaCashKarp54: return Ptr(new IntegratorImpl<runge_kutta_cash_karp54>);
	case RungeKuttaDopri5: return Ptr(new IntegratorImpl<runge_kutta_dopri5>);
	case RungeKuttaFehlberg78: return Ptr(new IntegratorImpl<runge_kutta_fehlberg78>);
	}
	return Ptr();
}
