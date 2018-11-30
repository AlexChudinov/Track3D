#include "stdafx.h"
#include "Integrator.h"

void IntegratorExplicitEuler::do_step(IntegratorState * s0, double t0, double dt)
{
	mStepper.do_step(IntegratorState::diff, s0, t0, s0, dt);
}

IntegratorState::Ptr & IntegratorState::operator*=(Ptr & s, double h)
{
	s->operator*=(h);
	return s;
}

IntegratorState::Ptr & IntegratorState::operator+=(Ptr & s0, const Ptr & s1)
{
	s0->operator+=(*s1);
	return s0;
}

void IntegratorState::diff(const IntegratorState & s, IntegratorState & ds, const double t)
{
	s.d(ds, t);
}
