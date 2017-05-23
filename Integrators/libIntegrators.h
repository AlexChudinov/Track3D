#pragma once
#ifndef _Integrators_c_interface_
#define _Integrators_c_interface_

#define MAX_STATE_SIZE 55 //It is hard embedded into code

#define INTEGRATORS_EXPORT extern "C" _declspec(dllexport)

enum STEPPER_TYPE_ID
{
	BOOST_EXPLICIT_EULER = 0,
	BOOST_MODIFIED_MIDPOINT = 1,
	BOOST_RUNGE_KUTTA4 = 2,
	BOOST_RUNGE_KUTTA_CASH_KARP54 = 3,
	BOOST_RUNGE_KUTTA_DOPRI5 = 4,
	BOOST_RUNGE_KUTTA_FEHLBERG78 = 5,
	//BOOST_ADAMS_BASHFORTH         = 6,
	//BOOST_ADAMS_BASHFORTH_MOULTON = 7,
	//BOOST_ADAMS_MOULTON           = 8,
	RUNGE_KUTTA4 = 9,
	RUNGE_KUTTA2 = 10,
	UNSUPPORTED_TYPE = 0xFF
};

INTEGRATORS_EXPORT void* create_integrator_interface
(
	unsigned long state_size,
	STEPPER_TYPE_ID stepper_id,
	void* pObj,
	void(__cdecl *diff_fun)(void*, const double*, double*, const double*)
);

INTEGRATORS_EXPORT void delete_integrator_interface(void* ptr_integrator_interface);

INTEGRATORS_EXPORT void do_integrator_step
(
	void* ptr_integrator_interface,
	double* state,
	const double* time,
	const double* dt
);


#endif // _Integrators_c_interface_
