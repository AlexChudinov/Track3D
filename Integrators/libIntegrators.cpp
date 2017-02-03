#include "libIntegrators.h"
#include "Integrators/Integrators.h"
#include <string>

void* create_integrator_interface
(
    unsigned long integrator_state_size,
    STEPPER_TYPE_ID stepper_id,
    const void* pObj, //big object owner of diff_fun and obs_fun
    //differentiation function diff_fun(pObj, state[in], dxdt[out], t)
    void (__cdecl *diff_fun)(const void*, const double*, double*, const double*)
)
{
    return integrator_state_size > MAX_STATE_SIZE ?
                nullptr
                :
                math::instantiate_states<MAX_STATE_SIZE>()
                    .instantiate_next
					(
						integrator_state_size, 
						math::STEPPER_TYPE_ID(stepper_id), 
						pObj, 
						diff_fun);
}

void delete_integrator_interface(void* ptr_integrator_interface)
{
    delete reinterpret_cast<math::base_integrator_interface*>(ptr_integrator_interface);
}

void do_integrator_step
(
    void* ptr_integrator_interface,
    double* state,
    const double* time,
    const double* dt
)
{
    reinterpret_cast<math::base_integrator_interface*>(ptr_integrator_interface)
        ->do_step(state, time, dt);
}
