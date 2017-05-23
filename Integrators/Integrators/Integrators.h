#pragma once
#ifndef _Integrators_
#define _Integrators_

#include <cstddef>
#include <linearAlgebra/vectorTemplate.h>
#include <string>
#include <boost/numeric/odeint.hpp>
#include <functional>

/**
 * AC 05/10/2016
 */

namespace math
{
    using namespace boost::numeric::odeint;

    //Differentiation function type
    typedef void (__cdecl * DFUN) (void*, const double*, double*, const double*);

    enum STEPPER_TYPE_ID
    {
        BOOST_EXPLICIT_EULER          = 0,
        BOOST_MODIFIED_MIDPOINT       = 1,
        BOOST_RUNGE_KUTTA4            = 2,
        BOOST_RUNGE_KUTTA_CASH_KARP54 = 3,
        BOOST_RUNGE_KUTTA_DOPRI5      = 4,
        BOOST_RUNGE_KUTTA_FEHLBERG78  = 5,
        //BOOST_ADAMS_BASHFORTH         = 6,
        //BOOST_ADAMS_BASHFORTH_MOULTON = 7,
        //BOOST_ADAMS_MOULTON           = 8,
        RUNGE_KUTTA4                  = 9,
        RUNGE_KUTTA2                  = 10,
        UNSUPPORTED_TYPE              = 0xFF
    };

    template
    <class state_type,
    class time_type,
    class dxdt_type,
    class dt_type,
    class dummy_algebra>
    struct math_runge_kutta4
    {
        template<class diff_fun> inline
        void do_step
        (
            diff_fun diff_fun_,
            state_type& state,
            time_type time,
            dt_type dt
        )
        {
            state_type k1, k2, k3, k4;
            dt_type dh = dt/2.;
            diff_fun_(state, k1, time);
            diff_fun_(state+k1*dh, k2, time+dh);
            diff_fun_(state+k2*dh, k3, time+dh);
            diff_fun_(state+k3*dt, k4, time+dt);
            state += (k1 + 2.*(k2 + k3) + k4)*dh/3.;
        }
    };

    template
    <class state_type,
    class time_type,
    class dxdt_type,
    class dt_type,
    class dummy_algebra>
    struct math_runge_kutta2
    {
        template<class diff_fun> inline
        void do_step
        (
            diff_fun diff_fun_,
            state_type& state,
            time_type time,
            dt_type dt
        )
        {
            state_type k1, k2;
            dt_type dh = dt/2.;
            diff_fun_(state, k1, time);
            diff_fun_(state+k1*dh, k2, time+dh);
            state += (k1 + k2)*dh;
        }
    };

    struct base_integrator_interface
    {
        virtual ~base_integrator_interface(){}
        virtual void do_step
        (
            double* state, const double* time, const double* dt
        ) {}

        virtual STEPPER_TYPE_ID get_stepper_type_id() const
        {
            return UNSUPPORTED_TYPE;
        }
    };

    template
    <
    size_t state_size,
    STEPPER_TYPE_ID stepper_id,
    template<class...> class stepper,
    class diff_fun
    >
    struct integrator_interface : base_integrator_interface
    {
        using state_type = vector_c<double, state_size>;
        using stepper_type
            = stepper<state_type, double, state_type, double, vector_space_algebra>;

        stepper_type stepper_;
        diff_fun diff_fun_;

        inline integrator_interface(diff_fun _diff_fun)
            : stepper_(), diff_fun_(_diff_fun)
        {}

        void do_step(double* state, const double* time, const double* dt)
        {
            state_type& ref_state
                = *reinterpret_cast<state_type*>(state);

            stepper_.do_step(std::cref(diff_fun_), ref_state, *time, *dt);
        }

        virtual STEPPER_TYPE_ID get_stepper_type_id() const { return stepper_id; }
    };

    //Creates new integrator depending on a state size
    template <size_t N>
    base_integrator_interface* New
    (
        STEPPER_TYPE_ID stepper_id,
        void* pObj, //big object owner of a diff_fun
        //differentiation function diff_fun(pObj, state[in], dxdt[out], t)
        DFUN diff_fun
    )
    {
        using state_type = vector_c<double, N>;

		auto diff_fun_wrapper = [pObj, diff_fun](const state_type& state, state_type& dxdt, double time)->void
		{
			(*diff_fun)(pObj, state.data(), dxdt.data(), &time);
		};

        #define DEF_CASE(id, name)\
        case id: return new math::integrator_interface \
			<N, id, name, decltype(diff_fun_wrapper)> (diff_fun_wrapper);

        switch(stepper_id)
        {
        DEF_CASE(BOOST_EXPLICIT_EULER, euler);
        DEF_CASE(BOOST_RUNGE_KUTTA4, runge_kutta4);
        DEF_CASE(BOOST_RUNGE_KUTTA_CASH_KARP54, runge_kutta_cash_karp54);
        DEF_CASE(BOOST_RUNGE_KUTTA_DOPRI5, runge_kutta_dopri5);
        DEF_CASE(BOOST_RUNGE_KUTTA_FEHLBERG78, runge_kutta_fehlberg78);
        DEF_CASE(RUNGE_KUTTA4, math_runge_kutta4);
        DEF_CASE(RUNGE_KUTTA2, math_runge_kutta2);
        //DEF_CASE(BOOST_ADAMS_BASHFORTH, adams_bashforth);
        //DEF_CASE(BOOST_ADAMS_BASHFORTH_MOULTON, adams_bashforth_moulton);
        //DEF_CASE(BOOST_ADAMS_MOULTON, adams_moulton);
        DEF_CASE(BOOST_MODIFIED_MIDPOINT, modified_midpoint);
        default:
            return nullptr;
        }

        #undef DEF_CASE
    }

    template<size_t N>
    struct instantiate_states
    {
        inline math::base_integrator_interface* instantiate_next
        (
            size_t n,
            STEPPER_TYPE_ID stepper_id,
            void* pObj,
            DFUN diff_fun
        ) const
        {
            return n == N ? New<N>(stepper_id, pObj, diff_fun)
                : instantiate_states<N-1>().instantiate_next(n,stepper_id,pObj,diff_fun);
        }
    };

    template<>
    struct instantiate_states<0>
    {
        inline math::base_integrator_interface* instantiate_next
        (
            size_t n,
            STEPPER_TYPE_ID stepper_id,
            void* pObj,
            DFUN diff_fun
        ) const
        {
            return nullptr;
        }
    };
}

#endif // !_Integrators_
