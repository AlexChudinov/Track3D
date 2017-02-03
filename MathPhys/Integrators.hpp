#pragma once
#ifndef _Integrators_
#define _Integrators_

#include <vector>
#include "Exception.h"

/**
 * AC 24/06/2016
 */

namespace math
{
	/**
	 * Runge-Kutta four order integrator template with a constant step h
	 * for system like this: y' = f(y,t).
	 */
	template<class _tY, class _tdY, class _tT, class _ttY>
	inline void rungeKuttaIntegrator4
		(
			_tT& t1,
			_tY& y1,
			_tdY(_tY::* f) (const _tT&) const,
			const _tT& t0,
			const _tY& y0,
			const _tT& h,
			const _tT& h2, //half step
			const _tT& h6 //one sixth of a step
		) throw()
	{
		_tdY k1, k2, k3, k4;
		_ttY y00 = y0;
		y1 = y00;
		k1 = (y1.*f)(t0);
		k2 = ((y1 = y00 + k1*h2).*f)(t0 + h2);
		k3 = ((y1 = y00 + k2*h2).*f)(t0 + h2);
		k4 = ((y1 = y00 + k3*h).*f)(t0 + h);
		t1 = t0 + h;
		y1 = y00 + h6*(k1 + 2. * k2 + 2. * k3 + k4);
	}

	/**
	* Second order Runge-Kutta
	*/
	template<class _tY, class _tdY, class _tT, class _ttY>
	inline void rungeKuttaIntegrator2
		(
			_tT& t1,
			_tY& y1,
			_tdY(_tY::* f) (const _tT&) const,
			const _tT& t0,
			const _tY& y0,
			const _tT& h,
			const _tT& h2 //half step
		) throw()
	{
		_ttY y00 = y0;
		_tdY k1, k2;
		y1 = y00;
		k1 = (y1.*f)(t0);
		k2 = ((y1 = y00 + k1*h2).*f)(t0 + h2);
		t1 = t0 + h;
		y1 = y00 + h2*(k1 + k2);
	}

	/**
	 * Predictor-corrector integrator
	 */
	template<class _tY, class _tdY, class _tT, class _ttY>
	inline void predictorCorrectorIntegrator
	(
		_tT& t1,
		_tY& y1,
		_tdY(_tY::* f) (const _tT&) const,
		const _tT& t0,
		const _tY& y0,
		const _tT& h,
		const _tT& h2 //half step
	) throw()
	{
		_ttY y00 = y0;
		_tdY k1, k2;
		y1 = y00;
		k1 = (y1.*f)(t0);
		k2 = ((y1 = y00 + k1*h2).*f)(t0 + h2);
		t1 = t0 + h;
		y1 = y00 + h*k2;
	}

	/**
	 * Integrators interface: parameters for integration and the integrating fun
	 */
	class IParams { public: virtual ~IParams() {} };
	class IIntegrator 
	{ 
	public: 
		virtual void doIntegratorStep(IParams*) const throw() = 0;
		virtual ~IIntegrator(){}
	};


	template<class _tY, class _tdY, class _tT>
	class Params : public IParams
	{
		typedef _tdY(_tY::*DF)(const _tT&) const;
	public:
		_tY& m_y;
		DF m_f;
		_tT& m_t;
		const _tT m_h, m_h2, m_h6;

		Params(_tY& y, DF f, _tT& t, const _tT& h)
			: m_y(y), m_f(f), m_t(t), m_h(h), m_h2(h / 2.0), m_h6(h / 6.0) 
		{}
	};

	template<class _tY, class _tdY, class _tT, class _ttY = _tdY>
	class IntegratorRK4 : public IIntegrator
	{
		typedef Params<_tY, _tdY, _tT> Pars;
	public:
		virtual void doIntegratorStep(IParams* params) const throw()
		{
			rungeKuttaIntegrator4<_tY,_tdY,_tT,_ttY>
			(
				((Pars*)params)->m_t,
				((Pars*)params)->m_y,
				((Pars*)params)->m_f,
				((Pars*)params)->m_t,
				((Pars*)params)->m_y,
				((Pars*)params)->m_h,
				((Pars*)params)->m_h2,
				((Pars*)params)->m_h6
			);
		}
	};

	template<class _tY, class _tdY, class _tT, class _ttY = _tdY>
	class IntegratorRK2 : public IIntegrator
	{
		typedef Params<_tY, _tdY, _tT> Pars;
	public:
		virtual void doIntegratorStep(IParams* params) const throw()
		{
			rungeKuttaIntegrator2<_tY, _tdY, _tT, _ttY>
				(
					((Pars*)params)->m_t,
					((Pars*)params)->m_y,
					((Pars*)params)->m_f,
					((Pars*)params)->m_t,
					((Pars*)params)->m_y,
					((Pars*)params)->m_h,
					((Pars*)params)->m_h2
				);
		}
	};

	template<class _tY, class _tdY, class _tT, class _ttY = _tdY>
	class IntegratorPC : public IIntegrator
	{
		typedef Params<_tY, _tdY, _tT> Pars;
	public:
		virtual void doIntegratorStep(IParams* params) const throw()
		{
			predictorCorrectorIntegrator<_tY, _tdY, _tT, _ttY>
				(
					((Pars*)params)->m_t,
					((Pars*)params)->m_y,
					((Pars*)params)->m_f,
					((Pars*)params)->m_t,
					((Pars*)params)->m_y,
					((Pars*)params)->m_h,
					((Pars*)params)->m_h2
				);
		}
	};

	//Integration scheme switching depending on a some condition "cond"
	class ISwitchDiffScheme
	{
	public:
		virtual ~ISwitchDiffScheme(){}
		virtual void doIntegratorStep() throw() = 0;
	};

	template<class _tY, class _tdY, class _tT>
	class SwitchDiffScheme : public ISwitchDiffScheme, public CBase
	{
		typedef std::size_t(_tY::*Cond)() const; //switching condition
		typedef std::vector<IParams*> _Pars;
		typedef _tdY(_tY::*_dFun) (const _tT&) const throw();
		typedef std::vector<_dFun> _DFuns;

		_tY& m_y;
		Cond m_cond;

		//We use a same integrator, but we are switching the params
		_Pars m_params;
		IIntegrator* m_integrator;
	public:
		SwitchDiffScheme
		(
			_tY& y, 
			const _DFuns& fs,
			_tT& t,
			const _tT& h, //Different step size for every integration scheme
			IIntegrator* integrator,
			Cond cond = nullptr
		)
			:
			m_y(y),
			m_cond(cond),
			m_params(fs.size()),
			m_integrator(integrator)
		{
			if(!integrator)
				throw Exception(toString() + ": from constructor:"
					+ " integrator was not set.");

			if (fs.empty())
				throw Exception(toString() + ": from constructor:"
					+ " functions array is empty.");

			for (std::size_t i = 0; i < fs.size(); ++i)
			{
				if (!fs[i])
					throw Exception(toString() + ": from constructor:"
						+ " null pointer in the functions array.");

				m_params[i] = new Params<_tY, _tdY, _tT>(y, fs[i], t, h);
			}
		}

		~SwitchDiffScheme()
		{
			for (_Pars::iterator it = m_params.begin(); it != m_params.end(); ++it)
				delete *it;
		}

		void doIntegratorStep() throw()
		{
			if (!m_cond)
				m_integrator->doIntegratorStep(m_params[0]);
			else
				m_integrator->doIntegratorStep(m_params.at((m_y.*m_cond)()));
		}

		std::string toString() const
		{
			return std::string("Class SwitchDiffScheme");
		}
	};

	//Params that are changed during one differenciation step
	template<class T, int N> class DiffParamOneType;
	template<class T> class DiffParamOneType<T, 1>;

	template<class T, int N>
	class DiffParamOneType : public DiffParamOneType<T,N-1>
	{
	protected:
		T t;
	public:

		template<int M> inline const T& param() const
		{ return DiffParamOneType<T,M>::t; }

		template<int M> inline T& param()
		{ return DiffParamOneType<T,M>::t; }

		DiffParamOneType<T, N> operator+(const DiffParamOneType<T, N>& p) const
		{
			DiffParamOneType<T, N> res;
			DiffParamOneType<T, N - 1> * pRes;
			const DiffParamOneType<T, N - 1> * pP;

			res.param<N>() = this->param<N>() + p.param<N>();

			pRes  =  reinterpret_cast<DiffParamOneType<T, N - 1>*>(&res);
			pP    =  reinterpret_cast<const DiffParamOneType<T, N - 1>*>(&p);
			*pRes = *reinterpret_cast<const DiffParamOneType<T, N - 1>*>(this) + *pP;

			return res;
		}

		template<class _tT>
		DiffParamOneType<T, N> operator*(const _tT& h) const
		{
			DiffParamOneType<T, N> res;
			DiffParamOneType<T, N - 1> *pRes;

			res.param<N>() = this->param<N>() * h;

			pRes = reinterpret_cast<DiffParamOneType<T, N - 1>*>(&res);
			*pRes = 
				*reinterpret_cast<const DiffParamOneType<T, N - 1>*>(this) * h;

			return res;
		}
	};

	template<class T>
	class DiffParamOneType<T, 1>
	{
	protected:
		T t;
	public:
		template<int N> T& param();
		template<int N> const T& param() const;

		template<>
		inline T& param<1>() 
		{ return t; }

		template<>
		inline const T& param<1>() const 
		{ return t; }

		DiffParamOneType<T, 1> operator+(const DiffParamOneType<T, 1>& p) const
		{
			DiffParamOneType<T, 1> res;
			res.param<1>() = this->param<1>() + p.param<1>();
			return res;
		}

		template<class _tT>
		DiffParamOneType<T, 1> operator*(const _tT& h) const
		{
			DiffParamOneType<T, 1> res;
			res.param<1>() = this->param<1>() * h;
			return res;
		}
	};

	template <class T, int N, class _tT>
	DiffParamOneType<T, N> operator*
		(
			const _tT& h,
			const DiffParamOneType<T, N>& p
		)
	{
		return p*h;
	}

	//Holds two types vector 3D and numeric
	template<class _V, int NV, class _N, int NN>
	class DiffStepParams
	{
		DiffParamOneType<_V, NV> m_vectors;
		DiffParamOneType<_N, NN> m_scalars;
	public:

		inline DiffParamOneType<_V, NV>& vectors() { return m_vectors; }
		inline const DiffParamOneType<_V, NV>& vectors() const { return m_vectors; }

		inline DiffParamOneType<_N, NN>& scalars() { return m_scalars; }
		inline const DiffParamOneType<_N, NN>& scalars() const { return m_scalars; }


		template<int N> inline
		_V& vec3D() { return m_vectors.param<N>(); }

		template<int N> inline
		const _V& vec3D() const { return m_vectors.param<N>(); }
		
		template<int N>
		_N& scalar() { return m_scalars.param<N>(); }

		template<int N> inline
		const _N& scalar() const { return m_scalars.param<N>(); }

		DiffStepParams<_V, NV, _N, NN> operator+
			(
				const DiffStepParams<_V, NV, _N, NN>& p
			) const
		{
			DiffStepParams<_V, NV, _N, NN> result;
			result.m_vectors = this->m_vectors + p.vectors();
			result.m_scalars = this->m_scalars + p.scalars();
			return result;
		}

		template<class _tT>
		DiffStepParams<_V, NV, _N, NN> operator*
			(
				const _tT& h
			) const
		{
			DiffStepParams<_V, NV, _N, NN> result;
			result.m_scalars = this->m_scalars * h;
			result.m_vectors = this->m_vectors * h;
			return result;
		}

	};

	template<class _V, int NV, class _N, int NN, class _tT>
	DiffStepParams<_V, NV, _N, NN> operator*
		(
			const _tT& h,
			const DiffStepParams<_V, NV, _N, NN>& p
		)
	{
		return p*h;
	}
}

#endif // !_Integrators_
