#pragma once
#ifndef _TRACKITEM_
#define _TRACKITEM_

#include <vector>
//#include "Symmetry.hpp"
#include "vector3d.hpp"
#include "../field_solver/MemoryPool.h"

namespace EvaporatingParticle
{

const ULONG ION_STATE_SIZE = 8;
const ULONG DROPLET_STATE_SIZE = 8;

struct CElem3D;
//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
struct CBaseTrackItem
{
  CBaseTrackItem()
    : nElemId(-1), time(0.)
  {
  }

  CBaseTrackItem(int nId, const Vector3D& p, const Vector3D& v, double t = 0.)
    : nElemId(nId), pos(p), vel(v), time(t)
  {
  }

  int         nElemId;  // mesh element which this item belongs to.

  Vector3D    pos,      // position and
              vel;      // velocity of the item.

  double      time;

  static CBaseTrackItem*  create(int nType, int nId, double fTime, double* pState);
  static CBaseTrackItem*  create(int nType);

  virtual CBaseTrackItem* copy() const = 0;
  virtual void            state(double* pState) const = 0;

  virtual double          get_mass() const = 0;
  virtual double          get_temp() const = 0;

  virtual void            save(CArchive& ar);
  virtual void            load(CArchive& ar);

  //[AC] Memory management
  static void operator delete(void* ptr, size_t n);
  virtual void deleteObj() = 0;
  //[/AC]
};

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
struct CIonTrackItem : public CBaseTrackItem, public BlockAllocator<CIonTrackItem>
{
	using CBaseTrackItem::operator delete;
  CIonTrackItem()
    : CBaseTrackItem(), temp(0), tempinf(0), unfragm(1)
  {
  }

  CIonTrackItem(int nId, const Vector3D& p, const Vector3D& v, double tmp, double unfrgm, double t = 0.)
    : CBaseTrackItem(nId, p, v, t),
      temp(tmp),
      tempinf(tmp),
      unfragm(unfrgm)
  {
  }

  double      temp,     // ion temperature.
              tempinf,  // steady-state (equilibrium) ion temperature.
              unfragm;  // part of unfragmented ions.

  virtual CBaseTrackItem* copy() const;
  virtual void            state(double* pState) const;

  virtual double          get_mass() const;
  virtual double          get_temp() const;

  virtual void            save(CArchive& ar);
  virtual void            load(CArchive& ar);

  //[AC] Memory management
  virtual void deleteObj();
  //[/AC]
};

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
struct CDropletTrackItem : public CBaseTrackItem, public BlockAllocator<CDropletTrackItem>
{
	using CBaseTrackItem::operator delete;
  CDropletTrackItem()
    : CBaseTrackItem(), temp(0), mass(0)
  {
  }

  CDropletTrackItem(int nId, const Vector3D& p, const Vector3D& v, double tmp, double m, double t = 0.)
    : CBaseTrackItem(nId, p, v, t),
      temp(tmp),
      mass(m)
  {
  }

  double      temp,   // evaporating droplet temperature.
              mass;   // variable droplet mass.

  virtual CBaseTrackItem* copy() const;
  virtual void            state(double* pState) const;

  virtual double          get_mass() const;
  virtual double          get_temp() const;

  virtual void            save(CArchive& ar);
  virtual void            load(CArchive& ar);

  //[AC] Memory management
  virtual void deleteObj();
  //[/AC]
};

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
struct CIntegrInterface
{
  CIntegrInterface(int nId, UINT nEnsInd, double phase, double curr)
    : nElemId(nId), nIndex(nEnsInd), fPhase(phase), fCurr(curr), fTionInf(0), bOk(true)
  {
  }

  int         nElemId;  // the mesh element which an item belongs to.
  UINT        nIndex;   // the ensemble index of an item (track specific).

  double      fPhase,   // the initial phase of a particle, either ion or droplet (track specific).
              fCurr,    // the current comprised inside the flux tube this ion belongs to (track specific).
              fTionInf; // the equilibrium ion temperature, (output).

  bool        bOk;      // a run-time flag intended to terminate a track (track-specific).
};

//---------------------------------------------------------------------------------------
// This is just a container to be filled by parameters of either CIonTrackItem or 
// CDropletTrackItem. Intended for usage in calculators and color contours.
//---------------------------------------------------------------------------------------
struct CTrackItem
{
  int         nElemId;  // mesh element which this item belongs to.

  Vector3D    pos,      // position and
              vel;      // velocity of the item.

  double      time,
              temp,     // ion temperature or droplet temperature.
              tempinf,  // steady-state ion temperature.
              unfragm,  // unfragmented part of ions.
              mass;     // variable droplet mass.
};

//---------------------------------------------------------------------------------------
// CTrack
//---------------------------------------------------------------------------------------
class CTrack : public std::vector<CBaseTrackItem*>
{
public:
  CTrack(int nType = ptIon, UINT nInd = 0, double fPhase = 0., double fCurr = 0.)
    : m_nType(nType), m_nIndex(nInd), m_fPhase(fPhase), m_fCurr(fCurr)
  {
  }

  enum  // Particle type
  {
    ptDroplet = 0,
    ptIon     = 1
  };

  int         get_type() const;
  void        set_type(int nType);

  UINT        get_index() const;
  void        set_index(UINT nIndex);

  double      get_phase() const;
  void        set_phase(double fPhase);

  double      get_current() const;
  void        set_current(double fCurr);

  void        get_track_item(size_t nInd, CTrackItem& item) const;  // no range control inside!

private:
  int         m_nType;
  UINT        m_nIndex;

  double      m_fPhase,
              m_fCurr;
};

typedef std::vector<CTrack> CTrackVector;

//---------------------------------------------------------------------------------------
// Inline implementation
//---------------------------------------------------------------------------------------
inline int CTrack::get_type() const
{
  return m_nType;
}

inline void CTrack::set_type(int nType)
{
  m_nType = nType;
}

inline UINT CTrack::get_index() const
{
  return m_nIndex;
}

inline void CTrack::set_index(UINT nIndex)
{
  m_nIndex = nIndex;
}

inline double CTrack::get_phase() const
{
  return m_fPhase;
}

inline void CTrack::set_phase(double fPhase)
{
  m_fPhase = fPhase;
}

inline double CTrack::get_current() const
{
  return m_fCurr;
}

inline void CTrack::set_current(double fCurr)
{
  m_fCurr = fCurr;
}

#if 0

struct CTrackItem
{
	CTrackItem(){}

	CTrackItem
	(
		const Vector3D& p, 
		const Vector3D& v, 
		double m, 
		double T, 
		double t = 0., 
		double phs = 0., 
		UINT ind = 0
	)
		:
		pos(p), vel(v), mass(m), 
		temp(T), time(t), phase(phs), 
		index(ind), curr(0)
	{
	}

  Vector3D  pos,
            vel;

  double    time,
            phase,  // phase now is a ion shift on a some period value

            mass,   // mass of the droplet.
            temp,   // temperature of the droplet.

            stationaryTemperature, //[AC] 22/07/2016

            charge,     // charge of the ion
            gasMolMass, // mass of the gas molecule
            kappa,      // ion mobility value
	      conc,       // unfragmented ion concentration	  

            diffCoef, // ion diffusion coefficient

            curr;   // full current which is comprised inside the current flux tube this item starts with and belongs to.

  UINT      index;  // index of the ensemble which this track belongs to.
};

	/**
	* Current Ansys mesh implementation.
	*/
	class AnsysMesh : public Mesh<CTracker, CElem3D, CNode3D>
	{
	public:
		AnsysMesh(const CTracker * mesh);
	};

	/**
	* AbsysField calculation
	*/
	class AnsysField : public IField
	{
		IMesh* m_mesh;
		const CTracker * m_tracker;
	public:
		AnsysField(const CTracker * tracker);

		~AnsysField();

		bool setNode(const V3D& pos);

		// [MS] 20-07-2016 added the ion velocity and current to support the Coulomb field in the axially-symmetric case.
		VV3D vectorFields(const V3D& pos, const V3D& vel, double time, double curr) const;


		Vd scalarFields(const V3D& pos, double time) const;

		/**
		 * Highest frequency of a field oscillation to adjust time step
		 */
		double fieldHighestFrequency() const;

		/**
		Get fragmentation energy
		*/
		double fragmentationEnergy() const;

		/**
		Get ion collision cross-section
		*/
		double ionCollisionCrossection() const;
	};

	/**
	 * Base class for different implementations of an ion moving
	 */
	/**
	* Switchs between high and low pressure cases
	*/
	#define LOW_PRESS  0UL
	#define HIGH_PRESS 1UL

	class BasePathIntegrator
	{
	protected:
		typedef math::DiffStepParams
			<
				Vector3D,	2,		//position and velocity
				double,		4		//temperature, diffusion 
									//[AC] 22/07/2016 and stationary temperature
									//[AC] 11/08/2016 ion's fragmentation probability
			> _DY;

		mutable bool m_isExist; //Checks wheither the particle is steel existing

		IField *m_field; //field interface
		CTrackItem m_item; //holds current item parameters
		math::IIntegrator* m_integrator;

		/**
		 * Calculate field values
		 */
		struct FieldValues
		{
			V3D E, U; //electrical field and gas flow
			double P, T, Cv; //pressure, temperature and gas heat capacitance
		};
		/**
		 * Calculates current field values
		 */
		FieldValues fieldValues(const V3D& pos, double time) const;

	public:

		/**
		 * Sets parameters for an differentiation step
		 */
		operator _DY() const;
		BasePathIntegrator& operator=(const _DY&);

		/**
		 * Inner item access
		 */
		const CTrackItem& item() const { return m_item; }
		CTrackItem& item() { return m_item; }

		/**
		 * Base construction
		 */
		BasePathIntegrator(const CTracker * tracker);

		virtual ~BasePathIntegrator();

		/**
		 * Calculates particle track saving system state every
		 * timeStepsPerSample stepof calculation
		 */
		virtual CTrack estimateTrack
		(
			int timeStepsPerSample,
			const CTrackItem& startConditions, //Start item conditions
			const bool * terminateFlag //If true then the calculation will be aborted
		) = 0;
	};

	class IonMobilityPathIntegrator : public BasePathIntegrator
	{
		math::ISwitchDiffScheme* m_switchDiffScheme; //Differentiation scheme

		/**
		 * Condition for the integrator switch between high and low pressure cases
		 */
		std::size_t switchIntegrator() const;

		//Estimates differentiatials at a current moment
		//for a high pressure case
		//[AC] 20/07/2016 one template parameter was added in order to switch between different ways
		//of ion mobility calculation
		template<class IonMobilityCalculator>
		_DY diffHP(const double& t0) const throw();

		//Estimates differentiatials at a current moment
		//for a low pressure case
		//[AC] 20/07/2016 one template parameter was added in order to switch between different ways
		//of ion mobility calculation
		template<class IonMobilityCalculator>
		_DY diffLP(const double& t0) const throw();

		/**
		 * Vector of functions that calculate derivatives
		 */
		typedef _DY(IonMobilityPathIntegrator::*dFun)(const double&) const throw();
		typedef std::vector<dFun> _DFuns;

		/**
		 * Integration step
		 */
		double m_h;

		/**
		* [AC] 19/07/2016 Time limit for an integration
		*/
		double m_maxTime;
	public:
		math::ISwitchDiffScheme* differScheme() { return m_switchDiffScheme; }
		double& h() { return m_h; }
		const double& h() const { return m_h; }
		double& maxTime() { return m_maxTime; }
		const double& maxTime() const { return m_maxTime; }

		IonMobilityPathIntegrator
		(
			const CTracker * tracker, //Mesh interface
			double h, //Time step of integration
			double maxTime //Limit time of integration
		);

		//Estimates ion track in a system
		virtual CTrack estimateTrack
		(
			int timeStepsPerSample,  //Number of time steps before item will be saved into track
			const CTrackItem& startConditions, //Start item conditions
			const bool * terminateFlag //If true then the calculation will be aborted
		);
	};

	/**
	 * Ion path integrator with an diffusion jumps after an each integration step
	 */
	class IonMobilityPathIntegratorDiffusionJumps 
		: public IonMobilityPathIntegrator
	{
	public:
		IonMobilityPathIntegratorDiffusionJumps
		(
			const CTracker * tracker,
			double h,
			double maxTime
		);

		virtual CTrack estimateTrack
		(
			int timeStepsPerSample,
			const CTrackItem& startConditions,
			const bool* terminateFlag
		);
	};

	template<class IonMobilityCalculator> BasePathIntegrator::_DY
	IonMobilityPathIntegrator::diffHP(const double & t0) const throw()
	{
		_DY diff;
		FieldValues fieldVals = fieldValues(item().pos, t0);
		if (!m_isExist)
			return diff;

		double mz = item().mass / item().charge;
		//current ion mobility
		IonMobilityCalculator kappaCalculator(item().kappa, fieldVals.T, fieldVals.P, item().temp);
		const double& kappa = kappaCalculator.ionMobility();
		//ion-gas relative velocity
		V3D vRel = item().vel - fieldVals.U;
		//time of relaxation
		double tau = kappa*mz;
		//ion-molecular reduced mass
		double mu = (item().mass*item().gasMolMass) / (item().mass + item().gasMolMass);
		//stationary temperature of an ion for a given relative ion-mol. velocity
		double Tinf = fieldVals.T + mu*(vRel&vRel) / 2 / fieldVals.Cv;
		//Stationary velocity
		V3D vinf = kappa*fieldVals.E + fieldVals.U;
		V3D dv = vinf - item().vel;
		double expH = (1. - std::exp(-m_h / tau));

		//Averaged over h step differential equations
		//It is supposed that E is a constant during the h interval	
		diff.vec3D<1>() = dv*expH / m_h; //dv
		diff.vec3D<2>() = vinf - tau*dv*expH / m_h; //dr

		diff.scalar<1>() = (Tinf - item().temp)*expH / m_h; //dT
															//diff coeff is proportional to a temperature
		diff.scalar<2>() = kappaCalculator.ionDiffusion(diff.scalar<1>(), fieldVals.T, item().charge);
		
		diff.scalar<3>() = (Tinf - item().stationaryTemperature) / m_h; //dTinf

		//[AC 11/08/2016] Calcultes ion's fragmentation probability
		//Activation energy of fragmentation
		double velocityFactor = sqrt(Const_Boltzmann*fieldVals.T / item().gasMolMass);
		double c_bar_gas =
			math::meanRelativeGasSpeed(vRel.length() / velocityFactor)*velocityFactor;
		//fragmentation constant
		double k = c_bar_gas * fieldVals.P / (Const_Boltzmann*fieldVals.T)
			* reinterpret_cast<AnsysField*>(m_field)->ionCollisionCrossection()
			* exp(-reinterpret_cast<AnsysField*>(m_field)->fragmentationEnergy() / (Const_Boltzmann * item().temp));
		diff.scalar<4>() = (exp(-k*m_h) - 1)*item().conc / m_h;

		return diff;
	}

	template<class IonMobilityCalculator> BasePathIntegrator::_DY
	IonMobilityPathIntegrator::diffLP(const double & t0) const throw()
	{
		//Check the particle is inside the system
		_DY diff;
		FieldValues fieldVals = fieldValues(item().pos, t0);
		if (!m_isExist) //Particle went out the system
			return diff;

		double mz = item().mass / item().charge;
		//current ion mobility
		IonMobilityCalculator kappaCalculator(item().kappa, fieldVals.T, fieldVals.P, item().temp);
		const double& kappa = kappaCalculator.ionMobility();
		//ion-gas relative velocity
		V3D vRel = item().vel - fieldVals.U;
		//time of relaxation
		double tau = kappa*mz;
		//ion-molecular reduced mass
		double mu = (item().mass*item().gasMolMass) / (item().mass + item().gasMolMass);
		//stationary temperature of an ion for a given relative ion-mol. velocity
		double Tinf = fieldVals.T + mu*(vRel&vRel) / 2 / fieldVals.Cv;

		diff.vec3D<1>() = fieldVals.E / mz - vRel / tau; //dv/dt
		diff.vec3D<2>() = item().vel; //dr/dt

		diff.scalar<1>() = (Tinf - item().temp) / tau; //dT/dt
		diff.scalar<2>() = kappaCalculator.ionDiffusion(diff.scalar<1>(), fieldVals.T, item().charge);

		diff.scalar<3>() = (Tinf - item().stationaryTemperature) / m_h; //dTinf

		//[AC 11/08/2016] Calcultes ion's fragmentation probability
	    //Activation energy of fragmentation
		double velocityFactor = sqrt(Const_Boltzmann*fieldVals.T / item().gasMolMass);
		double c_bar_gas =
			math::meanRelativeGasSpeed(vRel.length() / velocityFactor)*velocityFactor;
		//fragmentation constant
		double k = c_bar_gas * fieldVals.P / (Const_Boltzmann*fieldVals.T)
			* reinterpret_cast<AnsysField*>(m_field)->ionCollisionCrossection()
			* exp(-reinterpret_cast<AnsysField*>(m_field)->fragmentationEnergy() / (Const_Boltzmann * item().temp));
		diff.scalar<4>() = (exp(-k*m_h) - 1)*item().conc / m_h;

		return diff;
	}
#endif  // #if0

}

#endif // !_TRACKITEM_



