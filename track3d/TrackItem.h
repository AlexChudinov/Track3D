#pragma once
#ifndef _TRACKITEM_
#define _TRACKITEM_

#include <vector>
#include "vector3d.hpp"

//[AC] memory manager
#include <../utilities/MemoryPool.h>
//[/AC]

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
  virtual void            set_param(int nId, double fTime, double* pState) = 0;

  virtual double          get_mass() const = 0;
  virtual double          get_temp() const = 0;
  virtual double          get_ion_mob() const = 0;

  virtual void            save(CArchive& ar);
  virtual void            load(CArchive& ar);

  //[AC] memory manager
  static void operator delete(void* ptr, size_t n);
  virtual void deleteObj() = 0;
  //[/AC]
};

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
struct CIonTrackItem : public CBaseTrackItem, public BlockAllocator<CIonTrackItem>
{
	//[AC] memory manager
	using CBaseTrackItem::operator delete;
	//[/AC]
  CIonTrackItem()
    : CBaseTrackItem(), temp(0), tempinf(0), unfragm(1), mob(0)
  {
  }

  CIonTrackItem(int nId, const Vector3D& p, const Vector3D& v, double tmp, double unfrgm, double b, double t)
    : CBaseTrackItem(nId, p, v, t),
      temp(tmp),
      tempinf(tmp),
      unfragm(unfrgm),
      mob(b)
  {
  }

  double      temp,     // ion temperature.
              tempinf,  // steady-state (equilibrium) ion temperature.
              unfragm,  // part of unfragmented ions.
              mob;      // ion mobility, random diffusion velocity jumps support.

  virtual CBaseTrackItem* copy() const;
  virtual void            state(double* pState) const;
  virtual void            set_param(int nId, double fTime, double* pState);

  virtual double          get_mass() const;
  virtual double          get_temp() const;
  virtual double          get_ion_mob() const;

  virtual void            save(CArchive& ar);
  virtual void            load(CArchive& ar);

  virtual void deleteObj();
};

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
struct CDropletTrackItem : public CBaseTrackItem, public BlockAllocator<CDropletTrackItem>
{
	//[AC] memory manager
	using CBaseTrackItem::operator delete;
	//[/AC]
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
  virtual void            set_param(int nId, double fTime, double* pState);

  virtual double          get_mass() const;
  virtual double          get_temp() const;
  virtual double          get_ion_mob() const;

  virtual void            save(CArchive& ar);
  virtual void            load(CArchive& ar);

  virtual void deleteObj();
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
              fTionInf, // the equilibrium ion temperature, (output).
              fIonMob;  // random diffusion velocity jumps support.

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
              mass,     // variable droplet mass.
              mob;      // ion mobility at recalculated for current T and P.
};

//---------------------------------------------------------------------------------------
// CTrack
//---------------------------------------------------------------------------------------
class CTrack : public std::vector<CBaseTrackItem*>
{
public:
  CTrack(int nType = ptIon, UINT nInd = 0, double fPhase = 0., double fCurr = 0.)
    : m_nType(nType), m_nIndex(nInd), m_fPhase(fPhase), m_fCurr(fCurr), m_nRandSeed(0), m_nTermReason(ttrNone)
  {
  }

  enum  // Particle type
  {
    ptDroplet = 0,
    ptIon     = 1
  };

  enum // Track termination reason
  {
    ttrNone         = 0,
    ttrOvrTime      = 1,
    ttrRayleighLim  = 2,
    ttrFullyEvapor  = 3
  };

  int         get_type() const;
  void        set_type(int nType);

  UINT        get_index() const;
  void        set_index(UINT nIndex);

  double      get_phase() const;
  void        set_phase(double fPhase);

  double      get_current() const;
  void        set_current(double fCurr);

  UINT        get_rand_seed() const;
  void        set_rand_seed(UINT nSeed);

  int&        get_term_reason();
  int         get_term_reason() const;

  void        get_track_item(size_t nInd, CTrackItem& item) const;  // no range control inside!

private:
  int         m_nType;
  UINT        m_nIndex;

  double      m_fPhase,
              m_fCurr;

  UINT        m_nRandSeed;    // random diffusion support.
  int         m_nTermReason;  // track termination reason.
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

inline UINT CTrack::get_rand_seed() const
{
  return m_nRandSeed;
}

inline void CTrack::set_rand_seed(UINT nSeed)
{
  m_nRandSeed = nSeed;
}

inline int& CTrack::get_term_reason()
{
  return m_nTermReason;
}

inline int CTrack::get_term_reason() const
{
  return m_nTermReason;
}

} //namespace EvaporatingParticle

namespace EvaporatingParticle
{
	class CTracker;
}
namespace _new_version
{
	using EvaporatingParticle::Vector3D;
	using EvaporatingParticle::CTracker;

	enum TrackItemType
	{
		Droplet,
		Ion
	};

	struct TrackItemParams
	{
		virtual TrackItemType type() const = 0;

		//Particle position
		Vector3D mPos;
		//Particle velocity
		Vector3D mVel;
		//Calculation time
		double mTime;
	};

	class CBaseTrackItem
	{
	public:
		using ParamsPtr = std::unique_ptr<TrackItemParams>;
		using ItemPtr = std::unique_ptr<CBaseTrackItem>;

		static ItemPtr create(const TrackItemParams * params);

		virtual void do_step(CTracker * tracker) = 0;

		virtual ItemPtr clone() const = 0;

		static void operator delete(void* ptr, size_t n);
		virtual void deleteObj() = 0;
	};

	class CDropletTrackItem : public CBaseTrackItem, public ::BlockAllocator<CDropletTrackItem>
	{
	public:
		using CBaseTrackItem::operator delete;

		struct DropletParams : public TrackItemParams
		{
			virtual TrackItemType type() const;
		};

		CDropletTrackItem(const TrackItemParams * params);

		virtual void do_step(CTracker * tracker);

		virtual ItemPtr clone() const;

		virtual void deleteObj();

	private:

		DropletParams mParams;
	};
} //namespace new version

#endif // !_TRACKITEM_



