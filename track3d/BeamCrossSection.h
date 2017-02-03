#pragma once

#include "Elements.h"
#include "TrackItem.h"
#include "vector2d.hpp"
#include "CObject.h"

namespace EvaporatingParticle
{

typedef std::vector<Vector2D> CVector2DVec;
//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
class CBeamCrossSection
{
public:
  CBeamCrossSection(double fX);
  virtual ~CBeamCrossSection();

  void          init(const CTrackVector& vTracks);

  size_t        get_intersection_count() const;
  Vector2D      get_intersection(size_t nInd) const;

protected:
  void          collect_intersections(const CTrackVector& vTracks);
  void          compute_semi_axes();

public:
  double        m_fPosX,      // x-position of the cross-section.

                m_fSemiAxisA, // along local x (global y),
                m_fSemiAxisB; // along local y (global z).

  Vector2D      m_vCenter;    // average (y,z) of the points-intersections.
  double        m_fAverVx;    // average absolute value |Vx| of x-component of velocity for computation of the whole space charge.

private:
  CVector2DVec  m_vPoints;    // N = m_vPoints.size() is the count of trajectories intersecting the plane.
};

typedef std::vector<Vector3D> CInnerVertVec;
//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
class CSpaceChargeVolume
{
public:
  CSpaceChargeVolume(CBeamCrossSection* pIn, CBeamCrossSection* pOut, double fTrackCurrent);

  void                collect_inner_vert(const CBox& box, double fMeshStep);

  UINT                get_inner_vert_count() const;
  Vector3D            get_inner_vertex(UINT nInd) const;

  double              get_charge_per_vertex() const;

protected:
  void                compute_full_charge(double fCurrPerTrack);

  bool                inside(const Vector3D& vPos) const;

private:
  double              m_fFullCharge;

  double              m_fMinX,  // x-positions of the left and right ends of the elliptical cylindrical volume.
                      m_fMaxX;

  double              m_fAverVx;    // average absolute value |Vx| of x-component of velocity for computation of the whole space charge.

  UINT                m_nTrackCount;  // count of trajectories inside the volume.

  CBeamCrossSection*  m_pPlaneIn;
  CBeamCrossSection*  m_pPlaneOut;

// Vertices of the domain which are inside this volume. I need their full count because I only know the full charge in the volume.
  CInnerVertVec       m_vInnerVert;

  friend class CSpaceChargeDistrib;
};

class CBarnesHut;
typedef std::vector<CBeamCrossSection*> CBeamCSColl;
typedef std::vector<CSpaceChargeVolume> CElliptVolVec;
//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
class CSpaceChargeDistrib : public CObject
{
public:
  CSpaceChargeDistrib();
  virtual ~CSpaceChargeDistrib();

  enum // Pseudo-ions distribution type
  {
    distAlongTraject = 0,   // the ions are distributed along the ion trajectories of the previous iteration ...
    distInCubicMesh  = 1,   // ... or in the nodes of some constant virtual cubic mesh.
    distCount        = 2
  };

  int                 get_ion_distrib_type() const;
  DWORD_PTR           get_ion_distrib_type_ptr() const;
  void                set_ion_distrib_type(int nType);

  UINT                get_planes_count() const;       // a count of x = const planes which subdivide the domain into finite
  DWORD_PTR           get_planes_count_ptr() const;   // elliptical volumes using the CSpaceChargeDistrib class; the pseudo-charges
  void                set_planes_count(UINT nCount);  // are settled with a constant density inside the volumes.

  double              get_space_charge_step() const;
  DWORD_PTR           get_space_charge_step_ptr() const;
  void                set_space_charge_step(double fStep);

  Vector3D            radial_coulomb(const Vector3D& vPos, double fAxialVel) const;

  static const char*  get_distrib_type_name(int nDistType);

// Gabovich M.D. Influence of the space charge on the intense beams of charged particles propagation. 
// UFN vol. 56, issue 2, pp 215 - 256, June, 1955.
  static Vector3D     radial_coulomb_field(const Vector3D& vPos, double fAxialVel, double fCurr);

  void                save(CArchive& archive);
  void                load(CArchive& archive);

// Run-time, called from CTracker::create_BH_object(...):
  void                set_run_time_data(double fMinX, double fMaxX, double fCurrPerTrack);
  bool                set_BH_object(CBarnesHut* pObj);
  void                clear();

protected:
  void                set_default();

  int                 get_volume_index(const Vector3D& vPos) const;

// Interface for distribution in the cubic mesh nodes:
  void                init_planes();    // called only from constructor.

  bool                init_volumes();
  bool                collect_inner_vert();
  bool                set_charges();

// Interface for distribution along trajectories taken from the previous iteration:
  bool                set_charges_from_tracks();

private:
  int                 m_nDistribType;
  UINT                m_nPlanesCount; // count of x = const planes subdividing the domain into elliptical volumes.

  double              m_fMeshStep;    // the pseudo-charges are placed into the nodes of a constant cubical mesh.

  double              m_fMinX,
                      m_fMaxX,

                      m_fCurrPerTrack;  // current, CGSE, which every track carries; the iteration and symmetry are taken into account.

  CBeamCSColl         m_vCrossSections;
  CElliptVolVec       m_vElliptVolumes;
  CBarnesHut*         m_pBHObject;
};


//---------------------------------------------------------------------------------------
// Inline implementation
//---------------------------------------------------------------------------------------
inline size_t CBeamCrossSection::get_intersection_count() const
{
  return m_vPoints.size();
}

inline Vector2D CBeamCrossSection::get_intersection(size_t nInd) const
{
  return m_vPoints.at(nInd);
}

inline UINT CSpaceChargeVolume::get_inner_vert_count() const
{
  return m_vInnerVert.size();
}

inline Vector3D CSpaceChargeVolume::get_inner_vertex(UINT nInd) const
{
  size_t nCount = m_vInnerVert.size();
  return nInd < nCount ? m_vInnerVert.at(nInd) : Vector3D();
}

inline double CSpaceChargeVolume::get_charge_per_vertex() const
{
  size_t nCount = m_vInnerVert.size();
  return nCount > 0 ? m_fFullCharge / nCount : 0;
}

inline int CSpaceChargeDistrib::get_ion_distrib_type() const
{
  return m_nDistribType;
}

inline DWORD_PTR CSpaceChargeDistrib::get_ion_distrib_type_ptr() const
{
  return (DWORD_PTR)&m_nDistribType;
}

inline void CSpaceChargeDistrib::set_ion_distrib_type(int nType)
{
  m_nDistribType = nType;
}

inline UINT CSpaceChargeDistrib::get_planes_count() const
{
  return m_nPlanesCount;
}

inline DWORD_PTR CSpaceChargeDistrib::get_planes_count_ptr() const
{
  return (DWORD_PTR)&m_nPlanesCount;
}

inline void CSpaceChargeDistrib::set_planes_count(UINT nCount)
{
  m_nPlanesCount = nCount;
}

inline double CSpaceChargeDistrib::get_space_charge_step() const
{
  return m_fMeshStep;
}

inline DWORD_PTR CSpaceChargeDistrib::get_space_charge_step_ptr() const
{
  return (DWORD_PTR)&m_fMeshStep;
}

inline void CSpaceChargeDistrib::set_space_charge_step(double fStep)
{
  m_fMeshStep = fStep;
}

inline void CSpaceChargeDistrib::set_run_time_data(double fMinX, double fMaxX, double fCurrPerTrack)
{
  m_fMinX = fMinX;
  m_fMaxX = fMaxX;
  m_fCurrPerTrack = fCurrPerTrack;
}

};  // namespace EvaporatingParticle