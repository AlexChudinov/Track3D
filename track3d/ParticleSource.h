
#pragma once

#include "CObject.h"
#include "constant.hpp"
#include "Vector3D.hpp"
#include "Elements.h"
#include "math.h"

#include <vector>
#include <string>

class CArchive;

namespace EvaporatingParticle
{

struct CPhasePoint
{
  CPhasePoint(const Vector3D& p, const Vector3D& v, double t = 0, double phs = 0., double tmp = 300., double b = 0., UINT id = 0, size_t ide = 0)
    : pos(p), vel(v), time(t), phase(phs), temp(tmp), mob(b), ind(id), elem(ide)
  {
  }

  Vector3D pos,
           vel;

  double   time,
           phase, // for RF phases randomization.
           temp,  // gas temperature.
           mob;   // ion mobility (random diffusion velocity jump support).

  UINT     ind;   // ensemble index.
  size_t   elem;  // index of the mesh element.
};

typedef std::vector<CPhasePoint> CSrcData;

//-------------------------------------------------------------------------------------------------
// CSource - a store for particles source parameters, it is used in setting the initial conditions.
//-------------------------------------------------------------------------------------------------
class CSource
{
public:
  CSource();
  ~CSource();

  enum  // type of source
  {
    stCone    = 0,  // point and full cone angle; the cone axis is m_vDir.
    stSpot    = 1,  // a circular spot centered at m_vPos with a normal direction m_vDir.
    stRect    = 2,  // a rectangular spot centered at m_vPos, m_vDir local X, Y dimensions equal to m_fWidth and m_fHeight, respectively.
    stRing    = 3,  // a ring centered at m_vPos perpendicular to the axis of the system.
    stSphere  = 4,  // in fact, this is a hemisphere.
    stSelReg  = 5,  // particles start from points randomly distributed on the selected 2D region(s).
    stCount   = 6
  };

  const char* get_src_type_name(int nType) const;

  enum  // type of injection
  {
    itRandom   = 0,
    itHomogen  = 1,
    itCount    = 2
  };

  const char*         get_inj_type_name(int nType) const;

  bool                generate_initial_cond();

  void                get(UINT      nPartIndex,
                          Vector3D& vPos,
                          Vector3D& vVel, 
                          double&   fTime,
                          double&   fPhase,
                          double&   fTemp,
                          double&   fMob,
                          UINT&     nEnsIndex,
                          size_t&   nElemId) const;

  int                 get_src_type() const;
  DWORD_PTR           get_src_type_ptr() const;
  void                set_src_type(int nType);

  int                 get_src_inject_type() const;
  DWORD_PTR           get_src_inject_type_ptr() const;
  void                set_src_inject_type(int nType);

  Vector3D            get_src_pos() const;
  DWORD_PTR           get_src_pos_ptr() const;
  void                set_src_pos(const Vector3D& vPos);

  Vector3D            get_inject_dir() const;
  DWORD_PTR           get_inject_dir_ptr() const;
  void                set_inject_dir(const Vector3D& vDir);

  double              get_cone_angle() const;
  DWORD_PTR           get_cone_angle_ptr() const;
  void                set_cone_angle(double fAngDeg);

  double              get_radius() const;
  DWORD_PTR           get_radius_ptr() const;
  void                set_radius(double fRadius);

  double              get_ring_width() const;
  DWORD_PTR           get_ring_width_ptr() const;
  void                set_ring_width(double fRingWidth);

  double              get_rect_width() const;
  DWORD_PTR           get_rect_width_ptr() const;
  void                set_rect_width(double fWidth);

  double              get_rect_height() const;
  DWORD_PTR           get_rect_height_ptr() const;
  void                set_rect_height(double fHeight);

  double              get_abs_vel() const;
  DWORD_PTR           get_abs_vel_ptr() const;
  void                set_abs_vel(double fAbsVel);

  UINT                get_particles_count() const;
  DWORD_PTR           get_particles_count_ptr() const;
  void                set_particles_count(UINT nCount);

  UINT                get_pitch_count() const;
  DWORD_PTR           get_pitch_count_ptr() const;
  void                set_pitch_count(UINT nCount);

  UINT                get_azim_count() const;
  DWORD_PTR           get_azim_count_ptr() const;
  void                set_azim_count(UINT nCount);

  bool                get_use_initial_gas_vel() const;
  DWORD_PTR           get_use_initial_gas_vel_ptr() const;
  void                set_use_initial_gas_vel(bool bEnable);

  UINT                get_ensemble_size() const;
  DWORD_PTR           get_ensemble_size_ptr() const;
  void                set_ensemble_size(UINT nSize);

  UINT                get_random_seed() const;
  DWORD_PTR           get_random_seed_ptr() const;
  void                set_random_seed(UINT nSeed);

  CString             get_selected_rgn_names() const;
  DWORD_PTR           get_selected_rgn_names_ptr() const;

// The tracks have been read from disk, set the initial conditions in accordance with these data.
  void                set_data_from_tracks(); 

  void                set_default();

// Streams support:
  void                save(CArchive& ar);
  void                load(CArchive& ar);

// Run-time
  void                invalidate();

protected:
  void                add_particle(const Vector3D& vPos, 
                                   const Vector3D& vVel,
                                   double          fGasTemp,
                                   double          fGasPress,
                                   double          fPeriodRF,
                                   size_t          nElemInd,
                                   UINT            nEnsSize, 
                                   UINT            nEnsIndex);

// A triad of mutually orthogonal axes in the source spot plane.
  void                calc_loc_triad();

// Selected region stuff. In contrast to the other types of source, the points are settled at a region lying within the
// domain. If the domain has symmetry, the points must be reflected. In the case of a single XY symmetry plane each point
// produces a miror point. In the case of two symmetry planes one single point results in four points.
  bool                populate_regions();

  bool                get_faces();

// The total count of points is m_nCount, but how many must be settled right on the selected region depends on the symmetry.
  UINT                calc_points_count(bool& bDoNotReflect) const; // if m_nCount is too small, apply no reflection. 

  bool                calc_weights(double* pWeights) const;   // pWeights is an array of m_vFaces.size() length.

  void                populate_face(CFace* pFace, UINT nPntCount, bool bReflect = true);
  void                reflect(const Vector3D& vPos, const Vector3D& vVel, double fTemp, double fPress, size_t nElemId);

  UINT                choose_face(UINT nFirst, UINT nLast) const;   // randomly select a number from a variety of nCount numbers.

// Termination
  bool                terminate();

  int                 m_nType,
                      m_nInjType;

  Vector3D            m_vPos,         // cm, position of the point source; in the case of line source this is the start point of the line.
                      m_vDir,         // unit vector of cone axis (spot normal) direction.
                      m_vLocX,
                      m_vLocY;

  double              m_fRadius,      // radius of emitter, cm.
                      m_fConeAngle,   // full angle of the cone, radian.
                      m_fRingWidth,   // the outer radius of the ring is m_fRadius + m_fRingWidth, cm.
                      m_fAbsVel,      // cm/s.
                      m_fHeight,      // height and width of the rectangular spot in local Y and X directions, respectively;
                      m_fWidth,       // (definition of local directions m_vLocX and m_vLocY see in calc_loc_triad()).

// In the case of selected region type the points are shifted a bit from the pFace plane in the direction -pFace->norm.
                      m_vFaceOffset;  // cm.

  UINT                m_nPitchCount,  // number of angular intervals in a plane where m_vDir lies.
                      m_nAzimCount,   // number of angular intervals in the azimuthal plane.
                      m_nCount,

                      m_nEnsembleSize,  // count of ions in every ensemble.

                      m_nRandomSeed;

  bool                m_bUseInitialGasVel;

  CSrcData            m_vData;
  CStringVector       m_vSelRegNames;
  CFacesCollection    m_vFaces; // a run-time variable, collection of all faces from all the selected regions.

private:
  bool                m_bReady; // a run-time variable to follow changes in the user-specified data.
};

inline int CSource::get_src_type() const
{
  return m_nType;
}

inline DWORD_PTR CSource::get_src_type_ptr() const
{
  return (DWORD_PTR)&m_nType;
}

inline void CSource::set_src_type(int nType)
{
  if(m_nType != nType)
  {
    m_nType = nType;
    m_bReady = false;
  }
}

inline int CSource::get_src_inject_type() const
{
  return m_nInjType;
}

inline DWORD_PTR CSource::get_src_inject_type_ptr() const
{
  return (DWORD_PTR)&m_nInjType;
}

inline void CSource::set_src_inject_type(int nType)
{
  if(m_nInjType != nType)
  {
    m_nInjType = nType;
    m_bReady = false;
  }
}

inline Vector3D CSource::get_src_pos() const
{
  return m_vPos;
}

inline DWORD_PTR CSource::get_src_pos_ptr() const
{
  return (DWORD_PTR)&m_vPos;
}

inline void CSource::set_src_pos(const Vector3D& vPos)
{
  if((m_vPos - vPos).length() > Const_Almost_Zero)
  {
    m_vPos = vPos;
    m_bReady = false;
  }
}

inline Vector3D CSource::get_inject_dir() const
{
  return m_vDir;
}

inline DWORD_PTR CSource::get_inject_dir_ptr() const
{
  return (DWORD_PTR)&m_vDir;
}

inline void CSource::set_inject_dir(const Vector3D& vDir)
{
  Vector3D vUnitDir = vDir.normalized();
  if((vUnitDir - m_vDir).length() > Const_Almost_Zero)
  {
    m_vDir = vUnitDir;
    m_bReady = false;
  }
}

inline double CSource::get_cone_angle() const
{
  return m_fConeAngle * Const_RadianToDegree;
}

inline DWORD_PTR CSource::get_cone_angle_ptr() const
{
  return (DWORD_PTR)&m_fConeAngle;
}

inline void CSource::set_cone_angle(double fAngDeg)
{
  double fAngRad = fAngDeg * Const_DegreeToRadian;
  if(fabs(m_fConeAngle - fAngRad) > Const_Almost_Zero)
  {
    m_fConeAngle = fAngRad;
    m_bReady = false;
  }
}

inline double CSource::get_radius() const
{
  return m_fRadius;
}

inline DWORD_PTR CSource::get_radius_ptr() const
{
  return (DWORD_PTR)&m_fRadius;
}

inline void CSource::set_radius(double fRadius)
{
  if(m_fRadius != fRadius)
  {
    m_fRadius = fRadius;
    m_bReady = false;
  }
}

inline double CSource::get_ring_width() const
{
  return m_fRingWidth;
}

inline DWORD_PTR CSource::get_ring_width_ptr() const
{
  return (DWORD_PTR)&m_fRingWidth;
}

inline void CSource::set_ring_width(double fRingWidth)
{
  if(m_fRingWidth != fRingWidth)
  {
    m_fRingWidth = fRingWidth;
    m_bReady = false;
  }
}

inline double CSource::get_rect_width() const
{
  return m_fWidth;
}

inline DWORD_PTR CSource::get_rect_width_ptr() const
{
  return (DWORD_PTR)&m_fWidth;
}

inline void CSource::set_rect_width(double fWidth)
{
  if(m_fWidth != fWidth)
  {
    m_fWidth = fWidth;
    m_bReady = false;
  }
}

inline double CSource::get_rect_height() const
{
  return m_fHeight;
}

inline DWORD_PTR CSource::get_rect_height_ptr() const
{
  return (DWORD_PTR)&m_fHeight;
}

inline void CSource::set_rect_height(double fHeight)
{
  if(m_fHeight != fHeight)
  {
    m_fHeight = fHeight;
    m_bReady = false;
  }
}

inline double CSource::get_abs_vel() const
{
  return m_fAbsVel;
}

inline DWORD_PTR CSource::get_abs_vel_ptr() const
{
  return (DWORD_PTR)&m_fAbsVel;
}

inline void CSource::set_abs_vel(double fAbsVel)
{
  if(m_fAbsVel != fAbsVel)
  {
    m_fAbsVel = fAbsVel;
    m_bReady = false;
  }
}

inline UINT CSource::get_particles_count() const
{
  return m_nCount;
}

inline DWORD_PTR CSource::get_particles_count_ptr() const
{
  return (DWORD_PTR)&m_nCount;
}

inline void CSource::set_particles_count(UINT nCount)
{
  if(m_nCount != nCount)
  {
    m_nCount = nCount;
    m_bReady = false;
  }
}

inline UINT CSource::get_pitch_count() const
{
  return m_nPitchCount;
}

inline DWORD_PTR CSource::get_pitch_count_ptr() const
{
  return (DWORD_PTR)&m_nPitchCount;
}

inline void CSource::set_pitch_count(UINT nCount)
{
  if(m_nPitchCount != nCount)
  {
    m_nPitchCount = nCount;
    m_bReady = false;
  }
}

inline UINT CSource::get_azim_count() const
{
  return m_nAzimCount;
}

inline DWORD_PTR CSource::get_azim_count_ptr() const
{
  return (DWORD_PTR)&m_nAzimCount;
}

inline void CSource::set_azim_count(UINT nCount)
{
  if(m_nAzimCount != nCount)
  {
    m_nAzimCount = nCount;
    m_bReady = false;
  }
}

inline bool CSource::get_use_initial_gas_vel() const
{
  return m_bUseInitialGasVel;
}

inline DWORD_PTR CSource::get_use_initial_gas_vel_ptr() const
{
  return (DWORD_PTR)&m_bUseInitialGasVel;
}

inline void CSource::set_use_initial_gas_vel(bool bEnable)
{
  if(m_bUseInitialGasVel != bEnable)
  {
    m_bUseInitialGasVel = bEnable;
    m_bReady = false;
  }
}

inline UINT CSource::get_ensemble_size() const
{
  return m_nEnsembleSize;
}

inline DWORD_PTR CSource::get_ensemble_size_ptr() const
{
  return (DWORD_PTR)&m_nEnsembleSize;
}

inline void CSource::set_ensemble_size(UINT nSize)
{
  if(m_nEnsembleSize != nSize)
  {
    m_nEnsembleSize = nSize;
    m_bReady = false;
  }
}

inline UINT CSource::get_random_seed() const
{
  return m_nRandomSeed;
}

inline DWORD_PTR CSource::get_random_seed_ptr() const
{
  return (DWORD_PTR)&m_nRandomSeed;
}

inline void CSource::set_random_seed(UINT nSeed)
{
  if(m_nRandomSeed != nSeed)
  {
    m_nRandomSeed = nSeed;
    m_bReady = false;
  }
}

inline CString CSource::get_selected_rgn_names() const
{
  return CObject::compile_string(m_vSelRegNames);
}

inline DWORD_PTR CSource::get_selected_rgn_names_ptr() const
{
  return (DWORD_PTR)&m_vSelRegNames;
}

inline void CSource::invalidate()
{
  m_bReady = false;
}

};  // EvaporatingParticle
