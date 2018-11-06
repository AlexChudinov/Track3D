
#include "stdafx.h"

#include "float.h"
#include "BeamCrossSection.h"
#include "ParticleTracking.h"
#include "constant.hpp"
#include "math.h"

#include "BarnesHut.h"

namespace EvaporatingParticle
{


CBeamCrossSection::CBeamCrossSection(double fX)
  : m_fPosX(fX)
{
  m_fSemiAxisA = 0;
  m_fSemiAxisB = 0;
  m_fAverVx = 0;
}

CBeamCrossSection::~CBeamCrossSection()
{
  m_vPoints.clear();
}

void CBeamCrossSection::init(const CTrackVector& vTracks)
{
  m_vPoints.clear();
  collect_intersections(vTracks);
  compute_semi_axes();
}

void CBeamCrossSection::collect_intersections(const CTrackVector& vTracks)
{
  double fKsi, fdX;
  Vector2D vPoint;
  m_fAverVx = 0;
  size_t nTrackCount = vTracks.size();
  for(size_t i = 0; i < nTrackCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    size_t nIonCount = track.size();

    Vector3D vPrev = track.at(0)->pos, vCurr;
    for(size_t j = 1; j < nIonCount; j++)
    {
      vCurr = track.at(j)->pos;
      if((vCurr.x >= m_fPosX) && (vPrev.x < m_fPosX) || (vCurr.x < m_fPosX) && (vPrev.x >= m_fPosX))
      {
        fdX = vCurr.x - vPrev.x;
        if(fabs(fdX) < Const_Almost_Zero)
          continue;

        fKsi = (m_fPosX - vPrev.x) / fdX;
        vPoint.x = vPrev.y + fKsi * (vCurr.y - vPrev.y);
        vPoint.y = vPrev.z + fKsi * (vCurr.z - vPrev.z);
        m_vPoints.push_back(vPoint);

        m_fAverVx += fabs(track.at(j)->vel.x);
        break;
      }

      vPrev = vCurr;
    }
  }

  size_t nPointsCount = m_vPoints.size();
  if(nPointsCount != 0)
    m_fAverVx /= nPointsCount;
}

void CBeamCrossSection::compute_semi_axes()
{
  size_t nPointsCount = m_vPoints.size();
  if(nPointsCount == 0)
    return;

  m_vCenter = Vector2D();
  for(size_t i = 0; i < nPointsCount; i++)
  {
    const Vector2D& point = m_vPoints.at(i);
    m_vCenter += point;
  }

  double fRecPointCount = 1./ nPointsCount;
  m_vCenter *= fRecPointCount;

// Corrections for different types of symmetry:
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  int nSymPlanes = pObj->get_sym_plane();
  int nSymXY = nSymPlanes & CTracker::spXY;
  int nSymXZ = nSymPlanes & CTracker::spXZ;

  if(nSymXY)
    m_vCenter.y = 0;
  if(nSymXZ)
    m_vCenter.x = 0;

  double fdX, fdY, fSigmX = 0, fSigmY = 0;
  double fAbsXMax = -FLT_MAX, fAbsYMax = -FLT_MAX;
  for(size_t i = 0; i < nPointsCount; i++)
  {
    const Vector2D& point = m_vPoints.at(i);
    fdX = fabs(point.x - m_vCenter.x);
    if(fdX > fAbsXMax)
      fAbsXMax = fdX;

    fdY = fabs(point.y - m_vCenter.y);
    if(fdY > fAbsYMax)
      fAbsYMax = fdY;

    fSigmX += fdX * fdX;
    fSigmY += fdY * fdY;
  }

  fSigmX *= fRecPointCount;
  m_fSemiAxisA = 2.0 * sqrt(fSigmX);
  if(m_fSemiAxisA > fAbsXMax)
    m_fSemiAxisA = fAbsXMax;    // make sure that the elliptical volume is not greater than the real size of the beam in X direction ...

  fSigmY *= fRecPointCount;
  m_fSemiAxisB = 2.0 * sqrt(fSigmY);
  if(m_fSemiAxisB > fAbsYMax)
    m_fSemiAxisB = fAbsYMax;    // ... and also take the same precaution for Y direction.
}


//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
CSpaceChargeVolume::CSpaceChargeVolume(CBeamCrossSection* pIn, CBeamCrossSection* pOut, double fTrackCurrent)
  : m_pPlaneIn(pIn), m_pPlaneOut(pOut)
{
  m_fMinX = m_pPlaneIn->m_fPosX;
  m_fMaxX = m_pPlaneOut->m_fPosX;

  m_fAverVx = 0.5 * (m_pPlaneIn->m_fAverVx + m_pPlaneOut->m_fAverVx);
  m_nTrackCount = (m_pPlaneIn->get_intersection_count() + m_pPlaneOut->get_intersection_count()) / 2;

  compute_full_charge(fTrackCurrent);
}

void CSpaceChargeVolume::collect_inner_vert(const CBox& box, double fMeshStep)
{
  double dx = m_fMinX - box.vMin.x;
  UINT nx0 = dx > Const_Almost_Zero ? 1 + (UINT)floor(dx / fMeshStep) : (UINT)floor(dx / fMeshStep);
  dx = m_fMaxX - box.vMin.x;
  UINT nx1 = (UINT)floor(dx / fMeshStep);

  UINT ny1 = (UINT)floor((box.vMax.y - box.vMin.y) / fMeshStep);
  UINT nz1 = (UINT)floor((box.vMax.z - box.vMin.z) / fMeshStep);

  double x, y, z, ksi, a, b, ra2, rb2, ry2, rz2, R;
  for(UINT i = nx0; i <= nx1; i++)
  {
    x = box.vMin.x + i * fMeshStep;
    ksi = (x - m_fMinX) / (m_fMaxX - m_fMinX);

    a = m_pPlaneIn->m_fSemiAxisA + ksi * (m_pPlaneOut->m_fSemiAxisA - m_pPlaneIn->m_fSemiAxisA);
    b = m_pPlaneIn->m_fSemiAxisB + ksi * (m_pPlaneOut->m_fSemiAxisB - m_pPlaneIn->m_fSemiAxisB);
    if(a < fMeshStep)
      a = fMeshStep;
    if(b < fMeshStep)
      b = fMeshStep;

    R = max(a, b);

    Vector2D c = m_pPlaneIn->m_vCenter + ksi * (m_pPlaneOut->m_vCenter - m_pPlaneIn->m_vCenter);

    ra2 = 1. / (a * a);
    rb2 = 1. / (b * b);

    for(UINT j = 0; j < ny1; j++)
    {
      y = box.vMin.y + j * fMeshStep - c.x;   // c.x is in fact the y coordinate in the global c.s.
      if(fabs(y) > R)
        continue;

      ry2 = y * y * ra2;  // (y/a)^2
      if(ry2 > 1)
        continue;

      for(UINT k = 0; k < nz1; k++)
      {
        z = box.vMin.z + k * fMeshStep - c.y; // c.y is in fact the z coordinate in the global c.s.
        if(fabs(z) > R)
          continue;

        rz2 = z * z * rb2;  // (z/b)^2
        if(rz2 > 1)
          continue;

        if(ry2 + rz2 > 1)
          continue;

        m_vInnerVert.push_back(Vector3D(x, y + c.x, z + c.y));
      }
    }
  }
}

void CSpaceChargeVolume::compute_full_charge(double fCurrPerTrack)
{
  m_fFullCharge = 0;
  if(m_fAverVx < Const_Almost_Zero)
    return;

  m_fFullCharge = fCurrPerTrack * m_nTrackCount * (m_fMaxX - m_fMinX) / m_fAverVx;
}

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
CSpaceChargeDistrib::CSpaceChargeDistrib()
  : m_pBHObject(NULL)
{
  set_default();
}

CSpaceChargeDistrib::~CSpaceChargeDistrib()
{
  clear();
}

void CSpaceChargeDistrib::set_default()
{
  m_nDistribType = distInCubicMesh;
  m_nPlanesCount = 100;

  m_fMeshStep = 0.01;   // cm.
  m_fTimeStep = 3e-7;   // s, 0.3 mcs.
}

bool CSpaceChargeDistrib::set_BH_object(CBarnesHut* pObj)
{
  m_pBHObject = pObj;
  if(m_pBHObject == NULL)
    return false;

  switch(m_nDistribType)
  {
    case distAlongTraject:
    {
      if(!set_charges_from_tracks())
        return false;
      break;
    }
    case distInCubicMesh:
    {
      init_planes();
      if(!init_volumes())
        return false;
      if(!collect_inner_vert())
        return false;
      if(!set_charges())
        return false;
      break;
    }
  }

  return true;
}

void CSpaceChargeDistrib::init_planes()
{
  if(m_nPlanesCount < 1)
    return;

  double fX, fdX = (m_fMaxX - m_fMinX) / (m_nPlanesCount - 1);
  for(UINT i = 0; i < m_nPlanesCount; i++)
  {
    fX = m_fMinX + i * fdX;
    if(i == 0)
      fX += 1e-6;

    m_vCrossSections.push_back(new CBeamCrossSection(fX));
  }
}

bool CSpaceChargeDistrib::init_volumes()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  if(pObj->get_terminate_flag())
    return false;

  set_job_name("Building Space Charge distribution volume...");
  set_progress(0);

  m_vElliptVolumes.clear();

  CTrackVector& vTracks = pObj->get_tracks();
  size_t nTrackCount = vTracks.size();
  if(nTrackCount < 1)
    return false;

  size_t nPlanesCount = m_vCrossSections.size();
  if(nPlanesCount < 1)
    return false;

  CBeamCrossSection* pSecIn = m_vCrossSections.at(0);
  pSecIn->init(vTracks);

  for(size_t i = 1; i < nPlanesCount; i++)
  {
    CBeamCrossSection* pSecOut = m_vCrossSections.at(i);
    pSecOut->init(vTracks);   // find intersections with tracks and build parameters of the cross-section's ellipse.

    CSpaceChargeVolume vol(pSecIn, pSecOut, m_fCurrPerTrack);
    m_vElliptVolumes.push_back(vol);

    pSecIn = pSecOut;

    if(pObj->get_terminate_flag())
      return false;

    set_progress(int(0.5 + 100. * (i + 1) / nPlanesCount));
  }

  set_progress(100);
  return true;
}

bool CSpaceChargeDistrib::collect_inner_vert()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  if(pObj->get_terminate_flag())
    return false;

  set_job_name("Collecting inner vertices...");
  set_progress(0);
  
  CBox& box = pObj->get_box();

  double fVolMinX, fVolMaxX;
  size_t nVolumeCount = m_vElliptVolumes.size();
  for(size_t j = 0; j < nVolumeCount; j++)
  {
    CSpaceChargeVolume& vol = m_vElliptVolumes.at(j);
    vol.collect_inner_vert(box, m_fMeshStep);
    if(pObj->get_terminate_flag())
      return false;

    set_progress(int(0.5 + 100. * (j + 1) / nVolumeCount));
  }

  set_progress(100);
  return true;
}

bool CSpaceChargeDistrib::set_charges()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  if(pObj->get_terminate_flag())
    return false;

  set_job_name("Setting pseudo-charges into inner vertices...");
  set_progress(0);

  Vector3D vPos;
  double fCharge = 0;
  size_t nVolumeCount = m_vElliptVolumes.size(), nInnerVertCount;
  for(size_t i = 0; i < nVolumeCount; i++)
  {
    if(pObj->get_terminate_flag())
      return false;

    set_progress(int(0.5 + 100. * (i + 1) / nVolumeCount));

    const CSpaceChargeVolume& vol = m_vElliptVolumes.at(i);
    fCharge = vol.get_charge_per_vertex();
    if(fCharge < Const_Almost_Zero)
      continue;

// Correction for symmetry:
    bool bCorrZ = false, bCorrY = false;
    int nSymType = m_pBHObject->get_sym_type();
    switch(nSymType)
    {
      case CBarnesHut::symXYonly: fCharge /= 2; bCorrZ = true; break;
      case CBarnesHut::symXZonly: fCharge /= 2; bCorrY = true; break;
      case CBarnesHut::symBoth: fCharge /= 4; bCorrZ = true; bCorrY = true; break;
    }

    nInnerVertCount = vol.get_inner_vert_count();
    for(size_t j = 0; j < nInnerVertCount; j++)
    {
      vPos = vol.get_inner_vertex(j);
// Correction for symmetry:
      if((vPos.z < 0) && bCorrZ)
        vPos.z = -vPos.z;
      if((vPos.y < 0) && bCorrY)
        vPos.y = -vPos.y;

      m_pBHObject->add_particle(vPos, fCharge);
    }
  }

  set_progress(100);
  return true;
}

void CSpaceChargeDistrib::clear()
{
  size_t nCount = m_vCrossSections.size();
  for(size_t i = 0; i < nCount; i++)
    delete m_vCrossSections.at(i);

  m_vCrossSections.clear();
  m_vElliptVolumes.clear();
  m_pBHObject = NULL;
}

const char* CSpaceChargeDistrib::get_distrib_type_name(int nDistType)
{
  switch(nDistType)
  {
    case CSpaceChargeDistrib::distAlongTraject: return _T("Along Trajectories");
    case CSpaceChargeDistrib::distInCubicMesh: return _T("In Cubical Mesh Nodes");
  }

  return _T("None");
}

Vector3D CSpaceChargeDistrib::radial_coulomb_field(const Vector3D& vPos, double fAxialVel, double fCurr)
{
  double r = sqrt(vPos.y * vPos.y + vPos.z * vPos.z);
  if(r < Const_Almost_Zero) // note that if r == 0, fCurr must be zero, too, as I = I0 * (r / r0)^2.
    return Vector3D();

  double fEr = 2 * fCurr / (r * fAxialVel);
  double fOne_ovr_R = 1. / r;
  return Vector3D(0, vPos.y * fOne_ovr_R * fEr, vPos.z * fOne_ovr_R * fEr);
}

Vector3D CSpaceChargeDistrib::radial_coulomb(const Vector3D& vPos, double fAxialVel) const
{
  int nInd = get_volume_index(vPos);
  if(nInd < 0)
    return Vector3D(0, 0, 0);

  const CSpaceChargeVolume& vol = m_vElliptVolumes.at(nInd);

  double fdX = vol.m_fMaxX - vol.m_fMinX;
  if(fdX < Const_Almost_Zero)
    return Vector3D(0, 0, 0);

  double fKsi = (vPos.x - vol.m_fMinX) / fdX;

  double fTrackCount = (1 - fKsi) * vol.m_pPlaneIn->get_intersection_count() + fKsi * vol.m_pPlaneOut->get_intersection_count();
  double fFullCurr = m_fCurrPerTrack * fTrackCount;   // full current flowing through this volume.

  Vector2D vC(0, 0);  // = (1 - fKsi) * vol.m_pPlaneIn->m_vCenter + fKsi * vol.m_pPlaneOut->m_vCenter;
  Vector3D vLocPos(0, vPos.y - vC.x, vPos.z - vC.y);  // local position with respect to the local x-axis of the volume.
  double fR2 = vLocPos.sqlength();
  if(fR2 < Const_Almost_Zero)
    return Vector3D(0, 0, 0); // the radial Coulomb field is always zero at the axis of the beam.

  double fA = (1 - fKsi) * vol.m_pPlaneIn->m_fSemiAxisA + fKsi * vol.m_pPlaneOut->m_fSemiAxisA;
  double fB = (1 - fKsi) * vol.m_pPlaneIn->m_fSemiAxisB + fKsi * vol.m_pPlaneOut->m_fSemiAxisB;
  double fMaxR = 0.5 * (fA + fB);   // the beam is expected to have almost circular cross-section.
  double fMaxR2 = fMaxR * fMaxR;

  double fCurr = fFullCurr;
  if(fR2 < fMaxR2)
    fCurr *= (fR2 / fMaxR2);  // current comprised inside the flux tube this trajectory belongs to.

  return radial_coulomb_field(vLocPos, fAxialVel, fCurr);
}

int CSpaceChargeDistrib::get_volume_index(const Vector3D& vPos) const
{
  size_t nVolCount = m_vElliptVolumes.size();
  if(nVolCount < 1 || vPos.x < m_fMinX || vPos.x > m_fMaxX)
    return -1;

  double fdX = (m_fMaxX - m_fMinX) / nVolCount;
  int nInd = int((vPos.x - m_fMinX) / fdX);
  if(nInd >= nVolCount)
    return -1;

  return nInd;
}

// Setting charges along the trajectories kitchen.
bool CSpaceChargeDistrib::set_charges_from_tracks()
{
  set_job_name("Setting pseudo-charges along the trajectories...");
  set_progress(0);

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CTrackVector& vTracks = pObj->get_tracks();
  size_t nTracksCount = vTracks.size();

// Attention! In the current version we use vTracks = pObj->get_tracks() to set ions along the trajectories. Note that
// vTracks contain items that are separated by N * pObj->get_time_step(), where N = nCoeff:
  int nCoeff = pObj->get_output_engine().get_output_time_step() > pObj->get_time_step() ? int(pObj->get_output_engine().get_output_time_step() / pObj->get_time_step()) : 1;

// The value of charge to set to the Barnes-Hut object, the current per track will be conserved.
  double fCharge = m_fCurrPerTrack * nCoeff * pObj->get_time_step();

// Correction for symmetry:
  bool bCorrZ = false, bCorrY = false;
  int nSymType = m_pBHObject->get_sym_type();
  switch(nSymType)
  {
    case CBarnesHut::symXYonly: fCharge /= 2; bCorrZ = true; break;
    case CBarnesHut::symXZonly: fCharge /= 2; bCorrY = true; break;
    case CBarnesHut::symBoth:   fCharge /= 4; bCorrZ = true; bCorrY = true; break;
  }

  Vector3D vPos;
  for(size_t i = 0; i < nTracksCount; i++)
  {
    if(pObj->get_terminate_flag())
      return false;

    set_progress(int(0.5 + 100. * (i + 1) / nTracksCount));

    const CTrack& track = vTracks.at(i);

// Randomization of pseudo-charge positions. Track-specific coefficient based on the track.get_phase():
    double fKsi = track.get_phase() / Const_2PI;  // must be in the range from 0 to 1.

    CBaseTrackItem* pItem = NULL;
    CBaseTrackItem* pNext = NULL;
    size_t nItemsCount = track.size();
    for(size_t j = 0; j < nItemsCount; j++)
    {
      pItem = track.at(j);
      if(j < nItemsCount - 1)
      {
        pNext = track.at(j + 1);
        vPos = pItem->pos * (1 - fKsi) + pNext->pos * fKsi;
      }
      else
      {
        continue;
      }

// Correction for symmetry:
      if((vPos.z < 0) && bCorrZ)
        vPos.z = -vPos.z;
      if((vPos.y < 0) && bCorrY)
        vPos.y = -vPos.y;

      m_pBHObject->add_particle(vPos, fCharge);
    }
  }

  return true;
}

void CSpaceChargeDistrib::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  ar << m_nDistribType;
  ar << m_nPlanesCount;
  ar << m_fMeshStep;
}

void CSpaceChargeDistrib::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_nDistribType;
  ar >> m_nPlanesCount;
  ar >> m_fMeshStep;
}

};  // namespace EvaporatingParticle