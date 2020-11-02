
#include "stdafx.h"

#include "ParticleSource.h"
#include "ParticleTracking.h"
#include "Symmetry.hpp" // to be able to start particles in the virtually whole volume.
#include <algorithm>

namespace EvaporatingParticle
{


CSource::CSource()
{
  set_default();
}

CSource::~CSource()
{
  m_vData.clear();
}

void CSource::set_default()
{
  m_nType = stCone;
  m_nInjType = itHomogen;

  m_vPos = Vector3D(0, 0, 0);
  m_vDir = Vector3D(1, 0, 0);
  
  m_fRadius = 0.1;    // cm.
  m_fRingWidth = 0;   // cm, zero means a circular line (by default).

  m_fHeight = 0.06;   // cm, the letter-box slot height.
  m_fWidth = 0.2;     // cm, the letter-box slot width.

  m_fConeAngle = 1.0471975511965977461542144610932;
  m_fAbsVel = 100.;   // cm/s.

  m_nAzimCount = 0;
  m_nPitchCount = 0;
  m_nCount = 1 + m_nPitchCount * m_nAzimCount;
  m_nEnsembleSize = 1;
  m_nRandomSeed = 0;

  m_vFaceOffset = 0.002;  // cm, by default m_vFaceOffset is set to 0.05 mm.

  m_sFileName = CString("");

  m_bUseInitialGasVel = true;
  m_bReady = false;
}

void CSource::get(UINT nPartIndex,
                  Vector3D& vPos,
                  Vector3D& vVel,
                  double&   fTime,
                  double&   fPhase,
                  double&   fTemp,
                  double&   fIonMob,
                  UINT&     nEnsIndex,
                  size_t&   nElemId) const
{
  if(nPartIndex >= m_vData.size())
    return;

  const CPhasePoint& point = m_vData.at(nPartIndex);

  vPos = point.pos;
  vVel = point.vel;
  fTime = point.time;
  fPhase = point.phase;
  fTemp = point.temp;
  fIonMob = point.mob;
  nEnsIndex = point.ind;
  nElemId = point.elem;
}

bool CSource::gen_random_cone(double fPeriodRF, UINT nEnsSize)
{
  AfxMessageBox("Cone type of source with random distribution is not supported yet.");
  return false;
}

bool CSource::gen_homogen_cone(double fPeriodRF, UINT nEnsSize)
{
  CNode3D node;
  const CElem3D* pElem = NULL;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  Vector3D vPos = m_vPos, vVel, vAccel;
  CSymCorrData data = pObj->sym_corr_forward(vPos, vVel); // after this call vPos must be inside the ANSYS mesh.
// Check whether the modified position is really inside the simulation domain.
  if(!pObj->interpolate(vPos, 0, 0, node, pElem))
    return false;

  if(data.nSymFlag)
    pObj->sym_corr_back(vPos, vVel, vAccel, data);  // if position has been modified, get it back.

  UINT nEnsIndex = 0;
  vVel = m_fAbsVel * m_vDir;
  add_particle(vPos, vVel, node, fPeriodRF, pElem->nInd, nEnsSize, nEnsIndex);  // one particle starts exactly in m_vDir direction.
  nEnsIndex++;

  double fPitch, fAzim;
  double fPitchIncr = m_nPitchCount > 0 ? 0.5 * m_fConeAngle / m_nPitchCount : 0;
  double fAzimIncr = m_nAzimCount > 0 ? Const_2PI / m_nAzimCount : 0;
  Matrix3D mRotPitch, mRotAzim;

  UINT nMaxInd = m_nInjType == itRandom ? m_nCount : 1 + m_nPitchCount * m_nAzimCount;  // progress bar support.
  for(UINT i = 0; i < m_nPitchCount; i++)
  {
    if(pObj->get_terminate_flag())
      return terminate();

    fPitch = (i + 1) * fPitchIncr;
    mRotPitch = Matrix3D::rot(m_vLocY, fPitch);
    for(UINT j = 0; j < m_nAzimCount; j++)
    {
      fAzim = j * fAzimIncr;
      mRotAzim = Matrix3D::rot(m_vDir, fAzim);
      vVel = m_fAbsVel * (mRotAzim * (mRotPitch * m_vDir));
      add_particle(vPos, vVel, node, fPeriodRF, pElem->nInd, nEnsSize, nEnsIndex);
      nEnsIndex++;
// Progress bar support:
      pObj->set_progress(100 * (j + i * m_nAzimCount) / nMaxInd);
    }
  }

  return true;
}

bool CSource::gen_homogen_spot(double fPeriodRF, UINT nEnsSize)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  UINT nEnsIndex = 0;
  Vector3D vPos = m_vPos, vVel, vAccel, vShift;
  CSymCorrData data = pObj->sym_corr_forward(vPos, vVel);

  CNode3D node;
  const CElem3D* pElem = NULL;
  if(pObj->interpolate(vPos, 0, 0, node, pElem))
  {
    vVel = node.vel;
    if(data.nSymFlag)
      pObj->sym_corr_back(vPos, vVel, vAccel, data);

    if(!m_bUseInitialGasVel)
      vVel = m_fAbsVel * m_vDir;

    add_particle(vPos, vVel, node, fPeriodRF, pElem->nInd, nEnsSize, nEnsIndex);
    nEnsIndex++;
  }

  double fRad, fRadIncr = m_nPitchCount > 0 ? m_fRadius / m_nPitchCount : 0;
  double fAzim, fAzimIncr = m_nAzimCount > 0 ? Const_2PI / m_nAzimCount : 0;
  Matrix3D mRotAzim;
  UINT nMaxInd = 1 + m_nPitchCount * m_nAzimCount;
  for(UINT i = 0; i < m_nPitchCount; i++)
  {
    if(pObj->get_terminate_flag())
      return terminate();

    fRad = (i + 1) * fRadIncr;
    for(UINT j = 0; j < m_nAzimCount; j++)
    {
      fAzim = j * fAzimIncr;
      mRotAzim = Matrix3D::rot(m_vDir, fAzim);
      vShift = mRotAzim * m_vLocY;
      vPos = m_vPos + fRad * vShift;
      CSymCorrData data = pObj->sym_corr_forward(vPos, vVel);
      if(pObj->interpolate(vPos, 0, 0, node, pElem))
      {
        vVel = node.vel;
        if(data.nSymFlag)
          pObj->sym_corr_back(vPos, vVel, vAccel, data);

        if(!m_bUseInitialGasVel)
          vVel = m_fAbsVel * m_vDir;

        add_particle(vPos, vVel, node, fPeriodRF, pElem->nInd, nEnsSize, nEnsIndex);
        nEnsIndex++;
      }
// Progress bar support:
      pObj->set_progress(100 * (j + i * m_nAzimCount) / nMaxInd);
    }
  }

  return true;
}

bool CSource::gen_random_spot(double fPeriodRF, UINT nEnsSize)
{
  double fdX, fdY, fAmpl = 2 * m_fRadius;
  int nAttempt = 0;
  UINT nCount = m_nCount * nEnsSize;
  if(nCount == 0)
    return false;

  m_RndGen.seed(m_nRandomSeed);
  std::uniform_real_distribution<double> uniDistr(0.0, 1.0);

  Vector3D vPos = m_vPos, vVel, vAccel;
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  CNode3D node;
  const CElem3D* pElem = NULL;

  int nLoop = 1;
  UINT nEnsIndex = 0;
  while((m_vData.size() < nCount) && (nAttempt < RAND_MAX))
  {
    if((nLoop % 3 == 0) && pObj->get_terminate_flag())
      return terminate();

    fdX = -m_fRadius + fAmpl * uniDistr(m_RndGen);
    fdY = -m_fRadius + fAmpl * uniDistr(m_RndGen);

    vPos = m_vPos + fdX * m_vLocX + fdY * m_vLocY;
    if((vPos - m_vPos).length() > m_fRadius)
      continue;

    CSymCorrData data = pObj->sym_corr_forward(vPos, vVel);
    if(pObj->interpolate(vPos, 0, 0, node, pElem))
    {
      vVel = node.vel;
      if(data.nSymFlag)
        pObj->sym_corr_back(vPos, vVel, vAccel, data);

      if(!m_bUseInitialGasVel)
        vVel = m_fAbsVel * m_vDir;

      add_particle(vPos, vVel, node, fPeriodRF, pElem->nInd, nEnsSize, nEnsIndex);
      nEnsIndex++;
    }
// Progress bar support:
    pObj->set_progress(100 * m_vData.size() / m_nCount);
    nAttempt++;
    nLoop++;
  }

  return true;
}

bool CSource::gen_homogen_rect(double fPeriodRF, UINT nEnsSize)
{
  if(m_fHeight < Const_Almost_Zero || m_fWidth < Const_Almost_Zero)
    return false;

  Vector3D vPos, vVel, vAccel;
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  CNode3D node;
  const CElem3D* pElem = NULL;

  double fdX, fdY, fHalfWidth = 0.5 * m_fWidth, fHalfHeight = 0.5 * m_fHeight;
  UINT nCount = m_nCount * nEnsSize;

  double fRatio = m_fWidth / m_fHeight;
  double fNx = sqrt(fRatio * nCount);
  UINT Nx = UINT(0.5 + fNx);
  UINT Ny = UINT(0.5 + (double)nCount / fNx);

  double fOneOvrNx = Nx > 1 ? 1. / (Nx - 1) : 1;
  double fOneOvrNy = Ny > 1 ? 1. / (Ny - 1) : 1;
  UINT nEnsIndex = 0;
  for(UINT i = 0; i < Nx; i++)
  {
    if(pObj->get_terminate_flag())
      return terminate();

    fdX = -fHalfWidth + m_fWidth * i * fOneOvrNx;
    for(UINT j = 0; j < Ny; j++)
    {
      fdY = -fHalfHeight + m_fHeight * j * fOneOvrNy;

      vPos = m_vPos + fdX * m_vLocX + fdY * m_vLocY;
      CSymCorrData data = pObj->sym_corr_forward(vPos, vVel);
      if(pObj->interpolate(vPos, 0, 0, node, pElem))
      {
        vVel = node.vel;
        if(data.nSymFlag)
          pObj->sym_corr_back(vPos, vVel, vAccel, data);

        if(!m_bUseInitialGasVel)
          vVel = m_fAbsVel * m_vDir;

        add_particle(vPos, vVel, node, fPeriodRF, pElem->nInd, nEnsSize, nEnsIndex);
        nEnsIndex++;
      }
// Progress bar support:
      pObj->set_progress(100 * (j + i * Ny) / m_nCount);
    }
  }

  return true;
}

bool CSource::gen_random_rect(double fPeriodRF, UINT nEnsSize)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  double fdX, fdY, fHalfWidth = 0.5 * m_fWidth, fHalfHeight = 0.5 * m_fHeight;
  UINT nCount = m_nCount * nEnsSize;

  m_RndGen.seed(m_nRandomSeed);
  std::uniform_real_distribution<double> uniDistr(0.0, 1.0);

  CNode3D node;
  const CElem3D* pElem = NULL;
  int nLoop = 1, nAttempt = 0;
  Vector3D vPos, vVel, vAccel;
  UINT nEnsIndex = 0;

  while((m_vData.size() < nCount) && (nAttempt < RAND_MAX))
  {
    if((nLoop % 3 == 0) && pObj->get_terminate_flag())
      return terminate();

    fdX = -fHalfWidth + m_fWidth * uniDistr(m_RndGen);
    fdY = -fHalfHeight + m_fHeight * uniDistr(m_RndGen);

    vPos = m_vPos + fdX * m_vLocX + fdY * m_vLocY;
    CSymCorrData data = pObj->sym_corr_forward(vPos, vVel);
    if(pObj->interpolate(vPos, 0, 0, node, pElem))
    {
      vVel = node.vel;
      if(data.nSymFlag)
        pObj->sym_corr_back(vPos, vVel, vAccel, data);

      if(!m_bUseInitialGasVel)
        vVel = m_fAbsVel * m_vDir;

      add_particle(vPos, vVel, node, fPeriodRF, pElem->nInd, nEnsSize, nEnsIndex);
      nEnsIndex++;
    }
// Progress bar support:
    pObj->set_progress(100 * m_vData.size() / m_nCount);
    nAttempt++;
    nLoop++;
  }

  return true;
}

bool CSource::gen_homogen_cyl(double fPeriodRF, UINT nEnsSize)
{
  AfxMessageBox("Cylinder type of source with homogeneous distribution is not supported.");
  return false;
}

bool CSource::gen_random_cyl(double fPeriodRF, UINT nEnsSize)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  CNode3D node;
  const CElem3D* pElem = NULL;

  m_RndGen.seed(m_nRandomSeed);
  std::uniform_real_distribution<double> uniDistr(0.0, 1.0);

  int nLoop = 1;
  double fdPitch, fZ;
  Matrix3D mRotPitch;
  Vector3D vGlobZ(0, 0, 1), vLocDir, vPos, vVel, vAccel;

  UINT nEnsIndex = 0;
  UINT nCount = m_nCount * nEnsSize;
  while((m_vData.size() < nCount) && (nLoop < RAND_MAX))
  {
    if((nLoop % 3 == 0) && pObj->get_terminate_flag())
      return terminate();

    fdPitch = m_fConeAngle * (-0.5 + uniDistr(m_RndGen));
    mRotPitch = Matrix3D::rot(vGlobZ, fdPitch);
    fZ = m_fHeight * (-0.5 + uniDistr(m_RndGen));

    vLocDir = mRotPitch * m_vDir;
    vPos = m_vPos + m_fRadius * vLocDir + fZ * vGlobZ;
    CSymCorrData data = pObj->sym_corr_forward(vPos, vVel);
    if(pObj->interpolate(vPos, 0, 0, node, pElem))
    {
      vVel = node.vel;
      if(data.nSymFlag)
        pObj->sym_corr_back(vPos, vVel, vAccel, data);

      if(!m_bUseInitialGasVel)
        vVel = m_fAbsVel * m_vDir;

      add_particle(vPos, vVel, node, fPeriodRF, pElem->nInd, nEnsSize, nEnsIndex);
      nEnsIndex++;
    }
// Progress bar support:
    pObj->set_progress(100 * m_vData.size() / m_nCount);
    nLoop++;
  }

  return true;
}

bool CSource::gen_from_file(double fPeriodRF, UINT nEnsSize)
{
  std::vector<Vector3D> vPosColl;
  if(!read_pos_from_file(vPosColl))
    return false;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  CNode3D node;
  const CElem3D* pElem = NULL;
  Vector3D vPos, vVel, vAccel;

  UINT nEnsIndex = 0;
  m_nCount = vPosColl.size();
  for(UINT i = 0; i < m_nCount; i++)
  {
    vPos = vPosColl.at(i);
    CSymCorrData data = pObj->sym_corr_forward(vPos, vVel);
    if(pObj->interpolate(vPos, 0, 0, node, pElem))
    {
      vVel = node.vel;
      if(data.nSymFlag)
        pObj->sym_corr_back(vPos, vVel, vAccel, data);

      if(!m_bUseInitialGasVel)
        vVel = m_fAbsVel * m_vDir;

      add_particle(vPos, vVel, node, fPeriodRF, pElem->nInd, nEnsSize, nEnsIndex);
      nEnsIndex++;
    }
// Progress bar support:
    pObj->set_progress(100 * i / m_nCount);
  }

  return true;
}

bool CSource::read_pos_from_file(std::vector<Vector3D>& vPosColl)
{
  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (const char*)m_sFileName, (const char*)("r"));
  if(nErr != 0 || pStream == 0)
    return false;

  size_t nPntCount = 0;
  int nRes = fscanf_s(pStream, "%zd", &nPntCount);
  if(nRes == EOF || nRes == 0 || nPntCount == 0)
  {
    fclose(pStream);
    return false;
  }

  float x, y, z;
  Vector3D vNull(0, 0, 0);
  vPosColl.resize(nPntCount, vNull);
  for(size_t i = 0; i < nPntCount; i++)
  {
    nRes = fscanf_s(pStream, "%f, %f, %f", &x, &y, &z);
    if(nRes == EOF || nRes == 0)
    {
      fclose(pStream);
      return false;
    }

    vPosColl[i] = Vector3D(x, y, z);
  }

  fclose(pStream);
  return true;
}

bool CSource::generate_initial_cond()
{
  if(m_bReady)
    return true;  // if the m_vData array is filled and nothing has been changed since that time, use the old data.

  m_vData.clear();
  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->clear_selected_traject(); // the initial conditions have been changed, clear all track selections.

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  UINT nEnsSize = (pObj->get_particle_type() == CTrack::ptIon && !pObj->get_2D_flag())? m_nEnsembleSize : 1;
  double fPeriodRF = Const_2PI / pObj->get_rf_frequency();

  calc_loc_triad();

  switch(m_nType)
  {
    case stCone:
    {
      m_bReady = m_nInjType == itRandom ? gen_random_cone(fPeriodRF, nEnsSize) : gen_homogen_cone(fPeriodRF, nEnsSize);
      break;
    }
    case stSpot:
    {
      m_bReady = m_nInjType == itRandom ? gen_random_spot(fPeriodRF, nEnsSize) : gen_homogen_spot(fPeriodRF, nEnsSize);
      break;
    }
    case stRect:
    {
      m_bReady = m_nInjType == itRandom ? gen_random_rect(fPeriodRF, nEnsSize) : gen_homogen_rect(fPeriodRF, nEnsSize);
      break;
    }
    case stCylinder:
    {
      m_bReady = m_nInjType == itRandom ? gen_random_cyl(fPeriodRF, nEnsSize) : gen_homogen_cyl(fPeriodRF, nEnsSize);
      break;
    }
    case stSelReg:
    {
      m_RndGen.seed(m_nRandomSeed);
      m_bReady = populate_regions();
      break;
    }
    case stPosFromFile:
    {
      m_bReady = gen_from_file(fPeriodRF, nEnsSize);
      break;
    }
  }

  m_nCount = m_vData.size();
  return m_bReady;
}

bool CSource::terminate()
{
  m_vData.clear();
  m_bReady = false;
  return false;
}

void CSource::add_particle(const Vector3D& vPos,
                           const Vector3D& vVel,
                           const CNode3D&  node,
                           double          fPeriodRF,
                           size_t          nElemId,
                           UINT            nEnsSize,
                           UINT            nEnsIndex)
{
  double fStartTime, fPhase, fIonMob;
  double fCoeff = nEnsSize > 0 ? 1. / nEnsSize : 1;
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  std::uniform_real_distribution<double> uniDistr(0.0, 1.0);

  for(UINT i = 0; i < nEnsSize; i++)
  {
    fStartTime = double(i) * fCoeff * fPeriodRF;
    fPhase = Const_2PI * uniDistr(m_RndGen);
    fIonMob = pObj->get_ion_mob(node.press, node.temp);
    CPhasePoint pnt(vPos, vVel, fStartTime, fPhase, node.temp, fIonMob, nEnsIndex, nElemId);
    m_vData.push_back(pnt);
  }
}

void CSource::set_data_from_tracks()
{
  m_vData.clear();
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CTrackVector& vTracks = pObj->get_tracks();
  size_t nTrackCount = vTracks.size();
  for(size_t i = 0; i < nTrackCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    CBaseTrackItem* pItem = track.at(0);
    CPhasePoint pnt(pItem->pos, pItem->vel, pItem->time, track.get_phase(), pItem->get_temp(), track.get_index(), pItem->nElemId);
    m_vData.push_back(pnt);
  }

  m_bReady = true;
}

void CSource::calc_loc_triad()
{
  const Vector3D vGlobX(1, 0, 0), vGlobY(0, 1, 0), vGlobZ(0, 0, 1);
  Vector3D vLocY = m_vDir * vGlobX;
  double fLenY = vLocY.length();
  if(fLenY > Const_Almost_Zero)
  {
    m_vLocY = vLocY.normalized();
    m_vLocX = m_vLocY * m_vDir;
  }
  else
  {
    m_vLocX = vGlobY;
    m_vLocY = vGlobZ;
  }
}

const char* CSource::get_src_type_name(int nType) const
{
  switch(nType)
  {
    case stCone: return _T("Cone");
    case stSpot: return _T("Round Spot");
    case stRect: return _T("Rectangular Spot");
    case stCylinder: return _T("Cylinder (for 2D)");
    case stSelReg: return _T("Selected Regions");
    case stPosFromFile: return _T("Positions from File");
  }

  return _T("None");
}

const char* CSource::get_inj_type_name(int nType) const
{
  switch(nType)
  {
    case itRandom: return _T("Random");
    case itHomogen: return _T("Homogenious");
  }

  return _T("None");
}

// Selected region type of source.
bool CSource::populate_regions()
{
  if(!get_faces())
    return false;   // no faces.

  size_t nFaceCount = m_vFaces.size();
  double* pWeight = new double[nFaceCount];
  if(!calc_weights(pWeight))
  {
    delete[] pWeight;
    return false;   // the full area is almost zero.
  }

  if(m_nInjType == itHomogen)
  {
    AfxMessageBox("Homogeneous distribution at the selected regions is not supported.");
    return false;
  }

  bool bDoNotReflect = false;
  UINT nPointsCount = calc_points_count(bDoNotReflect);  // symmetry dependent.
  if(nPointsCount == 0)
    return false;

  UINT nCount;
  CFace* pFace = NULL;
  double fMeanCount, fRest = 0; // Normally, 0 < fRest < 1. If fRest becomes greater than 1, one more point must be added.
// Indices of triangles processed. If fRest > 1 the additional point will be settled in one of these triangles.
  std::vector<size_t> vTriInd;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    fMeanCount = pWeight[i] * nPointsCount;   // average numbers of points at every face.
    nCount = UINT(floor(fMeanCount));
    fRest += fMeanCount - nCount;

    pFace = m_vFaces.at(i);
    populate_face(pFace, nCount, !bDoNotReflect);
    vTriInd.push_back(i);

    if(fRest >= 1.)
    {
      UINT nFirst = vTriInd.at(0);
      UINT nLast = vTriInd.at(vTriInd.size() - 1);
      UINT nFace = choose_face(nFirst, nLast);    // randomly select a face from the processed faces ...
      if(nFace < nFaceCount)
        populate_face(m_vFaces.at(nFace), 1, !bDoNotReflect); // ... and add one more point to the selected face.

      vTriInd.clear();
      fRest -= 1.;
    }

    pObj->set_progress(100 * m_vData.size() / m_nCount);
    if(pObj->get_terminate_flag())
      return terminate();
  }

  delete[] pWeight;
  m_vFaces.clear();
  return true;
}

bool CSource::get_faces()
{
  m_vFaces.clear();
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CRegionsCollection vGlobReg = pObj->get_regions();
  size_t nRegCount = vGlobReg.size();
  for(size_t j = 0; j < nRegCount; j++)
  {
    CRegion* pReg = vGlobReg.at(j);
    if(std::find(m_vSelRegNames.begin(), m_vSelRegNames.end(), pReg->sName) != m_vSelRegNames.end())
    {
      size_t nFaceCount = pReg->vFaces.size();
      for(size_t i = 0; i < nFaceCount; i++)
        m_vFaces.push_back(pReg->vFaces.at(i));
    }
  }

  return m_vFaces.size() > 0;
}

bool CSource::calc_weights(double* pWeight) const
{
  double fArea, fFullArea = 0;
  size_t nFaceCount = m_vFaces.size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    CFace* pFace = m_vFaces.at(i);
    fArea = 0.5 * ((pFace->p1->pos - pFace->p0->pos) * (pFace->p2->pos - pFace->p0->pos)).length();
    pWeight[i] = fArea;
    fFullArea += fArea;
  }

  if(fFullArea < Const_Almost_Zero)
    return false;

  for(size_t i = 0; i < nFaceCount; i++)
    pWeight[i] /= fFullArea;

  return true;
}

UINT CSource::calc_points_count(bool& bDoNotReflect) const
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  int nSymmPlane = pObj->get_sym_plane();
  int nDiv = 1;
  if(nSymmPlane & CTracker::spXY)
    nDiv *= 2;
  if(nSymmPlane & CTracker::spXZ)
    nDiv *= 2;
  if(nSymmPlane & CTracker::spYZ)
    nDiv *= 2;

  bDoNotReflect = m_nCount < nDiv;
  if(bDoNotReflect)
    return m_nCount;

  return m_nCount / nDiv;
}

void CSource::populate_face(CFace* pFace, UINT nPntCount, bool bReflect)
{
  if(nPntCount == 0)
    return;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  CNode3D node;
  const CElem3D* pElem = NULL;
  double fAlpha, fBeta;
  Vector3D v0 = pFace->p0->pos;
  Vector3D e1 = pFace->p1->pos - v0;
  Vector3D e2 = pFace->p2->pos - v0;
  Vector3D vPos, vVel, vNorm = pFace->norm;

  std::uniform_real_distribution<double> uniDistr(0.0, 1.0);

  UINT nSuccess = 0;
  while(nSuccess < nPntCount)
  {
    fAlpha = uniDistr(m_RndGen);
    fBeta = uniDistr(m_RndGen);
    vPos = v0 + 0.5 * (fAlpha * e1 + fBeta * e2) - m_vFaceOffset * vNorm;

    if(pObj->interpolate(vPos, 0, 0, node, pElem))
    {
      vVel = m_bUseInitialGasVel ? node.vel : m_fAbsVel * m_vDir;
      add_particle(vPos, vVel, node, 0., pElem->nInd, 1, 0);
      nSuccess++;

      if(bReflect)
        reflect(vPos, vVel, node, pElem->nInd);
    }
  }
}

void CSource::reflect(const Vector3D& vPos, const Vector3D& vVel, const CNode3D& node, size_t nElemId)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  int nSymmPlane = pObj->get_sym_plane();

  Vector3D vReflPos = vPos, vReflVel = vVel;
  bool bSymXY = nSymmPlane & CTracker::spXY;
  if(bSymXY)
  {
    vReflPos.z = -vPos.z;
    vReflVel.z = -vVel.z;
    add_particle(vReflPos, vReflVel, node, 0., nElemId, 1, 0);
  }

  vReflPos = vPos;
  vReflVel = vVel;
  bool bSymXZ = nSymmPlane & CTracker::spXZ;
  if(bSymXZ)
  {
    vReflPos.y = -vPos.y;
    vReflVel.y = -vVel.y;
    add_particle(vReflPos, vReflVel, node, 0., nElemId, 1, 0);
  }

  vReflPos = vPos;
  vReflVel = vVel;
  if(bSymXY && bSymXZ)
  {
    vReflPos.z = -vPos.z;
    vReflPos.y = -vPos.y;
    vReflVel.z = -vVel.z;
    vReflVel.y = -vVel.y;
    add_particle(vReflPos, vReflVel, node, 0., nElemId, 1, 0);
  }
}

UINT CSource::choose_face(UINT nFirst, UINT nLast)
{
  UINT nCount = nLast - nFirst;
  std::uniform_real_distribution<double> uniDistr(0.0, 1.0);
  return nFirst + (UINT)(0.5 + nCount * uniDistr(m_RndGen));
}

// Streams support:
void CSource::save(CArchive& ar)
{
  const UINT nVersion = 5;  // 5 - m_sFileName, reading positions from a file support; 4 - m_vSelRegNames; 3 - m_fWidth and m_fHeight for the stRect type of the source.
  ar << nVersion;

  ar << m_nType;
  ar << m_nInjType;
  ar << m_nAzimCount;
  ar << m_nPitchCount;
  ar << m_nCount;
  ar << m_vPos.x;
  ar << m_vPos.y;
  ar << m_vPos.z;
  ar << m_vDir.x;
  ar << m_vDir.y;
  ar << m_vDir.z;
  ar << m_fRadius;
  ar << m_fRingWidth;
  ar << m_fConeAngle;
  ar << m_fAbsVel;
  ar << m_bUseInitialGasVel;
  ar << m_nEnsembleSize;
  ar << m_nRandomSeed;
  ar << m_fWidth;
  ar << m_fHeight;

  size_t nSelRegCount = m_vSelRegNames.size();
  ar << nSelRegCount;
  for(size_t j = 0; j < nSelRegCount; j++)
  {
    CString cName(m_vSelRegNames.at(j).c_str());
    ar << cName;
  }

  ar << m_sFileName;
}

void CSource::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_nType;
  ar >> m_nInjType;
  ar >> m_nAzimCount;
  ar >> m_nPitchCount;
  ar >> m_nCount;
  ar >> m_vPos.x;
  ar >> m_vPos.y;
  ar >> m_vPos.z;
  ar >> m_vDir.x;
  ar >> m_vDir.y;
  ar >> m_vDir.z;
  ar >> m_fRadius;
  ar >> m_fRingWidth;
  ar >> m_fConeAngle;
  ar >> m_fAbsVel;

  if(nVersion >= 1)
  {
    ar >> m_bUseInitialGasVel;
    ar >> m_nEnsembleSize;
  }

  if(nVersion >= 2)
  {
    ar >> m_nRandomSeed;
  }

  if(nVersion >= 3)
  {
    ar >> m_fWidth;
    ar >> m_fHeight;
  }

  if(nVersion >= 4)
  {
    size_t nSelRegCount;
    ar >> nSelRegCount;
    m_vSelRegNames.clear();
    for(size_t j = 0; j < nSelRegCount; j++)
    {
      CString cName;
      ar >> cName;
      std::string sName((const char*)cName);
      m_vSelRegNames.push_back(sName);
    }
  }

  if(nVersion >= 5)
  {
    ar >> m_sFileName;
  }

  m_bReady = false;
}

};