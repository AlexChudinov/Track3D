
#include "stdafx.h"

#include "float.h"
#include "Calculator.h"

#include "Tracker.hpp"
#include "ParticleTracking.h"
#include "mathematics.h"

#include "..\utilities\ParallelFor.h"

#include <algorithm>


namespace EvaporatingParticle
{

static const double scfGasConstantR = Const_Boltzmann * Const_Avogadro;

//-------------------------------------------------------------------------------------------------
// CCalcCollection.
//-------------------------------------------------------------------------------------------------
CCalcCollection::~CCalcCollection()
{
  clear_calcs();
}

void CCalcCollection::clear_calcs()
{
  for(size_t i = 0; i < size(); i++)
  {
    CCalculator* pObj = at(i);
    delete pObj;
  }

  clear();
}

void CCalcCollection::calculate()
{
  for(size_t i = 0; i < size(); i++)
  {
    CCalculator* pObj = at(i);
    pObj->do_calculate();
  }
}

bool CCalcCollection::sel_region_changed(CStringVector* pRegNames)
{
  size_t nCalcCount = size();
  for(size_t i = 0; i < nCalcCount; i++)
  {
    CCalculator* pCalc = at(i);
    if(pCalc->type() == CCalculator::ctSelRegions)
    {
      CSelectedRegionCalculator* pSelRegCalc = (CSelectedRegionCalculator*)pCalc;
      if(pSelRegCalc->get_sel_reg_names_ptr() == (DWORD_PTR)pRegNames)
      {
        pSelRegCalc->do_calculate();
        return true;
      }
    }
  }

  return false;
}

void CCalcCollection::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  size_t nObjCount = size();
  ar << nObjCount;
  for(size_t i = 0; i < nObjCount; i++)
  {
    CCalculator* pObj = at(i);
    int nType = pObj->type();
    ar << nType;
    pObj->save(ar);
  }
}

void CCalcCollection::load(CArchive& ar)
{
  clear_calcs();

  UINT nVersion;
  ar >> nVersion;

  size_t nObjCount;
  ar >> nObjCount;
  for(size_t i = 0; i < nObjCount; i++)
  {
    int nType;
    ar >> nType;
    CCalculator* pObj = CCalculator::create(nType);
    pObj->load(ar);
    push_back(pObj);
  }
}

//-------------------------------------------------------------------------------------------------
// CCalculator - the base class designed for calculations on arbitrary scalar or vector fields.
//-------------------------------------------------------------------------------------------------
CCalculator::CCalculator()
  : m_pObj(NULL)
{
  m_fResult = 0;
  m_nClcVar = clcAveTemp;
  m_bClcVarChanged = false;
  m_fCharLength = 0.06;   // d = 0.6 mm by default for Reynolds number calculation in capillary.
  m_nCalcDir = dirX;      // by default the sequence calculations are performed in X direction.
  m_bEnable = true;
}

CCalculator::~CCalculator()
{
}

CCalculator* CCalculator::create(int nType)
{
  switch(nType)
  {
    case ctPlaneYZ: return new CPlaneYZCalculator();
    case ctSelRegions: return new CSelectedRegionCalculator();
    case ctAlongLine: return new CLineCalculator();
    case ctTrackCalc: return new CTrackCalculator();
    case ctTrackCrossSect: return new CTrackCrossSectionCalculator();
    case ctAlongSelTracks: return new CSelectedTracksCalculator();
    case ctTrackFaceCross: return new CTackFaceCross();
    case ctEndTrackCalc: return new CEndTrackCalculator();
  }

  return NULL;
}

const char* CCalculator::calc_name(int nType)
{
  switch(nType)
  {
    case ctPlaneYZ: return _T("Calculator at a Domain Cross-Section");
    case ctSelRegions: return _T("Calculator at Selected Regions");
    case ctAlongLine: return _T("Line Calculator");
    case ctTrackCalc: return _T("Ion/Droplet Track Calculator");
    case ctTrackCrossSect: return _T("Track's Cross-Section Calculator");
    case ctAlongSelTracks: return _T("Forces at Selected Tracks");
    case ctTrackFaceCross: return _T("Faces Crossed by Tracks Calculator");
    case ctEndTrackCalc: return _T("Coordinates of End Points of Tracks");
  }

  return _T("");
}

const char* CCalculator::plane_norm_name(int nDir)
{
  switch(nDir)
  {
    case CCalculator::dirX: return _T("X");
    case CCalculator::dirY: return _T("Y");
    case CCalculator::dirZ: return _T("Z");
  }

  return _T("");
}

bool CCalculator::get_tracker_ptr()
{
  m_pObj = CParticleTrackingApp::Get()->GetTracker();
  if(m_pObj == NULL || !m_pObj->is_ready())
    return false;

  return true;
}

bool CCalculator::intersect(const CRay& ray, CIntersectColl& results) const
{
  results.clear();
  const CRegionsCollection& regions = m_pObj->get_regions();
  size_t nRegCount = regions.size();

  double fDist;
  for(size_t j = 0; j < nRegCount; j++)
  {
    CRegion* pReg = regions.at(j);
    pReg->intersect(ray, results);
  }

  bool bOk = results.size() > 0;
  if(bOk)
    std::sort(results.begin(), results.end());

  return bOk;
}

bool CCalculator::process_face(CFace* pFace, double& fSquare, double& fRes) const
{
  Vector3D e1 = pFace->p1->pos - pFace->p0->pos;
  Vector3D e2 = pFace->p2->pos - pFace->p0->pos;
  fSquare = 0.5 * (e1 * e2).length();
  if(fSquare < Const_Almost_Zero)
    return false;

  double fx, fy, fz;
  switch(m_nClcVar)
  {
    case clcMassFlow:
    case clcVolumeFlow:
    {
      fx = (pFace->p0->vel.x * pFace->p0->dens + pFace->p1->vel.x * pFace->p1->dens + pFace->p2->vel.x * pFace->p2->dens) * pFace->norm.x;
      fy = (pFace->p0->vel.y * pFace->p0->dens + pFace->p1->vel.y * pFace->p1->dens + pFace->p2->vel.y * pFace->p2->dens) * pFace->norm.y;
      fz = (pFace->p0->vel.z * pFace->p0->dens + pFace->p1->vel.z * pFace->p1->dens + pFace->p2->vel.z * pFace->p2->dens) * pFace->norm.z;
      fRes = Const_One_Third * (fx + fy + fz);
      if(m_nClcVar == clcVolumeFlow)
        fRes /= Const_Air_Dens;
      break;
    }
    case clcEnergyFlow:
    {
      fRes = 0;
      break;
    }
    case clcAvePress:
    {
      fRes = Const_One_Third * (pFace->p0->press + pFace->p1->press + pFace->p2->press);
      break;
    }
    case clcAveTemp:
    {
      fRes = Const_One_Third * (pFace->p0->temp + pFace->p1->temp + pFace->p2->temp);
      break;
    }
    case clcAveDens:
    {
      fRes = Const_One_Third * (pFace->p0->dens + pFace->p1->dens + pFace->p2->dens);
      break;
    }
    case clcAveVx:
    {
      fRes = Const_One_Third * (pFace->p0->vel.x + pFace->p1->vel.x + pFace->p2->vel.x);
      break;
    }
    case clcAveVy:
    {
      fRes = Const_One_Third * (pFace->p0->vel.y + pFace->p1->vel.y + pFace->p2->vel.y);
      break;
    }
    case clcAveVz:
    {
      fRes = Const_One_Third * (pFace->p0->vel.z + pFace->p1->vel.z + pFace->p2->vel.z);
      break;
    }
    case clcAveRe:
    {
      double fFlowRate = Const_One_Third * 
        (pFace->p0->dens * pFace->p0->vel.length() + pFace->p1->dens * pFace->p1->vel.length() + pFace->p2->dens * pFace->p2->vel.length());

      double fDynVisc = Const_One_Third * (pFace->p0->visc + pFace->p1->visc + pFace->p2->visc);
      fRes = m_fCharLength * fFlowRate / fDynVisc;
    }
  }

  return true;
}

void CCalculator::normalize_result(double fSumSquare, double fSumRes)
{
  switch(m_nClcVar)
  {
    case clcMassFlow:
    case clcVolumeFlow:
    {
      m_fResult = fSumRes;
      break;
    }
    case clcAvePress:
    case clcAveTemp:
    case clcAveDens:
    case clcAveVx:
    case clcAveVy:
    case clcAveVz:
    case clcAveRe:
    {
      if(fSumSquare > Const_Almost_Zero)
        m_fResult = fSumRes / fSumSquare;

      break;
    }
  }
}

void CCalculator::convert_result()
{
  double fMult = 1.;
  switch(m_nClcVar)
  {
    case clcMassFlow:   fMult = 1e+5 * CGS_to_SI_Weight; break;
    case clcVolumeFlow: fMult = 0.06; break;  // convertion from cm^3/s to L/min.
    case clcEnergyFlow: fMult = CGS_to_SI_Energy; break;
    case clcAvePress:   fMult = CGS_to_SI_Press;  break;
    case clcAveDens:    fMult = CGS_to_SI_Dens;   break;
    case clcAveVx:
    case clcAveVy:
    case clcAveVz:
    {
      fMult = CGS_to_SI_Vel;
      break;
    }
  }

  m_fResult *= fMult;

// Symmetry correction (mass and energy flows must be multiplied by 2 or 4 dependent on the type of symmetry).
  if((m_pObj != NULL) && (m_nClcVar == clcMassFlow || m_nClcVar == clcEnergyFlow || m_nClcVar == clcVolumeFlow))
  {
    double fSymMult = 1;
    int nSymPlanes = m_pObj->get_sym_plane();
    if(nSymPlanes  & CTracker::spXY)
      fSymMult *= 2;
    if(nSymPlanes  & CTracker::spXZ)
      fSymMult *= 2;

    m_fResult *= fSymMult;
  }
}

const char* CCalculator::units() const
{
  switch(m_nClcVar)
  {
    case clcMassFlow:   return _T("Mass Flow,  10^-5 kg/s");
    case clcVolumeFlow: return _T("Volume Flow,  L/min");
    case clcEnergyFlow: return _T("Energy Flow, J/s");
    case clcAvePress:   return _T("Pressure,  Pa");
    case clcAveTemp:    return _T("Temperature,  K");
    case clcAveDens:    return _T("Density,  kg/m3");
    case clcAveRe:      return _T("Reynolds Number");
    case clcAveVx:      return _T("Vx,  m/s");
    case clcAveVy:      return _T("Vy,  m/s");
    case clcAveVz:      return _T("Vz,  m/s");
  }

  return _T("");
}

const char* CCalculator::get_var_name(int nVar) const
{
  switch(nVar)
  {
    case clcMassFlow:   return _T("Mass Flow");
    case clcVolumeFlow: return _T("Volume Flow (1 atm, 293 K)");
    case clcEnergyFlow: return _T("Energy Flow");
    case clcAvePress:   return _T("Average Pressure");
    case clcAveTemp:    return _T("Average Temperature");
    case clcAveDens:    return _T("Average Density");
    case clcAveRe:      return _T("Average Reynolds Number");
    case clcAveVx:      return _T("Average X-component of Velocity");
    case clcAveVy:      return _T("Average Y-component of Velocity");
    case clcAveVz:      return _T("Average Z-component of Velocity");
  }

  return _T("");
}

bool CCalculator::get_update_flag() const
{
  return true;
}

void CCalculator::save(CArchive& ar)
{
  UINT nVersion = 2;    // 2 - m_nCalcDir.
  ar << nVersion;

  ar << m_nClcVar;
  ar << m_bEnable;

  ar << m_fCharLength;
  ar << m_nCalcDir;
}

void CCalculator::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_nClcVar;
  ar >> m_bEnable;

  if(nVersion >= 1)
    ar >> m_fCharLength;

  if(nVersion >= 2)
    ar >> m_nCalcDir;
}

//-------------------------------------------------------------------------------------------------
// CPlaneYZCalculator.
//-------------------------------------------------------------------------------------------------
CPlaneYZCalculator::CPlaneYZCalculator()
  : CCalculator()
{
  set_default();
}

CPlaneYZCalculator::~CPlaneYZCalculator()
{
  m_pObj = NULL;
  clear();
}

void CPlaneYZCalculator::run()
{
  terminate(false);
  do_sequence_calc();
}

void CPlaneYZCalculator::set_default()
{
  m_sCalcName = calc_name(ctPlaneYZ);

// All positions and lengths are in cm.
  m_CrossSect.set_plane_type(CDomainCrossSection::ptPlaneYZ);
  m_CrossSect.set_plane_origin(Vector3D(0.1, 0.0, 0.0));

  m_bReady = false; // force to build the plane mesh.

// Sequential calculations support (since CCalculator nVersion 2 the sequence calculations can be done in all three directions, X, Y and Z):
  m_fStartX = 0.;
  m_fEndX = 1.;
  m_nSeqCalcCount = 100;
}

void CPlaneYZCalculator::clear()
{
  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->set_cross_sections_array(NULL);
}

void CPlaneYZCalculator::do_calculate()
{
  m_fResult = 0;
  if(!m_bReady && !build_cs_mesh())
    return;

  CRegion* pReg = m_CrossSect.get_region();
  double fS, fSumS = 0, fRes, fSumRes = 0;
  size_t nFaceCount = pReg->vFaces.size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    CFace* pFace = pReg->vFaces.at(i);
    if(!process_face(pFace, fS, fRes))
      continue;   // the face square is less than Const_Almost_Zero.

    fSumRes += fS * fRes;
    fSumS += fS;
  }

  normalize_result(fSumS, fSumRes);
  convert_result();

  if(m_bEnable)
  {
    CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
    pDrawObj->set_cross_sections_array(&(pReg->vFaces));
    pDrawObj->draw();
  }
}

void CPlaneYZCalculator::update()
{
  if(!m_bReady && !build_cs_mesh())
    return;

  CRegion* pReg = m_CrossSect.get_region();
  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(m_bEnable)
    pDrawObj->set_cross_sections_array(&(pReg->vFaces));
  else
    pDrawObj->set_cross_sections_array(NULL);
}

bool CPlaneYZCalculator::build_cs_mesh()
{
  clear();
  if(!get_tracker_ptr())
    return false;

  switch(m_nCalcDir)
  {
    case dirX: m_CrossSect.set_plane_type(CDomainCrossSection::ptPlaneYZ); break;
    case dirY: m_CrossSect.set_plane_type(CDomainCrossSection::ptPlaneXZ); break;
    case dirZ: m_CrossSect.set_plane_type(CDomainCrossSection::ptPlaneXY); break;
  }

  m_CrossSect.invalidate();
  CRegion* pReg = m_CrossSect.get_region();  // here CDomainCrossSection::build_mesh() is called.
  if(pReg == NULL)
    return false;

  size_t nFacesCount = pReg->vFaces.size();
  for(size_t i = 0; i < nFacesCount; i++)
    pReg->vFaces.at(i)->norm = m_CrossSect.get_plane_norm();

  m_bReady = true;
  return true;
}

double CPlaneYZCalculator::get_plane_x() const
{
  switch(m_nCalcDir)
  {
    case CCalculator::dirX: return m_CrossSect.get_plane_origin().x;
    case CCalculator::dirY: return m_CrossSect.get_plane_origin().y;
    case CCalculator::dirZ: return m_CrossSect.get_plane_origin().z;
  }

  return 0;
}

void CPlaneYZCalculator::set_plane_x(double fX)
{
  switch(m_nCalcDir)
  {
    case CCalculator::dirX: 
    {
      if(m_CrossSect.get_plane_origin().x != fX)
      {
        m_CrossSect.set_plane_origin(Vector3D(fX, 0, 0));
        m_bReady = false;
        break;
      }
    }
    case CCalculator::dirY: 
    {
      if(m_CrossSect.get_plane_origin().y != fX)
      {
        m_CrossSect.set_plane_origin(Vector3D(0, fX, 0));
        m_bReady = false;
        break;
      }
    }
    case CCalculator::dirZ: 
    {
      if(m_CrossSect.get_plane_origin().z != fX)
      {
        m_CrossSect.set_plane_origin(Vector3D(0, 0, fX));
        m_bReady = false;
        break;
      }
    }
  }
}

//---------------------------------------------------------------------------------------
// Sequential calculations support:
//---------------------------------------------------------------------------------------
void CPlaneYZCalculator::do_sequence_calc()
{
  if(m_nSeqCalcCount == 0)
    return;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, m_sOutputFile.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  std::string cAction(_T("Calculating "));
  cAction += std::string(get_var_name(m_nClcVar));
  set_job_name(cAction.c_str());

  CString sCoordName = m_nCalcDir == dirX ? CString(_T("  x(mm),    ")) : (m_nCalcDir == dirY ? CString(_T("  y(mm),    ")) : CString(_T("  z(mm),    ")));
  std::string sHeader((const char*)sCoordName);
  sHeader += units();
  fputs(sHeader.c_str(), pStream);
  fputs("\n\n", pStream);

  double fCurrX = get_plane_x();
  double fStepX = (m_fEndX - m_fStartX) / m_nSeqCalcCount;

  for(UINT i = 0; i < m_nSeqCalcCount; i++)
  {
    set_progress(100 * i / m_nSeqCalcCount);
    if(get_terminate_flag())
      break;

    fCurrX = m_fStartX + i * fStepX;
    set_plane_x(fCurrX);
    do_calculate();

    fprintf(pStream, "%f,  %f\n", 10 * fCurrX, m_fResult);
  }

  fclose(pStream);
  if(get_terminate_flag())
    return;

  set_progress(100);
  set_plane_x(fCurrX);
  do_calculate();
}

//---------------------------------------------------------------------------------------
// Streamability:
//---------------------------------------------------------------------------------------
void CPlaneYZCalculator::save(CArchive& ar)
{
  UINT nVersion = 2;
  ar << nVersion;

  CCalculator::save(ar);

  double fCurrX = get_plane_x();
  ar << fCurrX;

  ar << m_fStartX;
  ar << m_fEndX;
  ar << m_nSeqCalcCount;

  CString cFileName(m_sOutputFile.c_str());
  ar << cFileName;
}

void CPlaneYZCalculator::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CCalculator::load(ar);

  double fCurrX;
  ar >> fCurrX;
  set_plane_x(fCurrX);
  if(nVersion < 2)
  {
    double fMeshSize;
    ar >> fMeshSize;
  }

  if(nVersion >= 1)
  {
    ar >> m_fStartX;
    ar >> m_fEndX;
    ar >> m_nSeqCalcCount;

    CString cFileName;
    ar >> cFileName;
    set_filename((const char*)cFileName);
  }
}

//-------------------------------------------------------------------------------------------------
// CSelectedRegionCalculator.
//-------------------------------------------------------------------------------------------------
CSelectedRegionCalculator::CSelectedRegionCalculator()
  : CCalculator()
{
  m_sCalcName = calc_name(ctSelRegions);
}

CSelectedRegionCalculator::~CSelectedRegionCalculator()
{
  m_pObj = NULL;
  clear();
}

void CSelectedRegionCalculator::run()
{
  terminate(false);
  do_calculate();
}

void CSelectedRegionCalculator::do_calculate()
{
  m_fResult = 0;
  if(!get_tracker_ptr())
    return;

  UINT nTotalFaceCount;
  CRegionsCollection vRegions;
  if(!get_sel_regions(vRegions, nTotalFaceCount))
    return;

  std::string cAction(_T("Calculation "));
  cAction += std::string(get_var_name(m_nClcVar));
  set_job_name(cAction.c_str());

  double fS, fSumS = 0, fRes, fSumRes = 0;
  size_t nRegCount = vRegions.size(), nFaceDone = 0;
  for(size_t j = 0; j < nRegCount; j++)
  {
    CRegion* pReg = vRegions.at(j);
    size_t nFaceCount = pReg->vFaces.size();
    for(size_t i = 0; i < nFaceCount; i++)
    {
      if(get_terminate_flag())
        return;

      if(nFaceDone % 5 == 0)
        set_progress(100 * nFaceDone / nTotalFaceCount);

      nFaceDone++;

      CFace* pFace = pReg->vFaces.at(i);
      if(!process_face(pFace, fS, fRes))
        continue;   // the face square is less than Const_Almost_Zero.

      fSumRes += fS * fRes;
      fSumS += fS;
    }
  }

  normalize_result(fSumS, fSumRes);
  convert_result();
}

void CSelectedRegionCalculator::update()
{
}

void CSelectedRegionCalculator::clear()
{
}

bool CSelectedRegionCalculator::get_sel_regions(CRegionsCollection& vRegions, UINT& nFaceCount) const
{
  const CRegionsCollection& vAllRegions = m_pObj->get_regions();
  size_t nRegCount = vAllRegions.size();
  nFaceCount = 0;
  for(size_t i = 0; i < nRegCount; i++)
  {
    CRegion* pReg = vAllRegions.at(i);
    if(std::find(m_vSelRegNames.begin(), m_vSelRegNames.end(), pReg->sName) != m_vSelRegNames.end())
    {
      vRegions.push_back(pReg);
      nFaceCount += pReg->vFaces.size();
    }
  }

  return vRegions.size() != 0;
}

void CSelectedRegionCalculator::save(CArchive& ar)
{
  UINT nVersion = 1;
  ar << nVersion;

  CCalculator::save(ar);

  size_t nRegCount = m_vSelRegNames.size();
  ar << nRegCount;
  for(size_t i = 0; i < nRegCount; i++)
  {
    CString cRegName(m_vSelRegNames.at(i).c_str());
    ar << cRegName;
  }
}

void CSelectedRegionCalculator::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CCalculator::load(ar);

  size_t nRegCount;
  ar >> nRegCount;
  for(size_t i = 0; i < nRegCount; i++)
  {
    CString cRegName;
    ar >> cRegName;
    m_vSelRegNames.push_back((const char*)cRegName);
  }
}

//-------------------------------------------------------------------------------------------------
// CLineCalculator.
//-------------------------------------------------------------------------------------------------
CLineCalculator::CLineCalculator()
  : CCalculator()
{
  m_sCalcName = calc_name(ctAlongLine);
  set_default();
}

CLineCalculator::~CLineCalculator()
{
  m_pObj = NULL;
  clear();
}

void CLineCalculator::run()
{
  terminate(false);
  do_calculate();
}

void CLineCalculator::set_default()
{
  m_nStepCount = 50;
  m_vLineStart = Vector3D(0, 0.001, 0.001);
  m_vLineEnd = Vector3D(5, 0.001, 0.001);
  m_nClcVar = lcTemp;

  m_fLineAverage = 0;
}

void CLineCalculator::do_calculate()
{
  if(!get_tracker_ptr())
    return;

  Vector3D vDir = m_vLineEnd - m_vLineStart;
  double fLineLen = vDir.length();
  if(m_nStepCount < 1 || fLineLen < Const_Almost_Zero)
    return;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, m_sOutputFile.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  std::string cAction(_T("Calculation "));
  cAction += std::string(get_var_name(m_nClcVar));
  set_job_name(cAction.c_str());

  std::string sHeader(_T("  x(mm),  y(mm),  z(mm),     "));
  sHeader += units();
  fputs(sHeader.c_str(), pStream);
  fputs("\n\n", pStream);

  Vector3D vPos;
  vDir.normalize();
  double fStep = fLineLen / m_nStepCount;
  double fPhase = m_pObj->get_enable_ansys_field() ? 0.0 : Const_Half_PI; // for RF-fields presentation only.
  const CElem3D* pElem = NULL;
  CNode3D node;
  m_fLineAverage = 0;
  for(UINT i = 0; i <= m_nStepCount; i++)
  {
    set_progress(100 * i / m_nStepCount);
    if(get_terminate_flag())
      break;

    vPos = m_vLineStart + vDir * (i * fStep);
    if(!m_pObj->interpolate(vPos, 0, fPhase, node, pElem))
      continue;

    assign_result(node);

    if(m_nClcVar == lcDens)
      fprintf(pStream, "%f, %f, %f, %e\n", 10 * vPos.x, 10 * vPos.y, 10 * vPos.z, m_fResult);
    else
      fprintf(pStream, "%f, %f, %f, %f\n", 10 * vPos.x, 10 * vPos.y, 10 * vPos.z, m_fResult);

    m_fLineAverage += m_fResult;
  }

  m_fLineAverage /= (m_nStepCount + 1);
  fclose(pStream);
}

void CLineCalculator::assign_result(const CNode3D& node)
{
  switch(m_nClcVar)
  {
    case lcPress: m_fResult = CGS_to_SI_Press * node.press; break;
    case lcTemp:  m_fResult = node.temp; break;
    case lcDens:  m_fResult = CGS_to_SI_Dens * node.dens; break;
    case lcVx:    m_fResult = CGS_to_SI_Vel * node.vel.x; break;
    case lcRe:    m_fResult = node.dens * node.vel.length() * m_fCharLength / node.visc; break;
    case lcEx:    m_fResult = CGS_to_SI_Voltage * (node.field.x + m_pObj->get_field_ptb().apply(node.pos).x); break;
    case lcEy:    m_fResult = CGS_to_SI_Voltage * (node.field.y + m_pObj->get_field_ptb().apply(node.pos).y); break;
    case lcEz:    m_fResult = CGS_to_SI_Voltage * (node.field.z + m_pObj->get_field_ptb().apply(node.pos).z); break;
    case lcPhiDC: m_fResult = node.phi; break;
    case lcRFx:
    {
      m_fResult = m_pObj->get_enable_ansys_field() ? m_pObj->get_rf_field(node, 0, Const_Half_PI).x : node.rf.x;
      m_fResult *= CGS_to_SI_Voltage;
      break;
    }
    case lcRFy:
    {
      m_fResult = m_pObj->get_enable_ansys_field() ? m_pObj->get_rf_field(node, 0, Const_Half_PI).y : node.rf.y;
      m_fResult *= CGS_to_SI_Voltage;
      break;
    }
    case lcRFz:
    {  
      m_fResult = m_pObj->get_enable_ansys_field() ? m_pObj->get_rf_field(node, 0, Const_Half_PI).z : node.rf.z;
      m_fResult *= CGS_to_SI_Voltage;
      break;
    }
  }
}

void CLineCalculator::update()
{
}

void CLineCalculator::clear()
{
}

const char* CLineCalculator::units() const
{
  switch(m_nClcVar)
  {
    case lcPress:   return _T("Pressure,  Pa");
    case lcTemp:    return _T("Temperature,  K");
    case lcDens:    return _T("Density,  kg/m3");
    case lcVx:      return _T("Vx,  m/s");
    case lcRe:      return _T("Reynolds Number");
    case lcEx:      return _T("DC Ex,  V/cm");
    case lcEy:      return _T("DC Ey,  V/cm");
    case lcEz:      return _T("DC Ez,  V/cm");
    case lcRFx:     return _T("RF Ex,  V/cm");
    case lcRFy:     return _T("RF Ey,  V/cm");
    case lcRFz:     return _T("RF Ez,  V/cm");
    case lcPhiDC:   return _T("Electric Potential, V");
  }
  return _T("");
}

const char* CLineCalculator::get_var_name(int nVar) const
{
  switch(nVar)
  {
    case lcPress:   return _T("Pressure");
    case lcTemp:    return _T("Temperature");
    case lcDens:    return _T("Density");
    case lcVx:      return _T("Velocity X");
    case lcRe:      return _T("Reynolds Number");
    case lcEx:      return _T("DC Field X");
    case lcEy:      return _T("DC Field Y");
    case lcEz:      return _T("DC Field Z");
    case lcRFx:     return _T("RF Field X");
    case lcRFy:     return _T("RF Field Y");
    case lcRFz:     return _T("RF Field Z");
    case lcPhiDC:   return _T("Electric Potential");
  }

  return _T("");
}

bool CLineCalculator::get_update_flag() const
{
  return false;
}

void CLineCalculator::save(CArchive& ar)
{
  UINT nVersion = 1;
  ar << nVersion;

  CCalculator::save(ar);

  CString cName(m_sOutputFile.c_str());
  ar << cName;
  
  ar << m_vLineStart.x;
  ar << m_vLineStart.y;
  ar << m_vLineStart.z;

  ar << m_vLineEnd.x;
  ar << m_vLineEnd.y;
  ar << m_vLineEnd.z;

  ar << m_nStepCount;
}

void CLineCalculator::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CCalculator::load(ar);

  CString cName;
  ar >> cName;
  m_sOutputFile = std::string((const char*)cName);

  ar >> m_vLineStart.x;
  ar >> m_vLineStart.y;
  ar >> m_vLineStart.z;

  ar >> m_vLineEnd.x;
  ar >> m_vLineEnd.y;
  ar >> m_vLineEnd.z;

  ar >> m_nStepCount;
}

//-------------------------------------------------------------------------------------------------
// CTrackCalculator - a class for calculating averaged ion parameters by track data.
//-------------------------------------------------------------------------------------------------
CTrackCalculator::CTrackCalculator()
  : CCalculator()
{
  m_sCalcName = calc_name(ctTrackCalc);
  set_default();
}

CTrackCalculator::~CTrackCalculator()
{
  m_pObj = NULL;
  clear();
}

void CTrackCalculator::run()
{
  terminate(false);
  do_sequence_calc();
}

void CTrackCalculator::set_default()
{
  m_bReady = false; // post-process fragmentation computation support.
  m_nClcVar = clcIonTemp;
  m_nCrossSectCount = 50;
  m_fPosCS = 0;
  m_fStartPos = 0;
  m_fEndPos = 1;

// Run-time:
  m_fIonTempEq = 0;
  m_fGasTemp = 0;

  m_fMaxDropDiam = 0;
  m_fMinDropTemp = 0;
  m_fMaxDropTemp = 0;
  m_fRayleigh = 0;

  m_fPartEvapor = 0;
  m_fPartRayleigh = 0;
  m_fPartHitWall = 0;
}

bool CTrackCalculator::valid_items(const CTrackItem& prev, const CTrackItem& curr, double& fKsi) const
{
  bool bValid = false;
  switch(m_nCalcDir)
  {
    case CCalculator::dirX:
    {
      bValid = (curr.pos.x >= m_fPosCS) && (prev.pos.x < m_fPosCS) || (curr.pos.x < m_fPosCS) && (prev.pos.x >= m_fPosCS);
      if(bValid)
      {
        double fdX = curr.pos.x - prev.pos.x;
        if(fabs(fdX) < Const_Almost_Zero)
          return false;

        fKsi = (m_fPosCS - prev.pos.x) / fdX;
      }
      break;
    }
    case CCalculator::dirY:
    {
      bValid = (curr.pos.y >= m_fPosCS) && (prev.pos.y < m_fPosCS) || (curr.pos.y < m_fPosCS) && (prev.pos.y >= m_fPosCS);
      if(bValid)
      {
        double fdY = curr.pos.y - prev.pos.y;
        if(fabs(fdY) < Const_Almost_Zero)
          return false;

        fKsi = (m_fPosCS - prev.pos.y) / fdY;
      }
      break;
    }
    case CCalculator::dirZ:
    {
      bValid = (curr.pos.z >= m_fPosCS) && (prev.pos.z < m_fPosCS) || (curr.pos.z < m_fPosCS) && (prev.pos.z >= m_fPosCS);
      if(bValid)
      {
        double fdZ = curr.pos.z - prev.pos.z;
        if(fabs(fdZ) < Const_Almost_Zero)
          return false;

        fKsi = (m_fPosCS - prev.pos.z) / fdZ;
      }
      break;
    }
  }

  return bValid;
}

void CTrackCalculator::do_calculate()
{
  if(!get_tracker_ptr())
    return;

  CTrackVector& vTracks = m_pObj->get_tracks();
  if(m_nClcVar == clcIonTemp)
    collect_elements();
  if((m_nClcVar == clcFragment) && !calc_fragmentation())
    return;

  m_fResult = m_fIonTempEq = m_fGasTemp = m_fMaxDropDiam = m_fMaxDropTemp = m_fRayleigh = 0;
  m_fMinDropTemp = FLT_MAX;
  m_nIntersectCount = 0;

  Vector3D vPos, vIonVel;
  size_t nTrackCount = vTracks.size();
  if(nTrackCount < 1)
    return;

  if(m_nClcVar == clcTerminated)
  {
    m_fResult = calc_term_tracks(); // part of terminated tracks.
    return;
  }

  CTrackItem curr, prev;
  double fKsi, fdX, fIonT, fGasT, fDropMass, fDropD, fDropTemp, fCrit;
  double fCurrIncr = m_pObj->get_full_current() / nTrackCount;
  for(size_t i = 0; i < nTrackCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    bool bCalcTime = track.get_type() == CTrack::ptDroplet ? true : (m_vReachedLastCS.size() == nTrackCount ? m_vReachedLastCS.at(i) : true);
    size_t nItemCount = track.size();

    for(size_t j = 1; j < nItemCount; j++)
    {
      track.get_track_item(j, curr);
      track.get_track_item(j - 1, prev);

      if(valid_items(prev, curr, fKsi))
      {
        m_nIntersectCount++;

        switch(m_nClcVar)
        {
// Ion parameters:
          case clcIonTemp:
          {
            if(track.get_type() != CTrack::ptIon)
              break;

            vPos = prev.pos + fKsi * (curr.pos - prev.pos);
            if(calc_gas_temp(vPos, fGasT))
            {
              m_fResult += prev.temp + fKsi * (curr.temp - prev.temp);
              m_fIonTempEq += prev.tempinf + fKsi * (curr.tempinf - prev.tempinf);
              m_fGasTemp += fGasT;
            }

            break;
          }
          case clcCurrent:
          {
            if(track.get_type() != CTrack::ptIon)
              break;

            m_fResult += fCurrIncr;
            break;
          }
          case clcFragment:
          {
            if(track.get_type() != CTrack::ptIon)
              break;

            m_fResult += prev.unfragm + fKsi * (curr.unfragm - prev.unfragm);
            break;
          }
// Droplet parameters:
          case clcDropDiameter:
          {
            if(track.get_type() != CTrack::ptDroplet)
              break;

            fDropMass = prev.mass + fKsi * (curr.mass - prev.mass);
            fDropD = m_pObj->get_particle_diameter(fDropMass);
            if(m_fMaxDropDiam < fDropD)
              m_fMaxDropDiam = fDropD;

            m_fResult += fDropD;

            fDropTemp = prev.temp + fKsi * (curr.temp - prev.temp);
            fCrit = m_pObj->get_particle_charge() / m_pObj->get_max_charge(fDropTemp, fDropD);
            m_fRayleigh += fCrit;
            break;
          }
          case clcDropTemp:
          {
            if(track.get_type() != CTrack::ptDroplet)
              break;

            fDropTemp = prev.temp + fKsi * (curr.temp - prev.temp);
            m_fResult += fDropTemp;

            if(m_fMinDropTemp > fDropTemp)
              m_fMinDropTemp = fDropTemp;
            if(m_fMaxDropTemp < fDropTemp)
              m_fMaxDropTemp = fDropTemp;

            break;
          }
          case clcTime:
          {
// If this function is called from do_sequence_calc() the vector m_vReachedLastCS contains "true" for those
// tracks which have reached m_fEndPos. It is only those tracks using which the time must be averaged.
            if(bCalcTime)
            {
              m_fResult += prev.time + fKsi * (curr.time - prev.time);
              break;
            }
          }
        }

        break;
      }
    }
  }

  convert_result();
  if(m_nClcVar == clcCurrent)
    return;

  if(m_nIntersectCount != 0)
  {
    double fNormCoeff = 1. / m_nIntersectCount;
    m_fResult *= fNormCoeff;
    if(m_nClcVar == clcFragment)
      m_fResult = 100 * (1.- m_fResult);

    m_fIonTempEq *= fNormCoeff;
    m_fGasTemp *= fNormCoeff;
    m_fRayleigh *= fNormCoeff;
  }
}

void CTrackCalculator::do_sequence_calc()
{
  if(!m_bEnable || !get_tracker_ptr())
    return;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, m_sOutputFile.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  std::string cAction(_T("Calculating "));
  cAction += std::string(get_var_name(m_nClcVar));
  set_job_name(cAction.c_str());

  std::string sHeader;
  std::string sUnits = m_nCalcDir == dirX ? _T("  x(mm),") : (m_nCalcDir == dirY ? _T("  y(mm),") : _T("  z(mm),"));
  switch(m_nClcVar)
  {
    case clcIonTemp:  sHeader = _T("      TI(K),      TIeq(K),    Tgas(K)"); break;
    case clcCurrent:  sHeader = _T("   current(nA),   trans(percent)"); break;
    case clcFragment: sHeader = _T("   fragmented(percent)"); break;
    case clcDropDiameter: sHeader = _T("  DropD(mcm),  MaxDropD(mcm), Raylegh(%)"); break;
    case clcDropTemp: sHeader = _T(" AverDropletTemp(K), Tmin(K),  Tmax(K)"); break;
    case clcTerminated: sHeader = _T(" Term(%), OvrTime(%), Evapor(%), Rayleigh(%), HitWall(%)"); break;
    case clcTime:     sHeader = _T("   time(ms)"); break;
  }

  fputs((sUnits + sHeader).c_str(), pStream);
  fputs("\n\n", pStream);

  if(m_nClcVar == clcTime)
    find_reached_last_cs();

  double fCurrPos = m_fPosCS;
  double fStep = m_nCrossSectCount > 1 ? (m_fEndPos - m_fStartPos) / (m_nCrossSectCount - 1) : 0;
  double fFullCurr = m_pObj->get_full_current() / Const_nA_to_CGSE;   // initial current in nA.
  for(UINT i = 0; i < m_nCrossSectCount; i++)
  {
    set_progress(100 * i / m_nCrossSectCount);
    if(get_terminate_flag())
      break;

    m_fPosCS = m_fStartPos + i * fStep;
    do_calculate();
    if(m_nIntersectCount == 0)
      continue;

    switch(m_nClcVar)
    {
      case clcIonTemp:  fprintf(pStream, "%f,  %f,  %f,  %f\n", 10 * m_fPosCS, m_fResult, m_fIonTempEq, m_fGasTemp); break;
      case clcCurrent:  fprintf(pStream, "%f,  %f,  %f\n", 10 * m_fPosCS, m_fResult, 100 * m_fResult / fFullCurr); break;
      case clcFragment: fprintf(pStream, "%f,  %f\n", 10 * m_fPosCS, m_fResult); break;
      case clcDropDiameter: fprintf(pStream, "%f,  %f,  %f,  %f\n", 10 * m_fPosCS, m_fResult, m_fMaxDropDiam, 100 * m_fRayleigh); break;
      case clcDropTemp: fprintf(pStream, "%f,  %f,  %f,  %f\n", 10 * m_fPosCS, m_fResult, m_fMinDropTemp, m_fMaxDropTemp); break;
      case clcTerminated: fprintf(pStream, "%f,  %f,  %f,  %f,  %f,  %f\n", 10 * m_fPosCS, m_fResult, m_fPartOvertime, m_fPartEvapor, m_fPartRayleigh, m_fPartHitWall); break;
      case clcTime:     fprintf(pStream, "%f,  %f\n", 10 * m_fPosCS, m_fResult); break;
    }
  }

  fclose(pStream);
  if(get_terminate_flag())
  {
    m_vReachedLastCS.clear();
    return;
  }

  set_progress(100);
  m_fPosCS = fCurrPos;
  m_vReachedLastCS.clear();
  do_calculate();
}

void CTrackCalculator::find_reached_last_cs()
{
  m_vReachedLastCS.clear();
  CTrackVector& vTracks = m_pObj->get_tracks();
  size_t nTrackCount = vTracks.size();
  m_vReachedLastCS.assign(nTrackCount, false);

  CTrackItem item;
  for(size_t i = 0; i < nTrackCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    size_t nItemCount = track.size();
    for(size_t j = 1; j < nItemCount; j++)
    {
      track.get_track_item(j, item);
      if(item.pos.x >= m_fEndPos)
      {
        m_vReachedLastCS[j] = true;
        break;
      }
    }
  }
}

void CTrackCalculator::update()
{
}

void CTrackCalculator::clear()
{
  m_vElements.clear();
}

const char* CTrackCalculator::units() const
{
  switch(m_nClcVar)
  {
    case clcIonTemp:  return _T("Aver. Ion Temperature,  K");
    case clcCurrent:  return _T("Ion Current,  nA");
    case clcFragment: return _T("Part of fragmented Ions, %");
    case clcDropDiameter: return _T("Aver. Droplet Diameter, mcm");
    case clcDropTemp: return _T("Aver. Droplet Temperature, K");
    case clcTerminated: return _T("Part of Terminated Tracks, %");
    case clcTime:     return _T("Moving Time, ms");
  }

  return _T("");
}

const char* CTrackCalculator::get_var_name(int nVar) const
{
  switch(nVar)
  {
    case clcIonTemp:  return _T("Average Ion Temperature");
    case clcCurrent:  return _T("Ion Current");
    case clcFragment: return _T("Part of fragmented Ions");
    case clcDropDiameter: return _T("Average Droplet Diameter");
    case clcDropTemp: return _T("Average Droplet Temperature");
    case clcTerminated: return _T("Part of Terminated Tracks");
    case clcTime:     return _T("Moving Time");
  }

  return _T("");
}

void CTrackCalculator::convert_result()
{
  switch(m_nClcVar)
  {
    case clcCurrent:  m_fResult /= Const_nA_to_CGSE; break;   // nA.
    case clcDropDiameter: m_fResult *= 10000; m_fMaxDropDiam *= 10000; break; // droplet diameter in mcm;
    case clcTime:     m_fResult *= 1000; break;               // milliseconds.
  }
}

static const double scfTempCoeff = 1. / Const_Boltzmann;

double CTrackCalculator::ion_temp(double fGasTemp, const Vector3D& vGasVel, const Vector3D& vIonVel)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  double fGasMass = pObj->get_molar_mass();   // mass of the air molecule, g.
  return fGasTemp + 0.2 * scfTempCoeff * fGasMass * (vIonVel - vGasVel).sqlength();
}

// Calculation of the stationary ion temperature from scratch:
bool CTrackCalculator::collect_elements()
{
  m_vElements.clear();
  const CElementsCollection& all_elems = m_pObj->get_elems();
  size_t nElemCount = all_elems.size();
  for(size_t i = 0; i < nElemCount; i++)
  {
    CElem3D* pElem = all_elems.at(i);
    if(m_fPosCS < pElem->box.vMin.x || m_fPosCS > pElem->box.vMax.x)
      continue;

    m_vElements.push_back(pElem);
  }

  return m_vElements.size() > 0;
}

bool CTrackCalculator::calc_ion_temp(const Vector3D& vIonPos, const Vector3D& vIonVel, double& fIonTemp, double& fGasTemp) const
{
  Vector3D vReflPos = vIonPos;
  Vector3D vReflVel = vIonVel;
  m_pObj->sym_corr_forward(vReflPos, vReflVel); // axial symmetry support.

  size_t nElemCount = m_vElements.size();
  for(size_t i = 0; i < nElemCount; i++)
  {
    CElem3D* pElem = m_vElements.at(i);
    if(!pElem->inside(vReflPos))
      continue;

    CNode3D node;
    pElem->interpolate(vReflPos, node);
    fIonTemp = ion_temp(node.temp, node.vel, vReflVel);
    fGasTemp = node.temp;
    return true;
  }

  return false;
}

bool CTrackCalculator::calc_gas_temp(const Vector3D& vIonPos, double& fGasTemp) const
{
  Vector3D vVel;
  Vector3D vReflPos = vIonPos;
  m_pObj->sym_corr_forward(vReflPos, vVel); 

  size_t nElemCount = m_vElements.size();
  for(size_t i = 0; i < nElemCount; i++)
  {
    CElem3D* pElem = m_vElements.at(i);
    if(!pElem->inside(vReflPos))
      continue;

    CNode3D node;
    pElem->interpolate(vReflPos, node);
    fGasTemp = node.temp;
    return true;
  }

  return false;
}

// Fragmentation in the post-process:
bool CTrackCalculator::calc_fragmentation()
{
  if(m_bReady)
    return true;

  if(m_pObj->get_particle_type() != CTrack::ptIon)
    return false;

  CIonTrackItem* pItem = NULL;
  double fPrevTime, fPrevProb, fCurrProb;
  CTrackVector& vTracks = m_pObj->get_tracks();
  size_t nTrackCount = vTracks.size();
  for(size_t i = 0; i < nTrackCount; i++)
  {
    CTrack& track = vTracks.at(i);
    size_t nItemCount = track.size();
    if(nItemCount < 2 || track.get_type() != CTrack::ptIon)
      continue;

    double fInt = 0;
    pItem = (CIonTrackItem*)track.at(0);
    fPrevTime = pItem->time;
    if(!get_fragm_probability(pItem, fPrevProb))
      continue;

    for(size_t j = 1; j < nItemCount; j++)
    {
      pItem = (CIonTrackItem*)track.at(j);
      if(!get_fragm_probability(pItem, fCurrProb))
        continue;

      fInt += 0.5 * (fPrevProb + fCurrProb) * (pItem->time - fPrevTime);
      pItem->unfragm = exp(-fInt);

      fPrevTime = pItem->time;
      fPrevProb = fCurrProb;
    }
  }

  m_bReady = true;
  return true;
}

bool CTrackCalculator::get_fragm_probability(CIonTrackItem* pItem, double& fFragmProb)
{
  CElem3D* pElem = pItem->nElemId >= 0 ? m_pObj->get_elems().at(pItem->nElemId) : NULL;
  if(pElem == NULL)
    return false;

  Vector3D vIonPos = pItem->pos, vIonVel = pItem->vel;
  m_pObj->sym_corr_forward(vIonPos, vIonVel);

  double pWeight[6];
  if(!pElem->coeff(vIonPos, pWeight))
    return false;

// Environment gas parameters:
  Vector3D vGasVel(0, 0, 0);
  float w, fPress = 0, fTemp = 0;
  const CNodesVector& vNodes = m_pObj->get_nodes();
  size_t nNodeCount = pElem->get_node_count();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    const CNode3D& node = vNodes.at(pElem->get_node_index(i));
    w = pWeight[i];

    vGasVel += w * node.vel;
    fPress += w * node.press;
    fTemp += w * node.temp;
  }

  double fNumDens = fPress * scfTempCoeff / fTemp;
  double fIonNeutrFreq = fNumDens * (vIonVel - vGasVel).length() * m_pObj->get_ion_cross_section();

  double fActEnergy = m_pObj->get_act_energy() / Const_Erg_to_EV; // activation energy in CGS.
  double fProbInt = CMath::energy_probability_integral(fActEnergy * scfTempCoeff / pItem->temp);

  fFragmProb = fIonNeutrFreq * fProbInt;

  return true;
}

double CTrackCalculator::calc_term_tracks()
{
  CTrackItem item;
  CTrackVector& vTracks = m_pObj->get_tracks();
  size_t nTrackCount = vTracks.size();
  if(nTrackCount == 0)
    return 0;

  double fT, fD, fCrit;
  UINT nEvapor = 0, nRayleigh = 0, nHitWall = 0, nOvrTime = 0;

  for(size_t i = 0; i < nTrackCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    bool bIsDroplet = track.get_type() == CTrack::ptDroplet;

    size_t nItemCount = track.size();
    if(nItemCount < 1)
      continue;

    size_t nLast = nItemCount - 1;
    track.get_track_item(nLast, item);
    bool bContinue =
      (m_nCalcDir == dirX) && (item.pos.x > m_fPosCS) || (m_nCalcDir == dirY) && (item.pos.y > m_fPosCS) || (m_nCalcDir == dirZ) && (item.pos.z > m_fPosCS);
    if(bContinue)
      continue;   // we are not interested in tracks that have reached the current cross-section.

    int nTermReason = track.get_term_reason();
    switch(nTermReason)
    {
      case CTrack::ttrOvrTime: nOvrTime++; break;
      case CTrack::ttrFullyEvapor: nEvapor++; break;
      case CTrack::ttrRayleighLim: nRayleigh++; break;
      case CTrack::ttrNone: nHitWall++; break;
    }
  }

  double fNrmCoeff = 1. / nTrackCount;
  m_fPartOvertime = 100 * nOvrTime * fNrmCoeff;
  m_fPartEvapor = 100 * nEvapor * fNrmCoeff;
  m_fPartRayleigh = 100 * nRayleigh * fNrmCoeff;
  m_fPartHitWall = 100 * nHitWall * fNrmCoeff;
  m_nIntersectCount = 1;  // arbitrary non-zero number, this value is not used for terminated tracks calculation.
  return m_fPartOvertime + m_fPartEvapor + m_fPartRayleigh + m_fPartHitWall;
}

void CTrackCalculator::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CCalculator::save(ar);

  CString cName(m_sOutputFile.c_str());
  ar << cName;
  
  ar << m_fPosCS;
  ar << m_fStartPos;
  ar << m_fEndPos;

  ar << m_nCrossSectCount;
}

void CTrackCalculator::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CCalculator::load(ar);

  CString cName;
  ar >> cName;
  m_sOutputFile = std::string((const char*)cName);

  ar >> m_fPosCS;
  ar >> m_fStartPos;
  ar >> m_fEndPos;

  ar >> m_nCrossSectCount;
}

//-------------------------------------------------------------------------------------------------
// CTrackCrossSectionCalculator - a class for output ion parameters distributions in a given cross-section.
//-------------------------------------------------------------------------------------------------
CTrackCrossSectionCalculator::CTrackCrossSectionCalculator()
  : CCalculator()
{
  m_sCalcName = calc_name(ctTrackCrossSect);
}

CTrackCrossSectionCalculator::~CTrackCrossSectionCalculator()
{
}

void CTrackCrossSectionCalculator::run()
{
  terminate(false);
  do_calculate();
}

void CTrackCrossSectionCalculator::default_filename()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  std::string cPath = COutputEngine::get_full_path(pObj->get_filename());
  std::string cName("cross_sect_");

  double fCrssSctX = m_Object.get_cross_sect_pos();
  double fHalf = fCrssSctX >= 0 ? 0.5 : -0.5;
  fCrssSctX = 0.001 * int(fHalf + 10000 * fCrssSctX);  // cross-section X in mm with accuracy 2 digits after decimal point.
  char buff[_CVTBUFSIZE];
  if(_gcvt(fCrssSctX, 5, buff) == NULL) // convertation double into array of char.
    return;

  std::string cX(buff);
  std::string cExt("_mm.dat");
  m_sOutputFile = cPath + cName + cX + cExt;
}

void CTrackCrossSectionCalculator::do_calculate()
{
  if(m_sOutputFile.empty())
    default_filename();

  m_Object.set_norm_dir(m_nCalcDir);

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, m_sOutputFile.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  m_Object.set_handlers(m_hJobNameHandle, m_hProgressBarHandle);

  CString sCoord = m_nCalcDir == dirX ? CString(" y(mm),  z(mm),  ") : (m_nCalcDir == dirY ? CString(" x(mm),  z(mm),  ") : CString(" x(mm),  y(mm),  "));
  CString sVar;
  int nVar = m_Object.get_var();
  switch(nVar)
  {
    case CColoredCrossSection::varIonTemp:  sVar = CString("Color_Ti,   Ti(K)"); break;
    case CColoredCrossSection::varStartRadius: sVar = CString("Color_R0,  R0(cm)"); break;
  }

  fprintf(pStream, "%s\n\n", (const char*)(sCoord + sVar));

  m_Object.draw();

  calculate_statistics();

  double fVal;
  Vector3D vPoint;
  COLORREF nColor;
  RGB_Color clr;
  size_t nCount = m_Object.get_points_count();
  for(size_t i = 0; i < nCount; i++)
  {
    m_Object.get_colored_point(i, vPoint, fVal, clr);

    switch(m_nCalcDir)
    {
      case dirX: fprintf(pStream, "%f, %f, %d, %f\n", 10 * vPoint.y, 10 * vPoint.z, RGB(clr.red, clr.green, clr.blue), fVal); break;
      case dirY: fprintf(pStream, "%f, %f, %d, %f\n", 10 * vPoint.x, 10 * vPoint.z, RGB(clr.red, clr.green, clr.blue), fVal); break;
      case dirZ: fprintf(pStream, "%f, %f, %d, %f\n", 10 * vPoint.x, 10 * vPoint.y, RGB(clr.red, clr.green, clr.blue), fVal); break;
    }
  }

  fclose(pStream);
  m_Object.set_handlers(NULL, NULL);
}

static const Vector3D vNull(0, 0, 0);

void CTrackCrossSectionCalculator::calculate_statistics()
{
  m_vCenter = vNull;
  m_vSigma = vNull;

  double fVal;
  RGB_Color clr;
  Vector3D vPoint;
  size_t nCount = m_Object.get_points_count();
  if(nCount == 0)
    return;

  for(size_t i = 0; i < nCount; i++)
  {
    m_Object.get_colored_point(i, vPoint, fVal, clr);
    m_vCenter += vPoint;
  }

  m_vCenter /= (double)nCount;

  Vector3D vDiff;
  for(size_t i = 0; i < nCount; i++)
  {
    m_Object.get_colored_point(i, vPoint, fVal, clr);
    vDiff = vPoint - m_vCenter;
    m_vSigma.x += vDiff.x * vDiff.x;
    m_vSigma.y += vDiff.y * vDiff.y;
    m_vSigma.z += vDiff.z * vDiff.z;
  }

  m_vSigma.x = sqrt(m_vSigma.x / nCount);
  m_vSigma.y = sqrt(m_vSigma.y / nCount);
  m_vSigma.z = sqrt(m_vSigma.z / nCount);
}

UINT CTrackCrossSectionCalculator::get_count() const
{
  return m_Object.get_points_count();
}

void CTrackCrossSectionCalculator::update()
{
}

void CTrackCrossSectionCalculator::clear()
{
}

const char* CTrackCrossSectionCalculator::units() const
{
  int nVar = m_Object.get_var();
  switch(nVar)
  {
    case CColoredCrossSection::varIonTemp:     return _T("Ion Temperature, K");
    case CColoredCrossSection::varStartRadius: return _T("Ion Start Radius, mm");
  }

  return _T("");
}

const char* CTrackCrossSectionCalculator::get_var_name(int nVar) const
{
  switch(nVar)
  {
    case CColoredCrossSection::varIonTemp:     return _T("Ion Temperature");
    case CColoredCrossSection::varStartRadius: return _T("Ion Start Radius");
  }

  return _T("");
}

void CTrackCrossSectionCalculator::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CCalculator::save(ar);

  CString cName(m_sOutputFile.c_str());
  ar << cName;

  m_Object.save(ar);
}

void CTrackCrossSectionCalculator::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CCalculator::load(ar);

  CString cName;
  ar >> cName;
  m_sOutputFile = std::string((const char*)cName);

  m_Object.load(ar);
}

//-------------------------------------------------------------------------------------------------
// CSelectedTracksCalculator.
//-------------------------------------------------------------------------------------------------
CSelectedTracksCalculator::CSelectedTracksCalculator()
{
  set_default();
}

CSelectedTracksCalculator::~CSelectedTracksCalculator()
{
}

void CSelectedTracksCalculator::set_default()
{
  m_sCalcName = calc_name(ctAlongSelTracks);

  m_bGasDrag = true;
  m_bDCField = true;
  m_bRFField = false;
  m_bClmb = true;

  m_nSelTrackId = -1;
  m_nSkipPoints = 10;
  default_folder();
}

void CSelectedTracksCalculator::default_folder()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  m_sOutputFolder = (COutputEngine::get_full_path(pObj->get_filename())).c_str();
}

void CSelectedTracksCalculator::run()
{
  if(!get_tracker_ptr())
    return;

  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  CIdsVector vIds = pDrawObj->get_sel_traject_ids();
  size_t nSelCount = vIds.size();
  if(nSelCount == 0)
  {
    AfxMessageBox("No selected trajectories.");
    return;
  }

  char buff[8];
  CString sCount(itoa(nSelCount, buff, 10));
  CString sJobName(_T("Calculating accelerations along selected track(s) "));
  CString sOf(_T(" of "));
  for(size_t i = 0; i < nSelCount; i++)
  {
    CString sCurr = itoa(i + 1, buff, 10);
    m_nSelTrackId = vIds.at(i);
    set_job_name((const char*)(sJobName + sCurr + sOf + sCount));

    do_calculate();
  }
}

void CSelectedTracksCalculator::do_calculate()
{
  CTrackVector& vTracks = m_pObj->get_tracks();
  if(m_nSelTrackId < 0 || m_nSelTrackId >= vTracks.size())
    return;

  char buff[8];
  CString sSelTrack(itoa(m_nSelTrackId, buff, 10)); 
  CString sFileName = m_sOutputFolder + CString(_T("\\")) + CString(_T("Track_#")) + sSelTrack + CString(_T(".csv"));

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (const char*)sFileName, (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  const CTrack& track = vTracks.at(m_nSelTrackId);
  if(m_pObj->get_particle_type() == CTrack::ptIon)
    calc_ion_accel(track, pStream);
  else
    calc_droplet_accel(track, pStream);

  fclose(pStream);
}

void CSelectedTracksCalculator::calc_ion_accel(const CTrack& track, FILE* pStream)
{
  double fTimeStep = m_pObj->get_time_step();
  double fChargeMassRatio = m_pObj->get_particle_charge() / m_pObj->get_ion_mass();

  CElementsCollection& elems = m_pObj->get_elems();
  size_t nElemsCount = elems.size();

  double fCurrent = track.get_current();
  double fPhase = track.get_phase();

  size_t nItemsCount = track.size();
  size_t nStep = 1 + m_nSkipPoints;
  CTrackItem item;
  CNode3D node;

  double fMob, fOneOvrTau, fPow, fExpCoeff, fCoeff;
  Vector3D vGasDrag, vDCField, vRFField, vClmb, vSum, vAccel;
// Header:
  CString s1(_T("time(ms),    x(mm),    y(mm),    z(mm),"));
  CString s2(_T("    GasDrag.x,    GasDrag.y,    GasDrag.z,    DCField.x,    DCField.y,    DCField.z,"));
  CString s3(_T("    RFField.x,    RFField.y,    RFField.z,       Clmb.x,       Clmb.y,       Clmb.z,"));
  CString s4(_T("        Sum.x,        Sum.y,        Sum.z\n\n"));
  fputs((const char*)(s1 + s2 + s3 + s4), pStream);

  const double fAccelCoeff = 1e-5;  // from cm/s^2 to mm/mcs^2.
  for(size_t i = 0; i < nItemsCount; i += nStep)
  {
    track.get_track_item(i, item);
    if(item.nElemId < 0 || item.nElemId >= nElemsCount)
      continue;

    CSymCorrData data = m_pObj->sym_corr_forward(item.pos, item.vel);

    const CElem3D* pElem = elems.at(item.nElemId);
    if(!m_pObj->interpolate(item.pos, item.time, fPhase, node, pElem))
      continue;

    fMob = m_pObj->get_ion_mob(node.press, node.temp);
    fOneOvrTau = fChargeMassRatio / fMob;   // 1 / tau.
    fPow = fTimeStep * fOneOvrTau;          // in fact, fPow = dt / tau.
    fExpCoeff = fPow < 0.01 ? fOneOvrTau : (1. - exp(-fPow)) / fTimeStep;
    fExpCoeff *= fAccelCoeff; // to get acceleration in mm/mcs^2.
    fCoeff = fMob * fExpCoeff;

    if(m_bGasDrag)
      vGasDrag = (node.vel - item.vel) * fExpCoeff;
    if(m_bDCField)
      vDCField = get_DC_field(node) * fCoeff;
    if(m_bRFField)
      vRFField = get_RF_field(node, item.time, fPhase) * fCoeff;
    if(m_bClmb)
      vClmb = space_charge_field(node, item.vel, fCurrent) * fCoeff;

    vSum = vGasDrag + vDCField + vRFField + vClmb;

    if(data.nSymFlag)
    {
      m_pObj->sym_corr_back(item.pos, item.vel, vAccel, data);
      m_pObj->sym_corr_back(vGasDrag, vDCField, vRFField, data);
      m_pObj->sym_corr_back(vClmb, vSum, vAccel, data);
    }

    fprintf(pStream, "%f, %f, %f, %f, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e\n", 
      1000 * item.time, 10 * item.pos.x, 10 * item.pos.y, 10 * item.pos.z,
      vGasDrag.x, vGasDrag.y, vGasDrag.z, vDCField.x, vDCField.y, vDCField.z, vRFField.x, vRFField.y, vRFField.z,
      vClmb.x, vClmb.y, vClmb.z, vSum.x, vSum.y, vSum.z);

    set_progress(100 * i / nItemsCount);
    if(get_terminate_flag())
      break;
  }
}

void CSelectedTracksCalculator::calc_droplet_accel(const CTrack& track, FILE* pStream)
{
}

Vector3D CSelectedTracksCalculator::gas_drag_accel(const Vector3D& vGasVel, const Vector3D& vIonVel, double fExpCoeff) const
{
  return vNull;
}

Vector3D CSelectedTracksCalculator::get_DC_field(const CNode3D& node) const
{
  Vector3D vE = vNull;
  if(m_pObj->get_enable_ansys_field())
  {
    if(m_pObj->get_enable_field())  // enable/disable Ansys DC field:
      vE += node.field;

    vE += m_pObj->get_field_ptb().apply(node.pos);  // DC perturbation field can be swiched on/off individually.
  }
  else
  {
    vE += (node.field + m_pObj->get_field_ptb().apply(node.pos));
  }

  return vE;
}

Vector3D CSelectedTracksCalculator::get_RF_field(const CNode3D& node, double fTime, double fPhase) const
{
  Vector3D vE = vNull;
  if(m_pObj->get_enable_ansys_field())
  {
    if(m_pObj->get_enable_rf())     // RF field:
      vE += m_pObj->get_rf_field(node, fTime, fPhase);
  }
  else
  {
    vE = node.rf;
  }

  return vE;
}

Vector3D CSelectedTracksCalculator::space_charge_field(const CNode3D& node, const Vector3D& vIonVel, double fCurr) const
{
  Vector3D vE = vNull;
  if(m_pObj->get_enable_coulomb())
  {
    if(m_pObj->get_axial_symm())
    {
      if(vIonVel.x > Const_Almost_Zero)
        vE += CSpaceChargeDistrib::radial_coulomb_field(node.pos, vIonVel.x, fCurr);
    }
    else
    {
      Vector3D vClmb = node.clmb;
      if(m_pObj->get_use_radial_coulomb() && (node.pos.x > m_pObj->get_radial_coulomb_trans()))
        vClmb.x = 0;  // This is a workaround. I hope to get rid of in later.

      vE += vClmb;
    }
  }

  return vE;
}

void CSelectedTracksCalculator::update()
{
}

void CSelectedTracksCalculator::clear()
{
}

void CSelectedTracksCalculator::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CCalculator::save(ar);

  ar << m_sOutputFolder;
  ar << m_nSkipPoints;

  ar << m_bGasDrag;
  ar << m_bDCField;
  ar << m_bRFField;
  ar << m_bClmb;
}

void CSelectedTracksCalculator::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CCalculator::load(ar);

  ar >> m_sOutputFolder;
  ar >> m_nSkipPoints;

  ar >> m_bGasDrag;
  ar >> m_bDCField;
  ar >> m_bRFField;
  ar >> m_bClmb;
}

//-------------------------------------------------------------------------------------------------
// CTackFaceCross: calculate faces crossed by tracks in the user-defined range of X
//-------------------------------------------------------------------------------------------------
CTackFaceCross::CTackFaceCross()
{
  set_default();
}

void CTackFaceCross::set_default()
{
  m_sCalcName = calc_name(ctTrackFaceCross);
  m_fStartX = 10.0; // cm
  m_fEndX = 12.0;
}

void CTackFaceCross::run()
{
  if(!get_tracker_ptr())
    return;

	static ThreadPool::Mutex mtx;
	ThreadPool::Locker lock(mtx);
	do_calculate();
}

void CTackFaceCross::do_calculate()
{
	EvaporatingParticle::CTrackVector& tracks = m_pObj->get_tracks();
	const CRegionsCollection& regions = m_pObj->get_regions();
	CFaceIndices newIdxs(tracks.size());
	ThreadPool::splitInPar
	(
		tracks.size(),
		[&](size_t trackId)->void
	{
		if (tracks[trackId].get_term_reason() == CTrack::ttrNone)
		{
			CTrack::const_reverse_iterator last = tracks[trackId].rbegin();
      last++; // we have to take the last but one item because the last item is out of the mesh (see Tracker.cpp, CTracker::main_thread_func)

      Vector3D vPos = (*last)->pos, vVel = (*last)->vel;
      m_pObj->sym_corr_forward(vPos, vVel);
			const CRay ray(vPos, vVel);

			double minDist = DBL_MAX, fDist;
			UINT faceId;
      CFace* pFace = NULL;
			for (size_t regId = 0; regId < regions.size(); ++regId)
			{
				bool ok = regions[regId]->intersect(ray, fDist, faceId);
				if (ok && fDist < minDist)
				{
					minDist = fDist;
          pFace = regions.at(regId)->vFaces.at(faceId);
          if((pFace->box.vMax.x >= m_fStartX) && (pFace->box.vMin.x <= m_fEndX))
          {
            newIdxs[trackId].nReg = regId;
            newIdxs[trackId].nFace = faceId;
          }
				}
			}
		}

	});

	CFaceIndices::iterator _End 
		= std::remove(newIdxs.begin(), newIdxs.end(), CRegFacePair());
	newIdxs.resize(std::distance(newIdxs.begin(), _End));
	CParticleTrackingApp::Get()->GetDrawObj()->append_sel_faces_ids(newIdxs);
}

void CTackFaceCross::update()
{
}

void CTackFaceCross::clear()
{
}

int CTackFaceCross::type() const
{
	return ctTrackFaceCross;
}

void CTackFaceCross::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CCalculator::save(ar);

  ar << m_fStartX;
  ar << m_fEndX;
}

void CTackFaceCross::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CCalculator::load(ar);

  ar >> m_fStartX;
  ar >> m_fEndX;
}

//-------------------------------------------------------------------------------------------------
// CEndTrackCalculator: 
// look for endpoints of the (droplet) tracks that have ended with either "CTrack::ttrRayleighLim" 
// reason or "CTrack::ttrFullyEvapor" reason and output their positions to the file.
//-------------------------------------------------------------------------------------------------
CEndTrackCalculator::CEndTrackCalculator()
{
  set_default();
}

void CEndTrackCalculator::do_calculate()
{
  CTrackVector& vTracks = m_pObj->get_tracks();
  size_t nTrackCount = vTracks.size();

  std::vector<Vector3D> vEndPoints;
  vEndPoints.reserve(nTrackCount);

  for(size_t i = 0; i < nTrackCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    if((track.get_term_reason() != CTrack::ttrRayleighLim) && (track.get_term_reason() != CTrack::ttrFullyEvapor))
      continue;

    size_t nItemCount = track.size();
    vEndPoints.push_back(track.at(nItemCount - 1)->pos);
  }

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (const char*)m_sOutFileName, (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  Vector3D vPos;
  size_t nEndPointCount = vEndPoints.size();
  fprintf(pStream, "%zd\n", nEndPointCount);
  for(size_t j = 0; j < nEndPointCount; j++)
  {
    vPos = vEndPoints.at(j);
    if(m_bSI)
      vPos *= CGS_to_SI_Len;

    fprintf(pStream, "%12.8lf, %12.8lf, %12.8lf\n", vPos.x, vPos.y, vPos.z);
  }

  fclose(pStream);
  vEndPoints.clear();
}

void CEndTrackCalculator::set_default()
{
  m_bSI = false;
  m_sCalcName = calc_name(type());

  if(!get_tracker_ptr())
    return;

  std::string cPath = COutputEngine::get_full_path(m_pObj->get_filename());
  m_sOutFileName = CString(cPath.c_str()) + CString("end_points.csv");
}

void CEndTrackCalculator::run()
{
  if(!get_tracker_ptr())
    return;

  terminate(false);
	do_calculate();
}

void CEndTrackCalculator::update()
{
}

void CEndTrackCalculator::clear()
{
}

int CEndTrackCalculator::type() const
{
	return ctEndTrackCalc;
}

void CEndTrackCalculator::save(CArchive& ar)
{
  UINT nVersion = 1;  // 1 - m_bSI;
  ar << nVersion;

  CCalculator::save(ar);

  ar << m_sOutFileName;
  ar << m_bSI;
}

void CEndTrackCalculator::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CCalculator::load(ar);

  ar >> m_sOutFileName;
  if(nVersion >= 1)
    ar >> m_bSI;
}

};  // namespace EvaporatingParticle
