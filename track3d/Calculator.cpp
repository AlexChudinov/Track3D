
#include "stdafx.h"

#include "float.h"
#include "Calculator.h"

#include "Tracker.hpp"
#include "ParticleTracking.h"

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
  }

  return NULL;
}

const char* CCalculator::calc_name(int nType)
{
  switch(nType)
  {
    case ctPlaneYZ: return _T("Calculator at YZ Plane");
    case ctSelRegions: return _T("Calculator at Selected Regions");
    case ctAlongLine: return _T("Line Calculator");
    case ctTrackCalc: return _T("Ion Track Calculator");
    case ctTrackCrossSect: return _T("Track's Cross-Section Calculator");
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
    case clcAveVx:      fMult = CGS_to_SI_Vel;    break;
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
    case clcAveVx:      return _T("Average Velocity X");
  }

  return _T("");
}

bool CCalculator::get_update_flag() const
{
  return true;
}

void CCalculator::save(CArchive& ar)
{
  UINT nVersion = 1;
  ar << nVersion;

  ar << m_nClcVar;
  ar << m_bEnable;

  ar << m_fCharLength;
}

void CCalculator::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_nClcVar;
  ar >> m_bEnable;

  if(nVersion >= 1)
    ar >> m_fCharLength;
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
  m_bTerminate = false;
  do_sequence_calc();
}

void CPlaneYZCalculator::set_default()
{
  m_sCalcName = calc_name(ctPlaneYZ);

// All positions and lengths are in cm.
  m_fPosX = 1.0;
  m_bReady = false; // force to build the plane mesh.
  m_fMeshSize = 0.05;

// Sequential calculations support:
  m_fStartX = 0.;
  m_fEndX = 1.;
  m_nSeqCalcCount = 100;
}

void CPlaneYZCalculator::clear()
{
  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->set_cross_sections_array(NULL);

  size_t nNodeCount = m_vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
    delete m_vNodes.at(i);

  size_t nFaceCount = m_vFaces.size();
  for(size_t j = 0; j < nFaceCount; j++)
    delete m_vFaces.at(j);

  m_vNodes.clear();
  m_vFaces.clear();
}

void CPlaneYZCalculator::do_calculate()
{
  m_fResult = 0;
  if(!m_bReady && !build_cs_mesh())
    return;

  Vector3D e1, e2, vVel, vNorm(1, 0, 0);
  double fS, fSumS = 0, fRes, fSumRes = 0;
  size_t nFaceCount = m_vFaces.size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    CFace* pFace = m_vFaces.at(i);
    if(!process_face(pFace, fS, fRes))
      continue;   // the face square is less than Const_Almost_Zero.

    fSumRes += fS * fRes;
    fSumS += fS;
  }

  normalize_result(fSumS, fSumRes);
  convert_result();
}

void CPlaneYZCalculator::update()
{
  if(!m_bReady && !build_cs_mesh())
    return;

  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(m_bEnable)
    pDrawObj->set_cross_sections_array(&m_vFaces);
  else
    pDrawObj->set_cross_sections_array(NULL);
}

bool CPlaneYZCalculator::build_cs_mesh()
{
  clear();
  if(!get_tracker_ptr())
    return false;

  CBox& box = m_pObj->get_box();

  double fEps = 0.25 * m_fMeshSize;
  double fYMin = box.vMin.y - fEps, fYMax = box.vMax.y + fEps, fZMin = box.vMin.z - fEps, fZMax = box.vMax.z + fEps;

  double fLenZ = fZMax - fZMin;
  UINT nCountIntervZ = UINT(fLenZ / m_fMeshSize);
  double fStepZ = nCountIntervZ > 0 ? fLenZ / nCountIntervZ : m_fMeshSize;
  UINT nNz = nCountIntervZ + 1;

  double fLenY = fYMax - fYMin;
  UINT nCountIntervY = UINT(fLenY / m_fMeshSize);
  double fStepY = nCountIntervY > 0 ? fLenY / nCountIntervY : m_fMeshSize;
  UINT nNy = nCountIntervY + 1;

// Allocate memory for pointers:
  UINT nPtrsCount = nNy * nNz;
  CNode3D** pNodePtrs = new CNode3D*[nPtrsCount];
  for(UINT k = 0; k < nPtrsCount; k++)
    pNodePtrs[k] = NULL;  // all pointers in the array are initially NULL.

  Vector3D vOrig(m_fPosX, fYMin, 0), vDir(0, 1, 0), vPos;
  for(UINT j = 0; j < nNz; j++)
  {
    vOrig.z = fZMin + j * fStepZ;
    CIntersectColl coll;
    CRay ray(vOrig, vDir);
    if(!intersect(ray, coll) || coll.size() % 2 != 0)   // there must be an even number of intersections.
      continue;

    UINT i0 = UINT((coll.at(0).pos.y - fYMin) / fStepY);  // index of the first, nearest to the ray origin, point of intersection.
    UINT i1 = UINT((coll.at(coll.size() - 1).pos.y - fYMin) / fStepY);  // index of the last, farthest from the ray origin, point of intersection.

    CNode3D node; // this is just a container for interpolated data.
    CElem3D* pElem = NULL;
    for(UINT i = i0; i <= i1; i++)
    {
      vPos = Vector3D(m_fPosX, fYMin + i * fStepY, vOrig.z);
      if(!m_pObj->interpolate(vPos, node, pElem))
        continue;

// Allocate a node inside the 3D mesh:
      UINT k = i + j * nNy;
      pNodePtrs[k] = new CNode3D(node.pos, node.vel, node.field, node.rf, node.dens, node.press, node.temp, node.visc, node.cond, node.cp);
// This index > 0 will be used in process_triangle() to distinguish inner and outer points. The points, which have been initially outer,
// cannot be used as centers of attraction of the other outer points as they are located too close to the boundary.
      pNodePtrs[k]->nInd = k + 1;
    }
  }

  triangulate(nNy, nNz, fStepY, fStepZ, fYMin, fZMin, pNodePtrs);

  for(UINT k = 0; k < nPtrsCount; k++)
  {
    CNode3D* pNode = pNodePtrs[k];
    if(pNode != NULL)
      m_vNodes.push_back(pNode);
  }

  delete[] pNodePtrs; // free only memory allocated for array of pointers, do not delete the pointers themselves.
  m_bReady = true;
  return true;
}

void CPlaneYZCalculator::triangulate(UINT nY, UINT nZ, double fdY, double fdZ, double fY0, double fZ0, CNode3D** pNodePtrs)
{
  UINT k0, k1, k2, k3;
  Vector3D v0(m_fPosX, 0, 0), v1(m_fPosX, 0, 0), v2(m_fPosX, 0, 0), v3(m_fPosX, 0, 0);

  for(UINT j = 0; j < nZ - 1; j++)
  {
    v0.z = fZ0 + j * fdZ;
    v1.z = v0.z;
    v2.z = v0.z + fdZ;
    v3.z = v2.z;

    for(UINT i = 0; i < nY - 1; i++)
    {
      v0.y = fY0 + i * fdY;
      v1.y = v0.y + fdY;
      v2.y = v1.y;
      v3.y = v0.y;

      k0 = i + j * nY;
      k1 = k0 + 1;
      k2 = k1 + nY;
      k3 = k0 + nY;

      process_quad(v0, v1, v2, v3, pNodePtrs[k0], pNodePtrs[k1], pNodePtrs[k2], pNodePtrs[k3]);
    }
  }
}

void CPlaneYZCalculator::process_quad(const Vector3D& v0, const Vector3D& v1, const Vector3D& v2, const Vector3D& v3,
  CNode3D*& p0, CNode3D*& p1, CNode3D*& p2, CNode3D*& p3)
{
  CNode3D* pNodes[4] = { p0, p1, p2, p3 };

  int nGoodNode = -1, nGoodNodesCount = 0;
  for(UINT i = 0; i < 4; i++)
  {
    if((pNodes[i] != NULL))
    {
      nGoodNode = i;
      nGoodNodesCount++;
    }
  }

  if(nGoodNodesCount == 0)
  {
    return;
  }
  else if(nGoodNodesCount == 1)
  {
    switch(nGoodNode) // make sure that there is no triangle, all 3 points of which lie out of the mesh.
    {
      case 0:
      case 2:
      {
        process_triangle(v0, v1, v2, p0, p1, p2);
        process_triangle(v0, v2, v3, p0, p2, p3);
        break;
      }
      case 1:
      case 3:
      {
        process_triangle(v0, v1, v3, p0, p1, p3);
        process_triangle(v1, v2, v3, p1, p2, p3);
        break;
      }
    }
  }
  else
  {
    process_triangle(v0, v1, v3, p0, p1, p3);
    process_triangle(v1, v2, v3, p1, p2, p3);
  }
}

void CPlaneYZCalculator::process_triangle(const Vector3D& v0, const Vector3D& v1, const Vector3D& v2, CNode3D*& p0, CNode3D*& p1, CNode3D*& p2)
{
  Vector3D vVert[3] = { v0, v1, v2 };
  CNode3D* pNodes[3] = { p0, p1, p2 };
  
  UINT i;
  int nGoodNode = -1;
  for(i = 0; i < 3; i++)
  {
    if((pNodes[i] != NULL) && (pNodes[i]->nInd > 0))
    {
      nGoodNode = i;
      break;
    }
  } 

  if(nGoodNode < 0)
    return;   // nothing to do, because all the nodes are out of the mesh.

  CNode3D node; // this is just a container for interpolated data.
  CElem3D* pElem = NULL;
  for(i = 0; i < 3; i++)
  {
    if(pNodes[i] == NULL)
    {
      Vector3D vDest = vVert[i];
      if(!move_node(pNodes[nGoodNode]->pos, vDest))
        return;   // unable to find intersection with the boundary, this triangle cannot be built.

      if(!m_pObj->interpolate(vDest, node, pElem))
        return;   // unable to interpolate data to the new position, this triangle cannot be built.

      CNode3D* pNode = new CNode3D(node.pos, node.vel, node.field, node.rf, node.dens, node.press, node.temp, node.visc, node.cond, node.cp);
      pNode->nInd = 0;   // prohibit this node ever be a center of attraction for other outer nodes.
      switch(i)
      {
        case 0: p0 = pNode; break;
        case 1: p1 = pNode; break;
        case 2: p2 = pNode; break;
      }
    }
  }

  CFace* pFace = new CFace(p0, p1, p2);
  pFace->norm = Vector3D(1, 0, 0);
  m_vFaces.push_back(pFace);
}

bool CPlaneYZCalculator::move_node(const Vector3D& vOrig, Vector3D& vDest) const
{
  Vector3D vDir = (vDest - vOrig).normalized();
  CIntersectColl coll;
  CRay ray(vOrig, vDir);
  if(!intersect(ray, coll))
    return false;

  const CIntersectPoint& point = coll.at(0);
  vDest = vOrig + vDir * (0.99 * point.dist);  // locate the destination a little bit out of the boundary, inside the 3D mesh.

  return true;
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

  std::string sHeader(_T("  x(mm),    "));
  sHeader += units();
  fputs(sHeader.c_str(), pStream);
  fputs("\n\n", pStream);

  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();

  double fCurrPos = m_fPosX;
  double fStepX = (m_fEndX - m_fStartX) / m_nSeqCalcCount;

  for(UINT i = 0; i < m_nSeqCalcCount; i++)
  {
    set_progress(100 * i / m_nSeqCalcCount);
    if(m_bTerminate)
      break;

    set_plane_x(m_fStartX + i * fStepX);
    do_calculate();

    fprintf(pStream, "%f,  %f\n", 10 * m_fPosX, m_fResult);

    if(m_bEnable)
    {
      pDrawObj->set_cross_sections_array(&m_vFaces);
      pDrawObj->draw();
    }
  }

  fclose(pStream);
  if(m_bTerminate)
    return;

  set_progress(100);
  set_plane_x(fCurrPos);
  do_calculate();
}

//---------------------------------------------------------------------------------------
// Streamability:
//---------------------------------------------------------------------------------------
void CPlaneYZCalculator::save(CArchive& ar)
{
  UINT nVersion = 1;
  ar << nVersion;

  CCalculator::save(ar);

  ar << m_fPosX;
  ar << m_fMeshSize;

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

  ar >> m_fPosX;
  ar >> m_fMeshSize;

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
  m_bTerminate = false;
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
      if(m_bTerminate)
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
  m_bTerminate = false;
  do_calculate();
}

void CLineCalculator::set_default()
{
  m_nStepCount = 50;
  m_vLineStart = Vector3D(0, 0.001, 0.001);
  m_vLineEnd = Vector3D(5, 0.001, 0.001);
  m_nClcVar = lcTemp;
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
  CElem3D* pElem = NULL;
  CNode3D node;
  for(UINT i = 0; i <= m_nStepCount; i++)
  {
    set_progress(100 * i / m_nStepCount);
    if(m_bTerminate)
      break;

    vPos = m_vLineStart + vDir * (i * fStep);
    if(!m_pObj->interpolate(vPos, node, pElem))
      continue;

    assign_result(node);

    if(m_nClcVar == lcDens)
      fprintf(pStream, "%f, %f, %f, %e\n", 10 * vPos.x, 10 * vPos.y, 10 * vPos.z, m_fResult);
    else
      fprintf(pStream, "%f, %f, %f, %f\n", 10 * vPos.x, 10 * vPos.y, 10 * vPos.z, m_fResult);
  }

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
    case lcEx:      return _T("Ex,  V/cm");
    case lcEy:      return _T("Ey,  V/cm");
    case lcEz:      return _T("Ez,  V/cm");
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
  m_bTerminate = false;
  do_sequence_calc();
}

void CTrackCalculator::set_default()
{
  m_nClcVar = clcIonTemp;
  m_nCrossSectCount = 50;
  m_fPosCS = 0;
  m_fStartPos = 0;
  m_fEndPos = 1;
}

void CTrackCalculator::do_calculate()
{
  if(!get_tracker_ptr())
    return;

  CTrackVector& vTracks = m_pObj->get_tracks();
  bool bOldIntegr = m_pObj->get_use_old_integrator();
  if(m_nClcVar == clcIonTemp)
    collect_elements();

  m_fResult = m_fIonTempEq = m_fGasTemp = 0;

  Vector3D vPos, vIonVel;
  size_t nIntersectCount = 0;
  size_t nTrackCount = vTracks.size();
  if(nTrackCount < 1)
    return;

  CTrackItem curr, prev;
  double fKsi, fdX, fIonT, fGasT;
  double fCurrIncr = m_pObj->get_full_current() / nTrackCount;
  for(size_t i = 0; i < nTrackCount; i++)
  {
    bool bCalcTime = m_vReachedLastCS.size() == nTrackCount ? m_vReachedLastCS.at(i) : true;  // for time calculation only.

    const CTrack& track = vTracks.at(i);
    size_t nItemCount = track.size();

    for(size_t j = 1; j < nItemCount; j++)
    {
      track.get_track_item(j, curr);
      track.get_track_item(j - 1, prev);
      if((curr.pos.x >= m_fPosCS) && (prev.pos.x < m_fPosCS) || (curr.pos.x < m_fPosCS) && (prev.pos.x >= m_fPosCS))
      {
        fdX = curr.pos.x - prev.pos.x;
        if(fabs(fdX) < Const_Almost_Zero)
          continue;

        fKsi = (m_fPosCS - prev.pos.x) / fdX;

        switch(m_nClcVar)
        {
          case clcIonTemp:
          {
            vPos = prev.pos + fKsi * (curr.pos - prev.pos);
            if(calc_gas_temp(vPos, fGasT))
            {
              m_fResult += prev.temp + fKsi * (curr.temp - prev.temp);
              m_fIonTempEq += prev.tempinf + fKsi * (curr.tempinf - prev.tempinf);
              m_fGasTemp += fGasT;
              nIntersectCount++;
            }

            break;
          }
          case clcCurrent:
          {
            m_fResult += fCurrIncr;
            break;
          }
          case clcTime:
          {
// If this function is called from do_sequence_calc() the vector m_vReachedLastCS contains "true" for those
// tracks which have reached m_fEndPos. It is only those tracks using which the time must be averaged.
            if(bCalcTime)
            {
              m_fResult += prev.time + fKsi * (curr.time - prev.time);
              nIntersectCount++;
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

  if(nIntersectCount != 0)
  {
    double fNormCoeff = 1. / nIntersectCount;
    m_fResult *= fNormCoeff;
    m_fIonTempEq *= fNormCoeff;
    m_fGasTemp *= fNormCoeff;
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

  std::string cAction(_T("Calculation "));
  cAction += std::string(get_var_name(m_nClcVar));
  set_job_name(cAction.c_str());

  std::string sHeader;
  switch(m_nClcVar)
  {
    case clcIonTemp: sHeader = _T("  x(mm),      TI(K),      TIeq(K),    Tgas(K)"); break;
    case clcCurrent: sHeader = _T("  x(mm),   current(nA),   trans(percent)"); break;
    case clcTime:    sHeader = _T("  x(mm),   time(ms)"); break;
  }

  fputs(sHeader.c_str(), pStream);
  fputs("\n\n", pStream);

  if(m_nClcVar == clcTime)
    find_reached_last_cs();

  double fCurrPos = m_fPosCS;
  double fStep = m_nCrossSectCount > 1 ? (m_fEndPos - m_fStartPos) / (m_nCrossSectCount - 1) : 0;
  double fFullCurr = m_pObj->get_full_current() / Const_nA_to_CGSE;   // initial current in nA.
  for(UINT i = 0; i < m_nCrossSectCount; i++)
  {
    set_progress(100 * i / m_nCrossSectCount);
    if(m_bTerminate)
      break;

    m_fPosCS = m_fStartPos + i * fStep;
    do_calculate();

    switch(m_nClcVar)
    {
      case clcIonTemp: fprintf(pStream, "%f,  %f,  %f,  %f\n", 10 * m_fPosCS, m_fResult, m_fIonTempEq, m_fGasTemp); break;
      case clcCurrent: fprintf(pStream, "%f,  %f,  %f\n", 10 * m_fPosCS, m_fResult, 100 * m_fResult / fFullCurr); break;
      case clcTime:    fprintf(pStream, "%f,  %f\n", 10 * m_fPosCS, m_fResult); break;
    }
  }

  fclose(pStream);
  if(m_bTerminate)
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
  m_vReachedLastCS.reserve(nTrackCount);

  CTrackItem item;
  for(size_t i = 0; i < nTrackCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    size_t nItemCount = track.size();

    bool bReached = false;
    for(size_t j = 1; j < nItemCount; j++)
    {
      track.get_track_item(j, item);
      if(item.pos.x >= m_fEndPos)
      {
        bReached = true;
        break;
      }
    }

    m_vReachedLastCS.push_back(bReached);
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
    case clcIonTemp:  return _T("Average Ion Temperature,  K");
    case clcCurrent:  return _T("Ion Current,  nA");
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
    case clcTime:     return _T("Moving Time");
  }

  return _T("");
}

void CTrackCalculator::convert_result()
{
  switch(m_nClcVar)
  {
    case clcCurrent:  m_fResult /= Const_nA_to_CGSE; break;   // nA.
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
  Vector3D vAccel;
  Vector3D vReflPos = vIonPos;
  Vector3D vReflVel = vIonVel;
  m_pObj->sym_corr(vReflPos, vReflVel, vAccel);

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
  Vector3D vVel, vAccel;
  Vector3D vReflPos = vIonPos;
  m_pObj->sym_corr(vReflPos, vVel, vAccel);

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
  m_bTerminate = false;
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

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, m_sOutputFile.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  m_Object.set_handlers(m_hJobNameHandle, m_hProgressBarHandle);

  int nVar = m_Object.get_var();
  switch(nVar)
  {
    case CColoredCrossSection::varIonTemp: fprintf(pStream, "%s\n\n", " z(mm),  y(mm),  Color_Ti,   Ti(K)"); break;
    case CColoredCrossSection::varStartRadius: fprintf(pStream, "%s\n\n", " z(mm),  y(mm),  Color_R0,  R0(cm)"); break;
  }

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
    fprintf(pStream, "%f, %f, %d, %f\n", 10 * vPoint.z, 10 * vPoint.y, RGB(clr.red, clr.green, clr.blue), fVal);
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
    m_vSigma.y += vDiff.y * vDiff.y;
    m_vSigma.z += vDiff.z * vDiff.z;
  }

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

};  // namespace EvaporatingParticle
