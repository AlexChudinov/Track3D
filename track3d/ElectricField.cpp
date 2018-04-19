
#include "stdafx.h"

#include <algorithm>
#include <exception>

#include "ElectricField.h"
#include "ParticleTracking.h"
#include "BarnesHut.h"

#include "../field_solver/DCMeshAdapter.h"


namespace EvaporatingParticle
{

//---------------------------------------------------------------------------------------
// CPotentialBoundCond
//---------------------------------------------------------------------------------------
CPotentialBoundCond::CPotentialBoundCond(BoundaryMesh::BoundaryType type, int val)
  : nType(type), nFixedValType(val)
{
  fStartX = 13.225;
  fStepX = 0.1;
  fEndX = 5.825;
}

const char* CPotentialBoundCond::get_bc_type_name(BoundaryMesh::BoundaryType nType)
{
  switch(nType)
  {
    case BoundaryMesh::FIXED_VAL: return _T("Fixed Value");
    case BoundaryMesh::ZERO_GRAD: return _T("Zero Gradient");
  }

  return _T("None");
}

const char* CPotentialBoundCond::get_fixed_value_name(int nType)
{
  switch(nType)
  {
    case fvPlusUnity:   return _T("+1");
    case fvMinusUnity:  return _T("-1");
    case fvStepLike:    return _T("Stepwise Potential");
    case fvCoulomb:     return _T("Coulomb Potential");
  }

  return _T("0");
}

void CPotentialBoundCond::save(CArchive& ar)
{
  UINT nVersion = 1;  // step-like potential is supported since version 1.
  ar << nVersion;

  int nIntType = (int)nType;
  ar << nIntType;
  ar << nFixedValType;

  CString cStr;
  size_t nCount = vRegNames.size();
  ar << nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    cStr = CString(vRegNames.at(i).c_str());
    ar << cStr;
  }

  ar << fStartX;
  ar << fStepX;
  ar << fEndX;
}

void CPotentialBoundCond::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  int nIntType;
  ar >> nIntType;
  nType = (BoundaryMesh::BoundaryType)nIntType;

  ar >> nFixedValType;

  CString cStr;
  size_t nCount;
  ar >> nCount;
  vRegNames.reserve(nCount);

  for(size_t i = 0; i < nCount; i++)
  {
    ar >> cStr;
    std::string sRegName((const char*)cStr);
    vRegNames.push_back(sRegName);
  }

  if(nVersion >= 1)
  {
    ar >> fStartX;
    ar >> fStepX;
    ar >> fEndX;
  }
}

//-------------------------------------------------------------------------------------------------
// CFieldDataCollection.
//-------------------------------------------------------------------------------------------------
CFieldDataColl::CFieldDataColl()
{
  m_nCurrField = -1;
}

CFieldDataColl::~CFieldDataColl()
{
  clear_fields();
}

void CFieldDataColl::clear_fields()
{
  for(size_t i = 0; i < size(); i++)
  {
    CElectricFieldData* pObj = at(i);
    delete pObj;
  }

  clear();
}

bool CFieldDataColl::calc_fields(bool bMirrorClmb)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  HWND hProgress, hJobName, hDlgWnd;
  pObj->get_handlers(hJobName, hProgress, hDlgWnd);

  if(!bMirrorClmb)
    clear_fields_in_nodes();  // do not clear fields from the nodes if Mirror Coulomb field is about to be calculated.

  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = at(i);
    bool bDoCalc = bMirrorClmb ? pData->get_type() == CElectricFieldData::typeMirror : pData->get_type() != CElectricFieldData::typeMirror;
    if(!bDoCalc)
      continue;

    pData->set_handlers(hJobName, hProgress, hDlgWnd);
    if(!pData->calc_field())
    {
      pData->set_handlers(NULL, NULL, NULL);
      return false;
    }

    pData->set_handlers(NULL, NULL, NULL);
  }

  return true;
}

static const Vector3D scvNull(0, 0, 0);

void CFieldDataColl::clear_fields_in_nodes()
{
  CNodesCollection& vNodes = CParticleTrackingApp::Get()->GetTracker()->get_nodes();

  CNode3D* pNode = NULL;
  size_t nCount = vNodes.size();
  for(size_t i = 0; i < nCount; i++)
  {
    pNode = vNodes.at(i);
    pNode->field = scvNull;
    pNode->rf = scvNull;
  }
}

bool CFieldDataColl::sel_region_changed(CStringVector* pRegNames)
{
  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = at(i);
    size_t nBoundCondCount = pData->get_bc_count();
    for(size_t j = 0; j < nBoundCondCount; j++)
    {
      CPotentialBoundCond* pBC = pData->get_bc(j);
      if(pRegNames == &(pBC->vRegNames))
      {
        pData->invalidate();
        return true;
      }
    }
  }

  return false;
}

bool CFieldDataColl::remove_bound_cond(CPotentialBoundCond* pBC)
{
  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = at(i);
    size_t nBoundCondCount = pData->get_bc_count();
    for(size_t j = 0; j < nBoundCondCount; j++)
    {
      if(pBC == pData->get_bc(j))
      {
        pData->remove_bc(j);
        return true;
      }
    }
  }

  return false;
}

void CFieldDataColl::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  size_t nCount = size();
  ar << nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = at(i);
    pData->save(ar);
  }
}

void CFieldDataColl::load(CArchive& ar)
{
  clear_fields();

  UINT nVersion;
  ar >> nVersion;

  size_t nCount;
  ar >> nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = new CElectricFieldData();
    pData->load(ar);
    push_back(pData);
  }

  m_nCurrField = size() > 0 ? 0 : -1;
}

//---------------------------------------------------------------------------------------
// CElectricFieldData - an auxiliary class for simulating scalable electric fields.
//---------------------------------------------------------------------------------------
CElectricFieldData::CElectricFieldData(int nType)
  : m_nType(nType)
{
  set_default();
}

CElectricFieldData::~CElectricFieldData()
{
  clear_bc();
}

void CElectricFieldData::set_default()
{
  m_bEnable = true;
  m_nCalcMethod = cmLaplacian3;   // the oldest and most tested method so far.

  m_fScale = 10 * SI_to_CGS_Voltage;  // 10 V in CGS.
  set_freq(1.0e+6);                   // 1 MHz by default.

  m_bNeedRecalc = true;
  m_nIterCount = 100;

// An attempt to get analytic field in the flatapole. Alpha version.
  m_bAnalytField = false;
  m_fRadius = 0.21;   // cm, inscribed radius of the flatapole electrodes.
  m_fLowLimX = 14.7;  // cm, an analytic formula will be used if m_fLowLimX < x < m_fHighLimX;
  m_fHighLimX = 1e+5;

// Default name of the field:
  char buff[4];
  int nFieldsCount = CParticleTrackingApp::Get()->GetFields()->size();
  m_sName = CString("Electric Field # ") + CString(itoa(nFieldsCount + 1, buff, 10));
}

bool CElectricFieldData::calc_field(bool bTest)
{
  if(!m_bEnable)
    return true;

  if(!m_bNeedRecalc)
    return get_result(bTest);

  switch(m_nCalcMethod)
  {
    case cmLaplacian3: return calc_lap3(bTest);
    case cmDirTessLap3: return calc_dirichlet_lap3(bTest);
    case cmFinVolSeidel: return calc_finite_vol_seidel(bTest);
  }

  return false;
}

bool CElectricFieldData::calc_lap3(bool bTest)
{
  try
  {
    m_bTerminate = false;
    CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
    const CNodesCollection& vNodes = pObj->get_nodes();
    CMeshAdapter mesh(pObj->get_elems(), vNodes);
    mesh.progressBar()->set_handlers(m_hJobNameHandle, m_hProgressBarHandle, m_hDlgWndHandle);

    if(!set_default_boundary_conditions(mesh) || !set_boundary_conditions(mesh))
      return false;

    size_t nNodesCount = vNodes.size();
    std::vector<double> field(nNodesCount, 0.0);

// Laplacian Solver:
    mesh.boundaryMesh()->applyBoundaryVals(field);
    CMeshAdapter::PScalFieldOp pOp = mesh.createOperator(CMeshAdapter::LaplacianSolver3);
    if(pOp == NULL)
      return false;

    set_job_name("Field calculation...");
    for(int i = 0; i < m_nIterCount; ++i)
    {
      if(m_bTerminate)
        return false;

      field = pOp->applyToField(field);
      set_progress(100 * i / m_nIterCount);
    }

// Gradient calculation:
    CMeshAdapter::PScalFieldOp pGrad = mesh.createOperator(CMeshAdapter::GradX);
    std::vector<double> dPhiDx(nNodesCount, 0.0);
    dPhiDx = pGrad->applyToField(field);

    pGrad = mesh.createOperator(CMeshAdapter::GradY);
    std::vector<double> dPhiDy(nNodesCount, 0.0);
    dPhiDy = pGrad->applyToField(field);

    pGrad = mesh.createOperator(CMeshAdapter::GradZ);
    std::vector<double> dPhiDz(nNodesCount, 0.0);
    dPhiDz = pGrad->applyToField(field);

// DEBUG  (MS 04-04-2018 Mirror Coulomb potential testing)
    CBarnesHut* pBHObj = pObj->get_BH_object();
    CNode3D* pNode = NULL;
// END DEBUG

    m_vField.resize(nNodesCount);
    for(size_t i = 0; i < nNodesCount; i++)
    {
      m_vField[i] = Vector3F(-(float)dPhiDx[i], -(float)dPhiDy[i], -(float)dPhiDz[i]);

// DEBUG  (MS 04-04-2018 Mirror Coulomb potential visualization).
      pNode = CParticleTrackingApp::Get()->GetTracker()->get_nodes().at(i);
      pNode->phi = (float)field[i];
// END DEBUG

// An attempt to get analytic field in the flatapole. Alpha version.
      if(m_bAnalytField && (m_nType == typeFieldRF))
        apply_analytic_field(vNodes.at(i)->pos, m_vField[i]);
    }

    m_bNeedRecalc = m_nType == typeMirror;  // false;
	}
  catch (const std::exception& ex)
  {
    AfxMessageBox(ex.what());
  }

  bool bRes = get_result(bTest);
  notify_scene();

  /*[AC 03.05.2017] Memory clean up */
  BlockPoolInterface::cleanUpEveryPool();
  /*[/AC]*/

  return bRes;
}

bool CElectricFieldData::calc_dirichlet_lap3(bool bTest)
{
  CDirichletTesselation* pTessObj = CParticleTrackingApp::Get()->GetDirichletTess();
  try
  {
    m_bTerminate = false;
    CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
    const CNodesCollection& vNodes = pObj->get_nodes();
    DCMeshAdapter mesh(pObj->get_elems(), vNodes, *pTessObj);
    mesh.progressBar()->set_handlers(m_hJobNameHandle, m_hProgressBarHandle, m_hDlgWndHandle);

    if(!set_default_boundary_conditions(mesh) || !set_boundary_conditions(mesh))
      return false;

    size_t nNodesCount = vNodes.size();
    std::vector<double> field(nNodesCount, 0.0);

// Laplacian Solver:
    mesh.boundaryMesh()->applyBoundaryVals(field);
    DCMeshAdapter::PScalFieldOp pOp = mesh.createOperator(DCMeshAdapter::LaplacianSolver, NULL);
    if(pOp == NULL)
      return false;

    set_job_name("Field calculation...");
    for(int i = 0; i < m_nIterCount; ++i)
    {
      if(m_bTerminate)
        return false;

      field = pOp->applyToField(field);
      set_progress(100 * i / m_nIterCount);
    }

// Gradient calculation:
    DCMeshAdapter::PScalFieldOp pGrad = mesh.createOperator(DCMeshAdapter::GradX, NULL);
    std::vector<double> dPhiDx(nNodesCount, 0.0);
    dPhiDx = pGrad->applyToField(field);

    pGrad = mesh.createOperator(DCMeshAdapter::GradY, NULL);
    std::vector<double> dPhiDy(nNodesCount, 0.0);
    dPhiDy = pGrad->applyToField(field);

    pGrad = mesh.createOperator(DCMeshAdapter::GradZ, NULL);
    std::vector<double> dPhiDz(nNodesCount, 0.0);
    dPhiDz = pGrad->applyToField(field);

// DEBUG  (MS 04-04-2018 Mirror Coulomb potential testing)
    CBarnesHut* pBHObj = pObj->get_BH_object();
    CNode3D* pNode = NULL;
// END DEBUG

    m_vField.resize(nNodesCount);
    for(size_t i = 0; i < nNodesCount; i++)
    {
      m_vField[i] = Vector3F(-(float)dPhiDx[i], -(float)dPhiDy[i], -(float)dPhiDz[i]);

// DEBUG  (MS 04-04-2018 Mirror Coulomb potential visualization).
      pNode = CParticleTrackingApp::Get()->GetTracker()->get_nodes().at(i);
      pNode->phi = (float)field[i];
// END DEBUG

// An attempt to get analytic field in the flatapole. Alpha version.
      if(m_bAnalytField && (m_nType == typeFieldRF))
        apply_analytic_field(vNodes.at(i)->pos, m_vField[i]);
    }

    m_bNeedRecalc = m_nType == typeMirror;  // false;
	}
  catch (const std::exception& ex)
  {
    AfxMessageBox(ex.what());
  }

  bool bRes = get_result(bTest);
  notify_scene();

  /*[AC 03.05.2017] Memory clean up */
  BlockPoolInterface::cleanUpEveryPool();
  /*[/AC]*/

  return bRes;
}

bool CElectricFieldData::calc_finite_vol_seidel(bool bTest)
{
  return true;
}


bool CElectricFieldData::get_result(bool bTest) const
{
  CNodesCollection& vNodes = CParticleTrackingApp::Get()->GetTracker()->get_nodes();
  size_t nNodeCount = vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
//    if(bTest)
//      vNodes.at(i)->phi = m_vField[i].length();
    if(m_nType == typeFieldDC)
      vNodes.at(i)->field += Vector3D(m_vField[i].x, m_vField[i].y, m_vField[i].z) * m_fScale;
    else if(m_nType == typeMirror)
      vNodes.at(i)->clmb += Vector3D(m_vField[i].x, m_vField[i].y, m_vField[i].z);
  }

  return true;
}

void CElectricFieldData::apply_analytic_field(const Vector3D& vPos, Vector3F& vField)
{
  if(vPos.x < m_fLowLimX || vPos.x > m_fHighLimX)
    return;

  double fR = sqrt(vPos.y * vPos.y + vPos.z * vPos.z);
  if(fR > m_fRadius)
    return;

  double fCoeff = 2. / (m_fRadius * m_fRadius);
  Vector3F vAnalytRF(0, -fCoeff * vPos.y, fCoeff * vPos.z);

  const double scfTransWidth = 0.2;
  double fLowLimMax = m_fLowLimX + scfTransWidth; // the field changes gradually from numeric to analytic; transition interval is 2 mm long.
  if((vPos.x >= m_fLowLimX) && (vPos.x <= fLowLimMax))
  {
    float fKsi = float((vPos.x - m_fLowLimX) / scfTransWidth);
    vField = vAnalytRF * fKsi + vField * (1 - fKsi);
  }

  vField = vAnalytRF;
}

bool CElectricFieldData::set_boundary_conditions(CMeshAdapter& mesh)
{
  set_job_name("Setting boundary conditions...");
  set_progress(0);

  CRegion* pReg = NULL;
  size_t i, j, nBoundCondCount = m_vBoundCond.size(), nRegCount, nSumRegCount = 0;
  for(i = 0; i < nBoundCondCount; i++)
    nSumRegCount += m_vBoundCond.at(i)->vRegNames.size();

  size_t k = 0;
  for(i = 0; i < nBoundCondCount; i++)
  {
    CPotentialBoundCond* pBC = m_vBoundCond.at(i);
    nRegCount = pBC->vRegNames.size();
    for(j = 0; j < nRegCount; j++)
    {
      if(m_bTerminate)
        return false;

      k++;
      const std::string& sName = pBC->vRegNames.at(j);
      pReg = get_region(sName);
      if(pReg == NULL)
        continue;

      set_boundary_values(mesh, pReg, pBC);

      set_progress(100 * k / nSumRegCount);
    }
  }

  return true;
}

bool CElectricFieldData::set_default_boundary_conditions(CMeshAdapter& mesh)
{
  set_job_name("Setting default boundary conditions...");
  set_progress(0);

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CRegionsCollection& vRegions = pObj->get_regions();
  size_t nRegCount = vRegions.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    if(m_bTerminate)
      return false;

    CRegion* pReg = vRegions.at(i);
    if(is_selected(pReg) || pReg->bCrossSection)
      continue;

    set_boundary_values(mesh, pReg);

    set_progress(100 * i / nRegCount);
  }

  return true;
}

void CElectricFieldData::set_boundary_values(CMeshAdapter& mesh, CRegion* pReg, CPotentialBoundCond* pBC)
{
  CFace* pFace = NULL;
  std::vector<UINT> vNodeIds;
  std::vector<Vector3D> vNorms;
  size_t nBndFacesCount = pReg->vFaces.size();
  if(nBndFacesCount == 0)
    return;

  for(size_t k = 0; k < nBndFacesCount; k++)
  {
    pFace = pReg->vFaces.at(k);
    vNodeIds.push_back(pFace->p0->nInd);
    vNorms.push_back(calc_norm(pFace->p0));

    vNodeIds.push_back(pFace->p1->nInd);
    vNorms.push_back(calc_norm(pFace->p1));

    vNodeIds.push_back(pFace->p2->nInd);
    vNorms.push_back(calc_norm(pFace->p2));
  }

  const std::string& sName = pReg->sName;
  mesh.boundaryMesh()->addBoundary(sName, vNodeIds, vNorms);

  BoundaryMesh::BoundaryType nType = pBC != NULL ? pBC->nType : BoundaryMesh::FIXED_VAL;
  mesh.boundaryMesh()->boundaryType(sName, nType);

  if(nType != BoundaryMesh::FIXED_VAL)
    return;

// There are two cases: 
// 1) Ordinary field specified by a single constant potential value at all nodes of a region (metallic electrode);
// 2) Mirror Coulomb field, boundary potential values for which can be different at different nodes of a region.
  double fVal = 0;
  if((pBC != NULL) && (nType == BoundaryMesh::FIXED_VAL))
  {
    if(pBC->nFixedValType == CPotentialBoundCond::fvCoulomb)  // mirror Coulomb field.
    {
      CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
      CBarnesHut* pBHObj = pObj->get_BH_object();

      const BoundaryMesh::SetLabels& vNodeLabels = mesh.boundaryMesh()->boundaryLabels(sName);
      size_t nRegNodeCount = vNodeLabels.size();
      std::vector<double> vBoundPhi;
      BoundaryMesh::SetLabels::const_iterator citer;
      UINT nNodeId = 0;
      Vector3D vPos;
      for(citer = vNodeLabels.cbegin(); citer != vNodeLabels.cend(); citer++)
      {
        nNodeId = *citer;
        vPos = pObj->get_nodes().at(nNodeId)->pos;
        fVal = -pBHObj->coulomb_phi(vPos);
        vBoundPhi.push_back(fVal);
      }

      mesh.boundaryMesh()->boundaryVals(sName, vBoundPhi);
      return;
    }

    switch(pBC->nFixedValType)
    {
      case CPotentialBoundCond::fvPlusUnity: fVal = 1.0; break;
      case CPotentialBoundCond::fvMinusUnity: fVal = -1.0; break;
      case CPotentialBoundCond::fvStepLike: fVal = step_potential(pBC, pReg->vFaces.at(0)->p0->pos); break;
    }
  }
    
  mesh.boundaryMesh()->boundaryVals(
	  sName, 
	  std::vector<double>(mesh.boundaryMesh()->patchSize(sName), fVal));
}

double CElectricFieldData::step_potential(CPotentialBoundCond* pBC, const Vector3D& vPos) const
{
  double fX = vPos.x;
  int nStepsCount = int(fabs(pBC->fEndX - pBC->fStartX) / pBC->fStepX);
  double fPhiStep = 1. / nStepsCount;
  if(pBC->fStartX < pBC->fEndX)
  {
    if(fX <= pBC->fStartX)
      return 0;
    else if(fX <= pBC->fEndX)
      return fPhiStep * int((fX - pBC->fStartX) / pBC->fStepX);
    else
      return 1;
  }
  else
  {
    if(fX >= pBC->fStartX)
      return 0;
    else if(fX >= pBC->fEndX)
      return fPhiStep * int((pBC->fStartX - fX) / pBC->fStepX);
    else
      return 1;
  }
}

CRegion* CElectricFieldData::get_region(const std::string& sName) const
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CRegionsCollection& vRegions = pObj->get_regions();
  size_t nRegCount = vRegions.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    CRegion* pReg = vRegions.at(i);
    if(sName == pReg->sName)
      return pReg;
  }

  return NULL;
}

Vector3D CElectricFieldData::calc_norm(CNode3D* pNode) const
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CRegionsCollection& vRegs = pObj->get_regions();
  size_t nRegCount = vRegs.size();

  UINT nReg, nFace;
  CFace* pFace = NULL;
  Vector3D vNorm(0, 0, 0);
  size_t nNbrCount = pNode->vNbrFaces.size();
  for(size_t i = 0; i < nNbrCount; i++)
  {
    nReg = pNode->vNbrFaces.at(i).nReg;
    if(nReg >= nRegCount)
      continue;

    nFace = pNode->vNbrFaces.at(i).nFace;
    if(nFace >= vRegs.at(nReg)->vFaces.size())
      continue;

    pFace = vRegs.at(nReg)->vFaces.at(nFace);
    vNorm += pFace->norm;
  }

  vNorm.normalize();
  return -vNorm;
}

bool CElectricFieldData::is_selected(CRegion* pReg) const
{
  size_t nBoundCondCount = m_vBoundCond.size();
  for(size_t i = 0; i < nBoundCondCount; i++)
  {
    CPotentialBoundCond* pBC = m_vBoundCond.at(i);
    size_t nRegCount = pBC->vRegNames.size();
    for(size_t j = 0; j < nRegCount; j++)
    {
      const std::string& sName = pBC->vRegNames.at(j);
      if(pReg->sName == sName)
        return true;
    }
  }

  return false;
}

const char* CElectricFieldData::get_field_type_name(int nType)
{
  switch(nType)
  {
    case typeFieldDC: return _T("DC Electric Field");
    case typeFieldRF: return _T("Radio-Frequency Field");
    case typeMirror:  return _T("Mirror Coulomb Field");
  }

  return _T("None");
}

const char* CElectricFieldData::get_calc_method_name(int nCalcMethod)
{
  switch(nCalcMethod)
  {
    case cmLaplacian3: return _T("Laplacian Solver #3");
    case cmDirTessLap3: return _T("Advanced Laplacian Solver");
    case cmFinVolSeidel:  return _T("Finite Volume (Gauss-Seidel)");
  }

  return _T("None");
}

void CElectricFieldData::clear_bc()
{
  size_t nCount = m_vBoundCond.size();
  for(size_t i = 0; i < nCount; i++)
  {
    CPotentialBoundCond* pObj = m_vBoundCond.at(i);
    delete pObj;
  }

  m_vBoundCond.clear();
}

void CElectricFieldData::add_bc()
{
  char buff[4];
  size_t nId = m_vBoundCond.size();
  CString cName = CString("Boundary Conditions # ") + CString(itoa(nId + 1, buff, 10));

  CPotentialBoundCond* pBC = new CPotentialBoundCond();
  pBC->sName = std::string(cName);
  m_vBoundCond.push_back(pBC);
}

void CElectricFieldData::remove_bc(size_t nId)
{
  if(nId < m_vBoundCond.size())
  {
    CPotentialBoundCond* pBC = m_vBoundCond.at(nId);
    m_vBoundCond.erase(m_vBoundCond.begin() + nId);
    delete pBC;
  }
}

void CElectricFieldData::notify_scene()
{
  CCrossSectColl* pColl = CParticleTrackingApp::Get()->GetPlanes();
  for(size_t i = 0; i < pColl->size(); i++)
    pColl->at(i)->invalidate();
}

void CElectricFieldData::save(CArchive& ar)
{
  UINT nVersion = 3;  // 3 - Calculation method; 2 - Analytic RF field kitchen, alpha version; 1 - m_bEnable and calculated
  ar << nVersion;

  ar << m_bEnable;
  ar << m_nType;
  ar << m_fScale;
  ar << m_fOmega;
  ar << m_nIterCount;
  ar << m_nCalcMethod;

  size_t nCount = m_vBoundCond.size();
  ar << nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    CPotentialBoundCond* pObj = m_vBoundCond.at(i);
    pObj->save(ar);
  }

  size_t nNodesCount = m_bNeedRecalc ? 0 : m_vField.size();
  ar << nNodesCount;
  for(size_t j = 0; j < nNodesCount; j++)
  {
    ar << m_vField[j].x;
    ar << m_vField[j].y;
    ar << m_vField[j].z;
  }

  ar << m_bAnalytField;
  ar << m_fRadius;      // inscribed radius of the flatapole electrodes.
  ar << m_fLowLimX;     // an analytic formula will be used if m_fLowLimX < x < m_fHighLimX;
  ar << m_fHighLimX;
  ar << m_sName;
}

void CElectricFieldData::load(CArchive& ar)
{
  clear_bc();

  UINT nVersion;
  ar >> nVersion;

  if(nVersion >= 1)
    ar >> m_bEnable;

  ar >> m_nType;
  ar >> m_fScale;
  ar >> m_fOmega;
  ar >> m_nIterCount;

  if(nVersion >= 3)
    ar >> m_nCalcMethod;

  size_t nCount;
  ar >> nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    add_bc();
    CPotentialBoundCond* pObj = m_vBoundCond.at(m_vBoundCond.size() - 1);
    pObj->load(ar);
  }

  if(nVersion >= 1)
  {
    size_t nNodesCount;
    ar >> nNodesCount;
    if(nNodesCount > 0)
      m_vField.resize(nNodesCount, Vector3F(0, 0, 0));

    for(size_t j = 0; j < nNodesCount; j++)
    {
      ar >> m_vField[j].x;
      ar >> m_vField[j].y;
      ar >> m_vField[j].z;
    }

    if(nNodesCount > 0)
      m_bNeedRecalc = false;
  }

  if(nVersion >= 2)
  {
    ar >> m_bAnalytField;
    ar >> m_fRadius;      // inscribed radius of the flatapole electrodes.
    ar >> m_fLowLimX;     // an analytic formula will be used if m_fLowLimX < x < m_fHighLimX;
    ar >> m_fHighLimX;
    ar >> m_sName;
  }
}

};  // namespace EvaporatingParticle.