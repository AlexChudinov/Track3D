
#include "stdafx.h"

#include <algorithm>
#include <exception>

#include "ElectricField.h"
#include "ParticleTracking.h"


namespace EvaporatingParticle
{

//---------------------------------------------------------------------------------------
// CPotentialBoundCond
//---------------------------------------------------------------------------------------
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
    case fvPlusUnity: return _T("+1");
    case fvMinusUnity: return _T("-1");
  }

  return _T("0");
}

void CPotentialBoundCond::save(CArchive& ar)
{
  UINT nVersion = 0;
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
}

//-------------------------------------------------------------------------------------------------
// CFieldDataCollection.
//-------------------------------------------------------------------------------------------------
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

bool CFieldDataColl::calc_fields()
{
  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = at(i);
    if(!pData->calc_field())
      return false;
  }

  return true;
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

  m_fScale = 10;    // V.
  set_freq(1.0e+6); // 1 MHz by default.

  m_bNeedRecalc = true;
  m_bScaleChanged = true;
  m_nIterCount = 100;
}

bool CElectricFieldData::calc_field()
{
  if(!m_bEnable && !m_bNeedRecalc)
    return true;

  try
  {
    m_bTerminate = false;
    CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
    CMeshAdapter mesh(pObj->get_elems(), pObj->get_nodes());
    mesh.progressBar()->set_handlers(m_hJobNameHandle, m_hProgressBarHandle, m_hDlgWndHandle);

    if(!set_default_boundary_conditions(mesh) || !set_boundary_conditions(mesh))
      return false;

    std::vector<double> field(pObj->get_nodes().size(), 0.0);
    mesh.boundaryMesh()->applyBoundaryVals(field);

	CMeshAdapter::PScalFieldOp pDx = mesh.createOperator(CMeshAdapter::GradX);
	CMeshAdapter::PScalFieldOp pDy = mesh.createOperator(CMeshAdapter::GradY);
	CMeshAdapter::PScalFieldOp pDz = mesh.createOperator(CMeshAdapter::GradZ);

    CMeshAdapter::PScalFieldOp pOp = mesh.createOperator(CMeshAdapter::LaplacianSolver1);
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

    std::vector<double> dPhiDx(pObj->get_nodes().size(), 0.0);
    dPhiDx = pDx->applyToField(field);
    std::vector<double> dPhiDy(pObj->get_nodes().size(), 0.0);
    dPhiDy = pDy->applyToField(field);
    std::vector<double> dPhiDz(pObj->get_nodes().size(), 0.0);
    dPhiDz = pDz->applyToField(field);
    for(size_t i = 0; i < pObj->get_nodes().size(); i++)
      dPhiDx[i] = sqrt(dPhiDx[i] * dPhiDx[i] + dPhiDy[i] * dPhiDy[i] + dPhiDz[i] * dPhiDz[i]);

    if(!m_bTerminate)
    {
      get_result(dPhiDx);
      notify_scene();

      m_bNeedRecalc = false;
		}
	}
  catch (const std::exception& ex)
  {
    AfxMessageBox(ex.what());
  }

  return true;
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
      return false;;

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

  double fVal = 0;
  if((pBC != NULL) && (nType == BoundaryMesh::FIXED_VAL))
    fVal = pBC->nFixedValType == CPotentialBoundCond::fvPlusUnity ? 1 : -1;
    
  mesh.boundaryMesh()->boundaryVals(
	  sName, 
	  std::vector<double>(mesh.boundaryMesh()->patchSize(sName), fVal));
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

void CElectricFieldData::get_result(const std::vector<double>& vFieldVals)
{
	const CNodesCollection& vNodes = CParticleTrackingApp::Get()->GetTracker()->get_nodes();
	for (size_t i = 0; i < vNodes.size(); i++)
	{
		CNode3D* pNode = vNodes[i];
		pNode->phi = (float)vFieldVals[i];
	}
}

const char* CElectricFieldData::get_field_type_name(int nType)
{
  switch(nType)
  {
    case typeFieldDC: return _T("DC Electric Field");
    case typeFieldRF: return _T("Radio-Frequency Field");
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
  UINT nVersion = 0;
  ar << nVersion;

  ar << m_nType;
  ar << m_fScale;
  ar << m_fOmega;
  ar << m_nIterCount;

  size_t nCount = m_vBoundCond.size();
  ar << nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    CPotentialBoundCond* pObj = m_vBoundCond.at(i);
    pObj->save(ar);
  }
}

void CElectricFieldData::load(CArchive& ar)
{
  clear_bc();

  UINT nVersion;
  ar >> nVersion;

  ar >> m_nType;
  ar >> m_fScale;
  ar >> m_fOmega;
  ar >> m_nIterCount;

  size_t nCount;
  ar >> nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    add_bc();
    CPotentialBoundCond* pObj = m_vBoundCond.at(m_vBoundCond.size() - 1);
    pObj->load(ar);
  }
}

};  // namespace EvaporatingParticle.