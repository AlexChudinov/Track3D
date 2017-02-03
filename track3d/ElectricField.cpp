
#include <exception>

#include "stdafx.h"
#include "ElectricField.h"
#include "ParticleTracking.h"


namespace EvaporatingParticle
{

//---------------------------------------------------------------------------------------
// CPotentialBoundCond
//---------------------------------------------------------------------------------------
const char* CPotentialBoundCond::get_bc_type_name(PotentialField::BOUNDARY_TYPE nType)
{
  switch(nType)
  {
    case PotentialField::FIXED_VAL: return _T("Fixed Value");
    case PotentialField::ZERO_GRAD: return _T("Zero Gradient");
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
  nType = (PotentialField::BOUNDARY_TYPE)nIntType;

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

void CFieldDataColl::calc_fields()
{
  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = at(i);
    pData->calc_field();
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
  m_fScale = 10;    // V.
  set_freq(1.0e+6); // 1 MHz by default.

  m_bNeedRecalc = true;
  m_bScaleChanged = true;
  m_nIterCount = 100;
}

void CElectricFieldData::calc_field()
{
	try
	{
		Mesh* pMesh = create_mesh();
		if (pMesh == NULL)
		{
			AfxMessageBox(_T("NULL pointer is returned by create_mesh()."));
			return;
		}

		PotentialField* pField = PotentialField::createZeros(pMesh);
		Mesh::free(pMesh);

    set_default_boundary_conditions(pField);
		set_boundary_conditions(pField);

    pField->applyBoundaryConditions();
		ScalarFieldOperator* pSolver = ScalarFieldOperator::create(pField, ScalarFieldOperator::LaplacianSolver);

		set_job_name("Calculating field...");
		for(UINT i = 0; i < m_nIterCount; i++)
		{
			set_progress(100 * i / m_nIterCount);
			if (m_bTerminate)
				break;

			pSolver->applyToField(pField);
		}

		if(!m_bTerminate)
    {
			if (get_result(pField))
				notify_scene();
			else
				AfxMessageBox(_T("Size of result differs from the nodes count."));
		}

		PotentialField::free(pField);
	}
	catch (const std::exception& ex)
	{
		AfxMessageBox(ex.what());
	}
}

Mesh* CElectricFieldData::create_mesh()
{
  set_job_name("Creating mesh...");
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  const CNodesCollection& vNodes = pObj->get_nodes();
  size_t nNodeCount = vNodes.size();

  CElementsCollection& vElems = pObj->get_elems();
  size_t nElemCount = vElems.size();

  Graph* g = Graph::create();

  Vector3D v;
  std::vector<V3D> vNodePos(nNodeCount);
  for(size_t j = 0; j < nNodeCount; j++)
  {
    v = vNodes.at(j)->pos;
    vNodePos[j].x = v.x;
    vNodePos[j].y = v.y;
    vNodePos[j].z = v.z;
  }

  CElem3D* pElem = NULL;
  for(size_t i = 0; i < nElemCount; i++)
  {
    set_progress(100 * i / nElemCount);
    if(m_bTerminate)
      break;

    pElem = vElems.at(i);
    const CNodesCollection& vElNodes = pElem->vNodes;
    switch(pElem->vNodes.size())
    {
      case 4: g->addTet(vElNodes[0]->nInd, vElNodes[1]->nInd, vElNodes[2]->nInd, vElNodes[3]->nInd); break;
      case 5: g->addPyr(vElNodes[0]->nInd, vElNodes[1]->nInd, vElNodes[2]->nInd, vElNodes[3]->nInd, vElNodes[4]->nInd); break;
      case 6: g->addWedge(vElNodes[0]->nInd, vElNodes[1]->nInd, vElNodes[2]->nInd, vElNodes[3]->nInd, vElNodes[4]->nInd, vElNodes[5]->nInd); break;
      case 8: g->addHexa(vElNodes[0]->nInd, vElNodes[1]->nInd, vElNodes[3]->nInd, vElNodes[2]->nInd, vElNodes[4]->nInd, vElNodes[5]->nInd, vElNodes[7]->nInd, vElNodes[6]->nInd); break;
    }
  }

  Mesh* pMesh = NULL;
  if(!m_bTerminate)
    pMesh = Mesh::create(g, vNodePos);

  Graph::free(g);
  return pMesh;
}

void CElectricFieldData::set_boundary_conditions(PotentialField* pField)
{
  set_job_name("Setting boundary conditions...");

  CRegion* pReg = NULL;
  size_t nBoundCondCount = m_vBoundCond.size();
  for(size_t i = 0; i < nBoundCondCount; i++)
  {
    CPotentialBoundCond* pBC = m_vBoundCond.at(i);
    size_t nRegCount = pBC->vRegNames.size();
    for(size_t j = 0; j < nRegCount; j++)
    {
      const std::string& sName = pBC->vRegNames.at(j);
      pReg = get_region(sName);
      if(pReg == NULL)
        continue;

      set_boundary_values(pField, pReg, pBC);
    }
  }
}

void CElectricFieldData::set_default_boundary_conditions(PotentialField* pField)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CRegionsCollection& vRegions = pObj->get_regions();
  size_t nRegCount = vRegions.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    CRegion* pReg = vRegions.at(i);
    if(is_selected(pReg) || pReg->bCrossSection)
      continue;

    set_boundary_values(pField, pReg);
  }
}

void CElectricFieldData::set_boundary_values(PotentialField* pField, CRegion* pReg, CPotentialBoundCond* pBC)
{
  CFace* pFace = NULL;
  std::vector<UINT> vNodeIds;
  std::vector<V3D> vNorms;
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
  pField->addBoundary(sName, vNodeIds, vNorms);

  PotentialField::BOUNDARY_TYPE nType = pBC != NULL ? pBC->nType : PotentialField::FIXED_VAL;
  pField->setBoundaryType(sName, nType);

  if(nType != PotentialField::FIXED_VAL)
    return;

  double fVal = 0;
  if((pBC != NULL) && (nType == PotentialField::FIXED_VAL))
    fVal = pBC->nFixedValType == CPotentialBoundCond::fvPlusUnity ? 1 : -1;
    
  pField->setBoundaryVal(sName, fVal);
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

V3D CElectricFieldData::calc_norm(CNode3D* pNode) const
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
  V3D v3d = { -vNorm.x, -vNorm.y, -vNorm.z };
  return v3d;
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

bool CElectricFieldData::get_result(PotentialField* pField) const
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CNodesCollection& vNodes = pObj->get_nodes();
  size_t nNodeCount = vNodes.size();
  const std::vector<double>& vRes = pField->getPotentialVals();
  if(vRes.size() != nNodeCount)
    return false;

  for(size_t i = 0; i < nNodeCount; i++)
  {
    CNode3D* pNode = vNodes.at(i);
    pNode->phi = (float)vRes.at(i);
  }

  return true;
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