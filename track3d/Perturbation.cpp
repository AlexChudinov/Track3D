
#include "stdafx.h"

#include "ParticleTracking.h"
#include "Perturbation.h"
#include <math.h>

namespace EvaporatingParticle
{

//-----------------------------------------------------------------------------
// CFieldPerturbation - base class for various field perturbations
//-----------------------------------------------------------------------------
CFieldPerturbation* CFieldPerturbation::create(int nType)
{
  switch(nType)
  {
    case CFieldPerturbation::ptbRing: return new CChargedRingPerturbation();
    case CFieldPerturbation::ptbStackOfRings: return new CStackRingPerturbation();
    case CFieldPerturbation::ptbUniform: return new CUniformAddField();
    case CFieldPerturbation::ptbDoubleLayer: return new CDoubleLayerField();
  }

  return NULL;
}

const char* CFieldPerturbation::perturbation_name(int nType)
{
  switch(nType)
  {
    case CFieldPerturbation::ptbRing: return _T("Single Charged Dielectric Ring");
    case CFieldPerturbation::ptbStackOfRings: return _T("Stack of Charged Rings");
    case CFieldPerturbation::ptbUniform: return _T("Uniform Field");
    case CFieldPerturbation::ptbDoubleLayer: return _T("Thin Charged Dielectric Film");
  }

  return _T("None");
}

void CFieldPerturbation::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  ar << m_bEnable;
}

void CFieldPerturbation::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_bEnable;
}

void CFieldPerturbation::invalidate_contours() const
{
  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->invalidate_contours();
}

//-----------------------------------------------------------------------------
// CFieldPtbCollection
//-----------------------------------------------------------------------------
CFieldPtbCollection::~CFieldPtbCollection()
{
  clear_perturbations();
}

void CFieldPtbCollection::clear_perturbations()
{
  for(size_t i = 0; i < size(); i++)
  {
    CFieldPerturbation* pObj = at(i);
    delete pObj;
  }

  clear();
}

void CFieldPtbCollection::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  size_t nObjCount = size();
  ar << nObjCount;
  for(size_t i = 0; i < nObjCount; i++)
  {
    CFieldPerturbation* pObj = at(i);
    int nType = pObj->type();
    ar << nType;
    pObj->save(ar);
  }
}

void CFieldPtbCollection::load(CArchive& ar)
{
  clear_perturbations();

  UINT nVersion;
  ar >> nVersion;

  size_t nObjCount;
  ar >> nObjCount;
  for(size_t i = 0; i < nObjCount; i++)
  {
    int nType;
    ar >> nType;
    CFieldPerturbation* pObj = CFieldPerturbation::create(nType);
    pObj->load(ar);
    push_back(pObj);
  }
}

//-----------------------------------------------------------------------------
// CChargedRingPerturbation - a uniformly charged dielectric ring
//-----------------------------------------------------------------------------
CChargedRingPerturbation::CChargedRingPerturbation()
  : CFieldPerturbation()
{
  m_nType = CFieldPerturbation::ptbRing;
  set_default();
}

CChargedRingPerturbation::~CChargedRingPerturbation()
{
  delete[] m_pCosPhi;
}

void CChargedRingPerturbation::set_default()
{
  m_bReady = false;
  m_fRadius = 0.03;   // cm, the default radius is 0.3 mm.
  m_fCharge = 1.0e+6 * Const_Charge_CGS;  // 10^6 elementary charges.
  m_vCenter = Vector3D(2., 0., 0.);  // 20 mm.

  m_nDivCount = 30;
  m_pCosPhi = new double[m_nDivCount];
  double fdPhi = get_dPhi();
  for(UINT i = 0; i < m_nDivCount; i++)
    m_pCosPhi[i] = cos(i * fdPhi);
}

bool CChargedRingPerturbation::prepare()
{
  m_fCoeff = get_dPhi() * m_fCharge / (Const_PI * m_fRadius * m_fRadius);
  m_bReady = true;
  return true;
}

static const Vector3D vNull(0, 0, 0);

Vector3D CChargedRingPerturbation::field_ptb(const Vector3D& vPos)
{
  if(!m_bEnable || fabs(m_fCharge) < Const_Almost_Zero)
    return vNull;

  if(!m_bReady && !prepare())
    return vNull;

  Vector3D vRelPos = vPos - m_vCenter;
  double x = vRelPos.x;

  Vector3D vRad = Vector3D(0., vRelPos.y, vRelPos.z);
  double r = vRad.length();
  bool bAtAxis = r < Const_Almost_Zero;
  if(!bAtAxis)
    vRad /= r;  // unit vector in the radial direction.

  bool bAtPlane = fabs(x) < Const_Almost_Zero;
  if(bAtPlane && bAtAxis)
    return vNull;

  double u = x / m_fRadius;
  double v = r / m_fRadius;
  double a = 1 + u * u + v * v;

  double dEx, dEr;
  calc_dEx_dEr(0, a, v, dEx, dEr);

  double Ex = 0.5 * dEx;
  double Er = 0.5 * dEr;
  UINT nCount = m_nDivCount - 1;
  for(UINT i = 1; i < nCount; i++)
  {
    calc_dEx_dEr(i, a, v, dEx, dEr);
    Ex += dEx;
    Er += dEr;
  }

  calc_dEx_dEr(nCount, a, v, dEx, dEr);
  Ex += 0.5 * dEx;
  Er += 0.5 * dEr;

  Ex *= u * m_fCoeff;
  Er *= m_fCoeff;

  return Vector3D(Ex, vRad.y * Er, vRad.z * Er);
}

Vector3D CChargedRingPerturbation::get_field(size_t nNodeId)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CNodesCollection& vNodes = pObj->get_nodes();
  return nNodeId < vNodes.size() ? field_ptb(vNodes.at(nNodeId)->pos) : vNull;
}

float CChargedRingPerturbation::get_phi(size_t nNodeId)
{
  return 0;
}

static const double scfEps = 1e-9;

void CChargedRingPerturbation::calc_dEx_dEr(UINT i, double a, double v, double& dEx, double& dEr) const
{
  double D = a - 2 * v * m_pCosPhi[i];
  if(fabs(D) < scfEps)
    D = scfEps;

  D = 1. / sqrt(D * D * D);

  dEx = D;
  dEr = (v - m_pCosPhi[i]) * D;
}

void CChargedRingPerturbation::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CFieldPerturbation::save(ar);

  ar << m_vCenter.x;
  ar << m_vCenter.y;
  ar << m_vCenter.z;

  ar << m_fRadius;
  ar << m_fCharge;
}

void CChargedRingPerturbation::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CFieldPerturbation::load(ar);

  ar >> m_vCenter.x;
  ar >> m_vCenter.y;
  ar >> m_vCenter.z;

  ar >> m_fRadius;
  ar >> m_fCharge;
}

//-----------------------------------------------------------------------------
// CStackRingPerturbation - a stack of charged dielectric rings
//-----------------------------------------------------------------------------
CStackRingPerturbation::CStackRingPerturbation()
  : CFieldPerturbation(),
    m_bReady(false)
{
  m_nType = CFieldPerturbation::ptbStackOfRings;
  set_default();
}

CStackRingPerturbation::~CStackRingPerturbation()
{
  m_vChargedRings.clear_perturbations();
  m_vRingCharges.clear();
}

void CStackRingPerturbation::set_default()
{
  m_nChargeDistrType = distrUniform;
  m_nRingsCount = 10;
  m_fSumCharge = 1e+7 * Const_Charge_CGS;
  m_fRadius = 0.03;

  m_vStackBeg = Vector3D(0.2, 0., 0.);
  m_vStackEnd = Vector3D(4.2, 0., 0.);
}

Vector3D CStackRingPerturbation::field_ptb(const Vector3D& vPos)
{
  if(!m_bEnable || fabs(m_fSumCharge) < Const_Almost_Zero)
    return vNull;

  if(!m_bReady && !prepare())
    return vNull;

  return m_vChargedRings.apply(vPos);
}

Vector3D CStackRingPerturbation::get_field(size_t nNodeId)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CNodesCollection& vNodes = pObj->get_nodes();
  return nNodeId < vNodes.size() ? field_ptb(vNodes.at(nNodeId)->pos) : vNull;
}

float CStackRingPerturbation::get_phi(size_t nNodeId)
{
  return 0;
}

bool CStackRingPerturbation::prepare()
{
  m_vChargedRings.clear_perturbations();
  if(m_nRingsCount == 0 || m_vStackBeg == m_vStackEnd)
    return false;

  calc_ring_charges();

  double fKsi;
  Vector3D vPos;
  m_vChargedRings.reserve(m_nRingsCount);
  for(UINT i = 0; i < m_nRingsCount; i++)
  {
    fKsi = m_nRingsCount > 1 ? double(i) / (m_nRingsCount - 1) : 0;
    vPos = (1. - fKsi) * m_vStackBeg + fKsi * m_vStackEnd;
    CChargedRingPerturbation* pRingPtb = new CChargedRingPerturbation();
    pRingPtb->set_ring_pos(vPos);
    pRingPtb->set_ring_charge(m_vRingCharges.at(i));
    pRingPtb->set_ring_radius(m_fRadius);
    m_vChargedRings.push_back(pRingPtb);
  }

  m_bReady = true;
  return true;
}

void CStackRingPerturbation::calc_ring_charges()
{
  m_vRingCharges.clear();
  if(m_nRingsCount == 0)
    return;

  m_vRingCharges.reserve(m_nRingsCount);
  switch(m_nChargeDistrType)
  {
    case distrUniform:
    {
      double fQ = m_fSumCharge / m_nRingsCount;
      for(UINT i = 0; i < m_nRingsCount; i++)
        m_vRingCharges.push_back(fQ);
    }
    case distrRegress:
    {
// Temporarily! TO DO: write a regressive charge distribution, when q(x) ~ dI/dx.
      double fQ = m_fSumCharge / m_nRingsCount;
      for(UINT i = 0; i < m_nRingsCount; i++)
        m_vRingCharges.push_back(fQ);
    }
  }
}

const char* CStackRingPerturbation::get_distr_type_name(int nType) const
{
  switch(nType)
  {
    case distrUniform: return _T("Uniform");
    case distrRegress: return _T("Regressive");
  }

  return _T("None");
}

void CStackRingPerturbation::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CFieldPerturbation::save(ar);

  ar << m_nRingsCount;
  ar << m_fRadius;
  ar << m_fSumCharge;
  ar << m_nChargeDistrType;

  ar << m_vStackBeg.x;
  ar << m_vStackBeg.y;
  ar << m_vStackBeg.z;

  ar << m_vStackEnd.x;
  ar << m_vStackEnd.y;
  ar << m_vStackEnd.z;
}

void CStackRingPerturbation::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CFieldPerturbation::load(ar);

  ar >> m_nRingsCount;
  ar >> m_fRadius;
  ar >> m_fSumCharge;
  ar >> m_nChargeDistrType;

  ar >> m_vStackBeg.x;
  ar >> m_vStackBeg.y;
  ar >> m_vStackBeg.z;

  ar >> m_vStackEnd.x;
  ar >> m_vStackEnd.y;
  ar >> m_vStackEnd.z;
}

//-----------------------------------------------------------------------------
// CUniformAddField - a simple uniform addition to the external DC field
//-----------------------------------------------------------------------------
CUniformAddField::CUniformAddField()
  : CFieldPerturbation()
{
  m_nType = CFieldPerturbation::ptbUniform;
  set_default();
}

CUniformAddField::~CUniformAddField()
{
}

void CUniformAddField::set_default()
{
  m_vAddEdc = Vector3D(200 * SI_to_CGS_ElecField, 0, 0); // by default this is 2 V/cm.
  m_fAddEdcBegX = 5.85;   // cm.
  m_fAddEdcEndX = 11.29;  // cm.
}

Vector3D CUniformAddField::field_ptb(const Vector3D& vPos)
{
  return m_bEnable && (vPos.x >= m_fAddEdcBegX) && (vPos.x <= m_fAddEdcEndX) ? m_vAddEdc : vNull;
}

void CUniformAddField::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CFieldPerturbation::save(ar);

  ar << m_vAddEdc.x;
  ar << m_vAddEdc.y;
  ar << m_vAddEdc.z;

  ar << m_fAddEdcBegX;
  ar << m_fAddEdcEndX;
}

void CUniformAddField::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CFieldPerturbation::load(ar);

  ar >> m_vAddEdc.x;
  ar >> m_vAddEdc.y;
  ar >> m_vAddEdc.z;

  ar >> m_fAddEdcBegX;
  ar >> m_fAddEdcEndX;
}

//-----------------------------------------------------------------------------
// CDoubleLayerField - a model of the field from a charged thin dielectric film.
//-----------------------------------------------------------------------------
CDoubleLayerField::CDoubleLayerField()
{
  m_nType = CFieldPerturbation::ptbDoubleLayer;
  set_default();
}

CDoubleLayerField::~CDoubleLayerField()
{
  m_vField.clear();
  m_vPhi.clear();
}

void CDoubleLayerField::set_default()
{
  m_bReady = false;
  m_fFilmDepth = 1e-4;  // 1 mcm.
  set_charge_srf_dens(100 * Const_Srf_Charge_Dens); // 100 nA*hour per square millimeter.
}

// This function is for backward compatibility only, it would be used in obsolete approach.
// The modern approach assumes preliminary field calculation at every mesh node.
Vector3D CDoubleLayerField::field_ptb(const Vector3D& vPos)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CRegionsCollection& regions = pObj->get_regions();
  size_t nRegCount = regions.size();
  if(nRegCount == 0)
    return vNull;

  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  CFaceIndices vSelFaces = pDrawObj->get_sel_faces_ids();
  size_t nFacesCount = vSelFaces.size();
  if(nFacesCount == 0)
    return vNull;

  CFace* pFace = NULL;
  CRegion* pReg = NULL;
  Vector3F vRes(0, 0, 0);
  float fPhi = 0;
  for(size_t j = 0; j < nFacesCount; j++)
  {
    const CRegFacePair& face = vSelFaces.at(j);
    if(face.nReg >= nRegCount)
      continue;

    pReg = regions.at(face.nReg);
    if(face.nFace >= pReg->vFaces.size())
      continue;

    pFace = pReg->vFaces.at(face.nFace);
    calc_field_from_face(pFace, vPos, vRes, fPhi);
  }

  return vRes;
}

Vector3D CDoubleLayerField::get_field(size_t nNodeId)
{
  if(!m_bReady || !m_bEnable)
    return vNull;

  return nNodeId < m_vField.size() ? m_fScale * m_vField.at(nNodeId) : vNull;
}

float CDoubleLayerField::get_phi(size_t nNodeId)
{
  if(!m_bReady || !m_bEnable)
    return 0;

  return nNodeId < m_vPhi.size() ? CGS_to_SI_Voltage * m_fScale * m_vPhi.at(nNodeId) : 0;
}

bool CDoubleLayerField::calc_field()
{
  m_vPhi.clear();
  m_vField.clear();
  m_bReady = false;
  if(!m_bEnable)
    return false;

  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  CFaceIndices vSelFaces = pDrawObj->get_sel_faces_ids();
  size_t nFacesCount = vSelFaces.size();
  if(nFacesCount == 0)
  {
    AfxMessageBox("Field calculation failed: no faces were selected.");
    return false;
  }

  CString sJobName = CString(_T("Calculating field from ")) + CString(perturbation_name(ptbDoubleLayer)) + CString(_T("..."));
  set_job_name((const char*)sJobName);

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CNodesCollection& vNodes = pObj->get_nodes();
  size_t nNodesCount = vNodes.size();
  m_vPhi.resize(nNodesCount, 0);
  m_vField.resize(nNodesCount, Vector3F(vNull));

  const CRegionsCollection& regions = pObj->get_regions();
  size_t nRegCount = regions.size();
  if(nRegCount == 0)
    return false;

  CFace* pFace = NULL;
  CRegion* pReg = NULL;
  Vector3D vPos;
  for(size_t i = 0; i < nNodesCount; i++)
  {
    vPos = vNodes.at(i)->pos;
    if(i % 100 == 0)
      set_progress(100 * (i + 1) / nNodesCount);

    if(get_terminate_flag())
      return false;

    for(size_t j = 0; j < nFacesCount; j++)
    {
      const CRegFacePair& face = vSelFaces.at(j);
      if(face.nReg >= nRegCount)
        continue;

      pReg = regions.at(face.nReg);
      if(face.nFace >= pReg->vFaces.size())
        continue;

      pFace = pReg->vFaces.at(face.nFace);
      calc_field_from_face(pFace, vPos, m_vField[i], m_vPhi[i]);
    }
  }

  set_progress(100);
  m_bReady = true;
  return true;
}

bool CDoubleLayerField::calc_field_from_face(CFace* pFace, const Vector3D& vPos, Vector3F& vField, float& fPhi) const
{
  Vector3D vC = (pFace->p0->pos + pFace->p1->pos + pFace->p2->pos) * Const_One_Third;
  Vector3D vR = vPos - vC;

  double fR2 = vR.sqlength();
  double fR = sqrt(fR2);
  double fR3 = fR2 * fR;
  if(fR3 < Const_Almost_Zero)
    return false;

  vR /= fR; // normalize the radius-vector.
  Vector3D vN = -pFace->norm; // we need normal vector directed inside the domain.
  double fDot = vR & vN;
  if(fDot < -Const_Almost_Zero)
    return true;  // nothing to add to the field.

  Vector3D e1 = pFace->p1->pos - pFace->p0->pos;
  Vector3D e2 = pFace->p2->pos - pFace->p0->pos;
  double fS = 0.5 * (e1 * e2).length(); // square of the face.
  if(fS < Const_Almost_Zero)
    return false;

  fPhi += fS * fDot / fR2;
  vField += Vector3F((3 * fDot * vR - vN) * fS / fR3);  // both vR and vN are normalized.
  return true;
}

void CDoubleLayerField::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CFieldPerturbation::save(ar);

  ar << m_fFilmDepth;
  ar << m_fChargeSrfDens;

  size_t nNodesCount = m_bReady && (m_vPhi.size() == m_vField.size()) ? m_vPhi.size() : 0;
  ar << nNodesCount;
  for(size_t i = 0; i < nNodesCount; i++)
  {
    ar << m_vPhi[i];
    ar << m_vField[i].x;
    ar << m_vField[i].y;
    ar << m_vField[i].z;
  }
}

void CDoubleLayerField::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CFieldPerturbation::load(ar);

  ar >> m_fFilmDepth;
  ar >> m_fChargeSrfDens;

  m_fScale = 2 * m_fFilmDepth * m_fChargeSrfDens;

  m_vPhi.clear();
  m_vField.clear();
  size_t nNodesCount;
  ar >> nNodesCount;
  m_vPhi.resize(nNodesCount, 0.0f);
  m_vField.resize(nNodesCount, Vector3F(0, 0, 0));
  for(size_t i = 0; i < nNodesCount; i++)
  {
    ar >> m_vPhi[i];
    ar >> m_vField[i].x;
    ar >> m_vField[i].y;
    ar >> m_vField[i].z;
  }
}


};  // namespace EvaporatingParticle