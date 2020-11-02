
#include "stdafx.h"

#include "ParticleTracking.h"
#include "Perturbation.h"
#include "BarnesHut.h"
#include "Primitives.h"
#include <math.h>

#include "../utilities/ParallelFor.h"

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
    case CFieldPerturbation::ptbFlatChannelRF: return new CAnalytRFField();
    case CFieldPerturbation::ptbCylSubstrateRF: return new CCurvedSubstrateRF();
    case CFieldPerturbation::ptbElliptSubstrRF: return new CEllipticalSubstrateRF();
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
    case CFieldPerturbation::ptbFlatChannelRF: return _T("Analytical RF Field, Flat Substrate");
    case CFieldPerturbation::ptbCylSubstrateRF: return _T("Analytical RF Field, Round Cylinder Substrate");
    case CFieldPerturbation::ptbElliptSubstrRF: return _T("Analytical RF Field, Elliptical Substrate");
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

bool CFieldPtbCollection::sel_region_changed(CStringVector* pRegNames)
{
  for(size_t i = 0; i < size(); i++)
  {
    if(at(i)->sel_region_changed(pRegNames))
      return true;
  }

  return false;
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
  CNodesVector& vNodes = pObj->get_nodes();
  return nNodeId < vNodes.size() ? field_ptb(vNodes.at(nNodeId).pos) : vNull;
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
  CNodesVector& vNodes = pObj->get_nodes();
  return nNodeId < vNodes.size() ? field_ptb(vNodes.at(nNodeId).pos) : vNull;
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
  set_charge_srf_dens(0.001 / Const_Srf_Charge_Dens); // 0.001 nA*hour per square millimeter.
  m_bEnableMultiThreading = true;
}

// This function is for backward compatibility only, it would be used in obsolete approach.
// The modern approach assumes preliminary field calculation at every mesh node.
Vector3D CDoubleLayerField::field_ptb(const Vector3D& vPos)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CRegionsCollection& regions = pObj->get_regions(false);
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
  CNodesVector& vNodes = pObj->get_nodes();
  size_t nNodesCount = vNodes.size();
  m_vPhi.resize(nNodesCount, 0);
  m_vField.resize(nNodesCount, Vector3F(vNull));

  const CRegionsCollection& regions = pObj->get_regions(false);
  size_t nRegCount = regions.size();
  if(nRegCount == 0)
    return false;

  if(m_bEnableMultiThreading)
  {
    ThreadPool::splitInPar(vNodes.size(),
	    [&](size_t i) 
    {
      Vector3D vPos = vNodes.at(i).pos;
      CFace* pFace = NULL;
      CRegion* pReg = NULL;
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
    },
	    static_cast<CObject*>(this));
  }
  else
  {
    CFace* pFace = NULL;
    CRegion* pReg = NULL;
    Vector3D vPos;
    for(size_t i = 0; i < nNodesCount; i++)
    {
      vPos = vNodes.at(i).pos;
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
  }

  m_bReady = !get_terminate_flag();
  return m_bReady;
}

void CDoubleLayerField::calc_field_from_face(CFace* pFace, const Vector3D& vPos, Vector3F& vField, float& fPhi) const
{
  double fS = pFace->square();
  if(fS < Const_Almost_Zero)
    return;

  Vector3D vN = -pFace->norm; // we need normal vector directed inside the domain.
  Vector3D vC = (pFace->p0->pos + pFace->p1->pos + pFace->p2->pos) * Const_One_Third;

  do_calc_field_from_face(vPos, vC, vN, fS, vField, fPhi);

// Symmetry correction. Note that vPos is position of a mesh node, i.e. it is always inside the domain:
  Vector3D vCr, vNr;
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  int nSymType = pObj->get_symmetry_type();
  switch(nSymType)
  {
    case CBarnesHut::symXYonly:
    {
      vCr = Vector3D(vC.x, vC.y, -vC.z);
      vNr = Vector3D(vN.x, vN.y, -vN.z);
      do_calc_field_from_face(vPos, vCr, vNr, fS, vField, fPhi);
      break;
    }
    case CBarnesHut::symXZonly:
    {
      vCr = Vector3D(vC.x, -vC.y, vC.z);
      vNr = Vector3D(vN.x, -vN.y, vN.z);
      do_calc_field_from_face(vPos, vCr, vNr, fS, vField, fPhi);
      break;
    }
    case CBarnesHut::symBoth:
    {
      vCr = Vector3D(vC.x, vC.y, -vC.z);
      vNr = Vector3D(vN.x, vN.y, -vN.z);
      do_calc_field_from_face(vPos, vCr, vNr, fS, vField, fPhi);

      vCr = Vector3D(vC.x, -vC.y, vC.z);
      vNr = Vector3D(vN.x, -vN.y, vN.z);
      do_calc_field_from_face(vPos, vCr, vNr, fS, vField, fPhi);

      vCr = Vector3D(vC.x, -vC.y, -vC.z);
      vNr = Vector3D(vN.x, -vN.y, -vN.z);
      do_calc_field_from_face(vPos, vCr, vNr, fS, vField, fPhi);
      break;
    }
  }
}

void CDoubleLayerField::do_calc_field_from_face(const Vector3D& vPos, const Vector3D& vC, const Vector3D& vN, double fS, Vector3F& vField, float& fPhi) const
{
  Vector3D vR = vPos - vC;

  double fR2 = vR.sqlength();
  double fR = sqrt(fR2);
  double fR3 = fR2 * fR;
  if(fR3 < Const_Almost_Zero)
    return;

  vR /= fR; // normalize the radius-vector.
  double fDot = vR & vN;
  if(fDot < -Const_Almost_Zero)
    return;  // nothing to add to the field.

  fPhi += fS * fDot / fR2;
  vField += Vector3F((3 * fDot * vR - vN) * fS / fR3);  // both vR and vN are normalized.
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

  m_bReady = true;
}

//-----------------------------------------------------------------------------
// CAnalytRFField - a model of the RF field in a narrow channel from a set of 
//                  narrow stripes (MEMS devices).
//-----------------------------------------------------------------------------
CAnalytRFField::CAnalytRFField()
  : CFieldPerturbation(), m_fStripeWidth(0), m_fGapWidth(0), m_pA(NULL), m_pB(NULL), m_pLmb(NULL), m_pRes(NULL), m_pBoundShape(NULL)
{
  m_nType = CFieldPerturbation::ptbFlatChannelRF;
  set_default();
}
 
CAnalytRFField::~CAnalytRFField()
{
  delete_arrays();
  if(m_pBoundShape != NULL)
    delete m_pBoundShape;
}

Vector3D CAnalytRFField::field_on_the_fly(const Vector3D& vPos, double fTime, double fPhase)
{
  if(!m_bReady)   // m_bEnable is checked before call (in CTracker::get_ion_accel()).
    return vNull;

  if(!m_pBoundShape->contains(vPos))
    return vNull;

  Vector2D v = world_to_model_coord(vPos);

  UINT ix, jy;
  double ksix, ksiy;
  if(!cell_indices(v, ix, ksix, jy, ksiy))
    return vNull;

  Vector2D r00 = m_pRes[ix][jy];
  bool bInsX = ix < m_nNx - 1;
  Vector2D r10 = bInsX ? m_pRes[ix + 1][jy] : m_pRes[0][jy];  // periodicity in x.
  bool bInsY = jy < m_nNy - 1;
  Vector2D r01 = bInsY ? m_pRes[ix][jy + 1] : Vector2D(0, 0); // zero field at large y.

  Vector2D r11;
  if(bInsX && bInsY)
    r11 = m_pRes[ix + 1][jy + 1];
  else if(bInsY)
    r11 = m_pRes[0][jy + 1];
  else
    r11 = Vector2D(0, 0);

  Vector2D vRes = r00 * (1 - ksix - ksiy + ksix*ksiy) + r10 * (1 - ksiy) * ksix + r01 * (1 - ksix) * ksiy + r11 * ksix * ksiy;

  Vector3D vLocX, vLocZ(0, 0, 1);
  Vector3D vLocY = m_pBoundShape->get_loc_normal(vPos);
  switch(m_nTransDir)
  {
    case strOrtFlowDir:
    {
      vLocX = vLocY * vLocZ;
      break;
    }
    case strAlongFlowDir:
    {
      vLocX = vLocZ;
      break;
    }
  }

  Vector3D vF = vRes.x * vLocX + vRes.y * vLocY;

  vF *= m_fAmpl * sin(m_fOmega * fTime + fPhase);
  return vF;
}

Vector2D CAnalytRFField::world_to_model_coord(const Vector3D& vPos) const
{
  Vector2D v = m_pBoundShape->model_coord(vPos);
  UINT n = UINT(v.x / m_fWaveLength);
  if(n > 0)
    v.x -= n * m_fWaveLength;

  return v;
}

bool CAnalytRFField::prepare()
{
  if(m_bReady)
    return true;

  if(!calc_bounding_box())
    return false;

  if(!calc_coeff())
  {
    AfxMessageBox(_T("CAnalytRFField::prepare(): Bad perturbation parameters. Simulation terminated."));
    CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
    pObj->terminate();
    return false;
  }

  Vector2D v;
  for(UINT i = 0; i < m_nNx; i++)
  {
    for(UINT j = 0; j < m_nNy; j++)
    {
      v = model_coord(i, j);
      m_pRes[i][j] = field_in_model_coord(v.x, v.y);
    }
  }

// DEBUG
  test_output();
// END DEBUG

  m_bReady = true;
  return true;
}

void CAnalytRFField::set_default()
{
  m_nVersion = 2;   // backward compatibility.

  m_bReady = false;
  m_bOnTheFly = true;

  m_bRegVisible = true;
  m_nMergeOpt = CSelectedAreas::optSubst;

  set_stripe_width(0.0024); // 24 mcm.
  set_gap_width(0.0016);    // 16 mcm.

  m_fChannelHeight = 0.01;  // 100 mcm.

  set_rf_freq(2.0e+7);                  // radio frequency, 20 MHz.
  m_fAmpl = 150.0 * SI_to_CGS_Voltage;  // amplitude, 150 V in CGSE.

  m_nDecompCount = 5;
  m_nTransDir = strOrtFlowDir;

// Dimensions of the result array:
  m_nNx = 30;
  m_nNy = 30;
}

void CAnalytRFField::allocate_arrays()
{
  delete_arrays();
  if(m_nDecompCount == 0)
    return;

  m_pA = new double[m_nDecompCount];
  m_pB = new double[m_nDecompCount];
  m_pLmb = new double[m_nDecompCount];

  m_pRes = new Vector2D*[m_nNx];
  for(UINT i = 0; i < m_nNx; i++)
    m_pRes[i] = new Vector2D[m_nNy];
}

void CAnalytRFField::delete_arrays()
{
  if(m_pA != NULL)
    delete[] m_pA;

  if(m_pB != NULL)
    delete[] m_pB;

  if(m_pLmb != NULL)
    delete[] m_pLmb;

  m_pA = NULL;
  m_pB = NULL;
  m_pLmb = NULL;

  if(m_pRes != NULL)
  {
    for(UINT i = 0; i < m_nNx; i++)
      delete[] m_pRes[i];

    delete[] m_pRes;
    m_pRes = NULL;
  }
}

bool CAnalytRFField::calc_coeff()
{
  if(m_fWaveLength < Const_Almost_Zero || m_nDecompCount == 0)
    return false;

  UINT k;
  allocate_arrays();
  for(UINT i = 0; i < m_nDecompCount; i++)
  {
    k = 2 * i + 1;

    m_pLmb[i] = Const_2PI * k / m_fWaveLength;

    double fArg = m_pLmb[i] * m_fChannelHeight;
    m_pB[i] = fArg < 10.0 ? -cosh(fArg) / sinh(fArg) : -1.0;

    m_pA[i] = coeff_Fourier(k);
  }

  return true;
}

double CAnalytRFField::coeff_Fourier(UINT k)
{
  double fCoeff = 4 * m_fWaveLength / (Const_PI * Const_PI * k * k * m_fGapWidth);
  double fArg = Const_PI * k * m_fGapWidth / m_fWaveLength;
  return fCoeff * sin(fArg);
}

Vector2D CAnalytRFField::field_in_model_coord(double x, double y) const
{
  Vector2D vRes(0, 0);
  if(x < 0 || x > m_fWaveLength || y < 0 || y > m_fChannelHeight) 
    return vRes;

  double argx, argy, sinx, cosx, sinhy, coshy, coeff;
  for(UINT i = 0; i < m_nDecompCount; i++)
  {
    argx = m_pLmb[i] * x;
    sinx = sin(argx);
    cosx = cos(argx);

    argy = m_pLmb[i] * y;
    sinhy = sinh(argy);
    coshy = cosh(argy);

    coeff = m_pA[i] * m_pLmb[i];
    vRes.x -= coeff * cosx * (coshy + m_pB[i] * sinhy);
    vRes.y -= coeff * sinx * (sinhy + m_pB[i] * coshy);
  }

  return vRes;
}

Vector2D CAnalytRFField::model_coord(UINT ix, UINT jy) const
{
  double x = ix * m_fWaveLength / m_nNx;          // periodicity in x, f(x) == f(x + m_fWaveLength).
  double y = jy * m_fChannelHeight / (m_nNy - 1); // at jy == m_nNy - 1, y must be equal to m_fChannelHeight.
  return Vector2D(x, y);
}

bool CAnalytRFField::cell_indices(const Vector2D& vModPos, UINT& ix, double& ksix, UINT& jy, double& ksiy) const
{
  if(vModPos.x < 0 || vModPos.y < 0)
    return false;

  double fStepX = m_fWaveLength / m_nNx;
  ix = UINT(vModPos.x / fStepX);
  if(ix >= m_nNx)
    return false;

  ksix = vModPos.x - ix * fStepX;

  double fStepY = m_fChannelHeight / m_nNy;
  jy = UINT(vModPos.y / fStepY);
  if(jy >= m_nNy)
    return false;

  ksiy = vModPos.y - jy * fStepY;
  return true;
}

bool CAnalytRFField::calc_bounding_box()
{
  if(m_pBoundShape != NULL)
    delete m_pBoundShape;

  m_pBoundShape = new CPrimitiveBox(m_vRegNames, m_fChannelHeight);
  return m_pBoundShape->is_ready();
}

const char* CAnalytRFField::get_trans_dir_name(int nDir) const
{
  switch(nDir)
  {
    case strOrtFlowDir: return _T(" Orthogonal to Flow");
    case strAlongFlowDir: return _T(" Parallel to Flow");
  }
  return _T(" ");
}

bool CAnalytRFField::sel_region_changed(CStringVector* pRegNames)
{
  if((DWORD_PTR)pRegNames == get_region_names_ptr())
  {
    m_bReady = false;
    return true;
  }

  return false;
}

void CAnalytRFField::save(CArchive& ar)
{
  UINT nVersion = 3;  // since 3 the CCurvedSubstrateRF::save() saves the x0, y0 and R to the stream.
  ar << nVersion;

  CFieldPerturbation::save(ar);

  ar << m_fStripeWidth;
  ar << m_fGapWidth;

  ar << m_fChannelHeight;

  ar << m_fFreq;
  ar << m_fAmpl;

  ar << m_nDecompCount;

  ar << m_nTransDir;

  size_t nRegNamesCount = m_vRegNames.size();
  ar << nRegNamesCount;
  for(size_t i = 0; i < nRegNamesCount; i++)
  {
    CString sName(m_vRegNames.at(i).c_str());
    ar << sName;
  }
}

void CAnalytRFField::load(CArchive& ar)
{
  ar >> m_nVersion;

  CFieldPerturbation::load(ar);

  ar >> m_fStripeWidth;
  ar >> m_fGapWidth;

  ar >> m_fChannelHeight;
  if(m_nVersion < 2)
  {
    double fChannelLength, fChannelWidth;
    ar >> fChannelLength;
    ar >> fChannelWidth;
  }

  ar >> m_fFreq;
  ar >> m_fAmpl;

  ar >> m_nDecompCount;

  if(m_nVersion < 2)
  {
    Vector3D vStartPoint;
    ar >> vStartPoint.x;
    ar >> vStartPoint.y;
    ar >> vStartPoint.z;
  }

  if(m_nVersion >= 1)
    ar >> m_nTransDir;

  if(m_nVersion >= 2)
  {
    size_t nRegNamesCount;
    ar >> nRegNamesCount;
    if(nRegNamesCount > 0)
      m_vRegNames.reserve(nRegNamesCount);

    for(size_t i = 0; i < nRegNamesCount; i++)
    {
      CString sName;
      ar >> sName;
      std::string stdName((const char*)sName);
      m_vRegNames.push_back(stdName);
    }
  }
}

bool CAnalytRFField::test_output()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  std::string cPath = COutputEngine::get_full_path(pObj->get_filename());
  std::string cName("Analyt_RF_Field");
  std::string cExt(".csv");
  std::string cFileName = cPath + cName + cExt;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return false;

  double fStepX = m_fWaveLength / m_nNx;
  double fStepY = m_fChannelHeight / (m_nNy - 1);
  float x = 0, y = 0, f = 0;
// For Origin we need only potential data, the X and Y values are assumed to be equidistant.
  for(UINT j = 0; j < m_nNy; j++)
  {
    y = j * fStepY;
    for(UINT i = 0; i <= m_nNx; i++)
    {
      x = i * fStepX;
      f = phi_in_model_coord(x, y);
      if(i == m_nNx)
        fprintf_s(pStream, "%7.4f\n", f);
      else
        fprintf_s(pStream, "%7.4f,", f);
    }
  }

  fclose(pStream);
  return true;
}

double CAnalytRFField::phi_in_model_coord(double x, double y) const
{
  double phi = 0, fArgX, fArgY;
  for(UINT i = 0; i < m_nDecompCount; i++)
  {
    fArgX = m_pLmb[i] * x;
    fArgY = m_pLmb[i] * y;
    phi += m_pA[i] * sin(fArgX) * (cosh(fArgY) + m_pB[i] * sinh(fArgY));
  }

  return phi;
}

//-----------------------------------------------------------------------------
// CCurvedSubstrateRF - a model of the RF field from a set of narrow stripes
//                     (MEMS devices). Round cylinder shape of the substrate.
//-----------------------------------------------------------------------------
CCurvedSubstrateRF::CCurvedSubstrateRF()
{
  m_nType = CFieldPerturbation::ptbCylSubstrateRF;

  set_default();

  m_fx0 = 0;
  m_fy0 = 3.004;
  m_fRad = 3.0;
}

bool CCurvedSubstrateRF::calc_bounding_box()
{
  if(m_pBoundShape != NULL)
    delete m_pBoundShape;

  m_pBoundShape = new CPrimCylSector(m_vRegNames, m_fx0, m_fy0, m_fRad, m_fChannelHeight);
  return m_pBoundShape->is_ready();
}

void CCurvedSubstrateRF::save(CArchive& ar)
{
  CAnalytRFField::save(ar);

  UINT nVersion = 0;
  ar << nVersion;

  ar << m_fx0;
  ar << m_fy0;
  ar << m_fRad;
}

void CCurvedSubstrateRF::load(CArchive& ar)
{
  CAnalytRFField::load(ar);
  UINT nParentVersion = get_version();
  if(nParentVersion < 3)
    return;   // backward compatibility;

  UINT nVersion;
  ar >> nVersion;

  ar >> m_fx0;
  ar >> m_fy0;
  ar >> m_fRad;
}

//-----------------------------------------------------------------------------
// CEllipticalSubstrateRF - Elliptical cylinder shape of the substrate.
//-----------------------------------------------------------------------------
CEllipticalSubstrateRF::CEllipticalSubstrateRF()
{
  m_nType = CFieldPerturbation::ptbElliptSubstrRF;

  set_default();

  m_fx0 = 0;
  m_fy0 = 1.304;
  m_fa = 2.0;
  m_fb = 1.3;
}

bool CEllipticalSubstrateRF::calc_bounding_box()
{
  if(m_pBoundShape != NULL)
    delete m_pBoundShape;

  CEllipseData ellipse(m_fx0, m_fy0, m_fa, m_fb);
  m_pBoundShape = new CEllipticalCylSector(m_vRegNames, ellipse, m_fChannelHeight, 5e-6);
  return m_pBoundShape->is_ready();
}

void CEllipticalSubstrateRF::save(CArchive& ar)
{
  CAnalytRFField::save(ar);

  UINT nVersion = 0;
  ar << nVersion;

  ar << m_fx0;
  ar << m_fy0;
  ar << m_fa;
  ar << m_fb;
}

void CEllipticalSubstrateRF::load(CArchive& ar)
{
  CAnalytRFField::load(ar);

  UINT nVersion;
  ar >> nVersion;

  ar >> m_fx0;
  ar >> m_fy0;
  ar >> m_fa;
  ar >> m_fb;
}

};  // namespace EvaporatingParticle