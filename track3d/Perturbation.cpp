
#include "stdafx.h"

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
  m_vAddEdc = Vector3D(2000 * SI_to_CGS_ElecField, 0, 0); // by default this is 20 V/cm.
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

};  // namespace EvaporatingParticle