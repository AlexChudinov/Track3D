
#pragma once

#include "Elements.h"

namespace EvaporatingParticle
{

//-----------------------------------------------------------------------------
// CFieldPerturbation - base class for various field perturbations
//-----------------------------------------------------------------------------
class CFieldPerturbation : public CObject
{
public:
  CFieldPerturbation()
    : m_bEnable(true)
  {
  }

  virtual ~CFieldPerturbation()
  {
  }

  enum
  {
    ptbRing         = 0,  // uniformly charged ring of finite radius and negligible width.
    ptbStackOfRings = 1,  // a stack of equidistant charged rings; the rings can carry different charges.
    ptbUniform      = 2,  // simple uniform addition dE = const to the external DC field.
    ptbDoubleLayer  = 3,  // a model of field produced by a thin charged dielectric film.
    ptbCount        = 4
  };

  static 
  CFieldPerturbation* create(int nType);

  static
  const char*         perturbation_name(int nType);

  virtual Vector3D    field_ptb(const Vector3D& vPos) = 0;

  virtual bool        calc_field() { return true; }
  virtual Vector3D    get_field(size_t nNodeId) { return Vector3D(); }
  virtual float       get_phi(size_t nNodeId) { return 0; }

  int                 type() const;
  const char*         name() const;

  bool                get_enable() const;
  DWORD_PTR           get_enable_ptr() const;

  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

protected:
  bool                m_bEnable;
  int                 m_nType;

  void                invalidate_contours() const;
};

//-----------------------------------------------------------------------------
// CFieldPtbCollection - collection of field perturbations
//-----------------------------------------------------------------------------
class CFieldPtbCollection : public std::vector<CFieldPerturbation*>
{
public:
  virtual ~CFieldPtbCollection();

  void                clear_perturbations();
  Vector3D            apply(const Vector3D& vPos) const;

  void                save(CArchive& ar);
  void                load(CArchive& ar);
};

//-----------------------------------------------------------------------------
// CChargedRingPerturbation - a uniformly charged dielectric ring
//-----------------------------------------------------------------------------
class CChargedRingPerturbation : public CFieldPerturbation
{
public:
  CChargedRingPerturbation();
  virtual ~CChargedRingPerturbation();

  virtual Vector3D    field_ptb(const Vector3D& vPos);

  virtual Vector3D    get_field(size_t nNodeId);
  virtual float       get_phi(size_t nNodeId);

  Vector3D            get_ring_pos() const;
  DWORD_PTR           get_ring_pos_ptr() const;
  void                set_ring_pos(const Vector3D& vC);

  double              get_ring_radius() const;
  DWORD_PTR           get_ring_radius_ptr() const;
  void                set_ring_radius(double fR);

  double              get_ring_charge() const;
  DWORD_PTR           get_ring_charge_ptr() const;
  void                set_ring_charge(double fCharge);  // CGSE.

  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

protected:
  void                set_default();
  bool                prepare();

  double              get_dPhi() const;

// u = x/R,  v = r/R,  a = 1 + u*u + v*v.
  void                calc_dEx_dEr(UINT i, double a, double v, double& dEx, double& dEr) const;

private:
  Vector3D            m_vCenter;

  double              m_fRadius,
                      m_fCharge;    // CGSE.

  double*             m_pCosPhi;
  UINT                m_nDivCount;

// Run-time:
  bool                m_bReady;
  double              m_fCoeff;
};

//-----------------------------------------------------------------------------
// CStackRingPerturbation - a stack of charged dielectric rings
//-----------------------------------------------------------------------------
class CStackRingPerturbation : public CFieldPerturbation
{
public:
  CStackRingPerturbation();
  virtual ~CStackRingPerturbation();

  enum  // charge distribution type.
  {
    distrUniform = 0,
    distrRegress = 1,
    distrCount   = 2
  };

  virtual Vector3D    field_ptb(const Vector3D& vPos);

  virtual Vector3D    get_field(size_t nNodeId);
  virtual float       get_phi(size_t nNodeId);

  UINT                get_rings_count() const;
  DWORD_PTR           get_rings_count_ptr() const;
  void                set_rings_count(UINT nCount);

  Vector3D            get_stack_beg_pos() const;
  DWORD_PTR           get_stack_beg_pos_ptr() const;
  void                set_stack_beg_pos(const Vector3D& vBeg);

  Vector3D            get_stack_end_pos() const;
  DWORD_PTR           get_stack_end_pos_ptr() const;
  void                set_stack_end_pos(const Vector3D& vEnd);

  double              get_ring_radius() const;
  DWORD_PTR           get_ring_radius_ptr() const;
  void                set_ring_radius(double fR);

  double              get_sum_charge() const;
  DWORD_PTR           get_sum_charge_ptr() const;
  void                set_sum_charge(double fCharge);  // CGSE.

  int                 get_charge_distr_type() const;
  DWORD_PTR           get_charge_distr_type_ptr() const;
  void                set_charge_distr_type(int nType);

  const char*         get_distr_type_name(int nType) const;

  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

protected:
  void                set_default();
  bool                prepare();

  void                calc_ring_charges();

private:
  CFieldPtbCollection m_vChargedRings;
  std::vector<double> m_vRingCharges;

  UINT                m_nRingsCount;

  double              m_fRadius,
                      m_fSumCharge;

  int                 m_nChargeDistrType;

  Vector3D            m_vStackBeg,
                      m_vStackEnd;
// Run-time:
  bool                m_bReady;
};

//-----------------------------------------------------------------------------
// CUniformAddField - a simple uniform addition to the external DC field
//-----------------------------------------------------------------------------
class CUniformAddField : public CFieldPerturbation
{
public:
  CUniformAddField();
  virtual ~CUniformAddField();

  virtual Vector3D    field_ptb(const Vector3D& vPos);

  virtual Vector3D    get_field(size_t nNodeId) { return m_vAddEdc; }
  virtual float       get_phi(size_t nNodeId)   { return 0; }

  Vector3D            get_add_Edc() const;
  DWORD_PTR           get_add_Edc_ptr() const;
  void                set_add_Edc(const Vector3D& vE);

  double              get_add_Edc_beg_x() const;
  DWORD_PTR           get_add_Edc_beg_x_ptr() const;
  void                set_add_Edc_beg_x(double fX);

  double              get_add_Edc_end_x() const;
  DWORD_PTR           get_add_Edc_end_x_ptr() const;
  void                set_add_Edc_end_x(double fX);

  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

protected:
  void                set_default();

private:
  Vector3D            m_vAddEdc;

  double              m_fAddEdcBegX,  // the additional DC field is applied at m_fAddEdcBegX < x < m_fAddEdcEndX.
                      m_fAddEdcEndX;
};

//-----------------------------------------------------------------------------
// CDoubleLayerField - a model of the field from a charged thin dielectric film.
//-----------------------------------------------------------------------------
class CDoubleLayerField : public CFieldPerturbation
{
public:
  CDoubleLayerField();
  virtual ~CDoubleLayerField();

  virtual Vector3D      field_ptb(const Vector3D& vPos);

  virtual bool          calc_field();
  virtual Vector3D      get_field(size_t nNodeId);
  virtual float         get_phi(size_t nNodeId);

  double                get_film_depth() const;
  DWORD_PTR             get_film_depth_ptr() const;
  void                  set_film_depth(double fD);

  double                get_charge_srf_dens() const;
  DWORD_PTR             get_charge_srf_dens_ptr() const;
  void                  set_charge_srf_dens(double fD);

  virtual void          save(CArchive& ar);
  virtual void          load(CArchive& ar);

protected:
  void                  set_default();

  bool                  calc_field_from_face(CFace* pFace, const Vector3D& vPos, Vector3F& vField, float& fPhi) const;

private:
  double                m_fFilmDepth,
                        m_fChargeSrfDens;

// Run-time:
  std::vector<float>    m_vPhi;   // for field visualization.
  std::vector<Vector3F> m_vField;

  float                 m_fScale;

  bool                  m_bReady;
};
//-----------------------------------------------------------------------------
// Inline implementation:
//-----------------------------------------------------------------------------
inline int CFieldPerturbation::type() const
{
  return m_nType;
}

inline const char* CFieldPerturbation::name() const
{
  return CFieldPerturbation::perturbation_name(m_nType);
}

inline bool CFieldPerturbation::get_enable() const
{
  return m_bEnable;
}

inline DWORD_PTR CFieldPerturbation::get_enable_ptr() const
{
  return (DWORD_PTR)&m_bEnable;
}

//-----------------------------------------------------------------------------
// CFieldPtbCollection
//-----------------------------------------------------------------------------
inline Vector3D CFieldPtbCollection::apply(const Vector3D& vPos) const
{
  Vector3D vRes;
  size_t nPtbCount = size();
  for(size_t i = 0; i < nPtbCount; i++)
    vRes += at(i)->field_ptb(vPos);

  return vRes;
}

//-----------------------------------------------------------------------------
// CChargedRingPerturbation
//-----------------------------------------------------------------------------
inline Vector3D CChargedRingPerturbation::get_ring_pos() const
{
  return m_vCenter;
}

inline DWORD_PTR CChargedRingPerturbation::get_ring_pos_ptr() const
{
  return (DWORD_PTR)&m_vCenter;
}

inline void CChargedRingPerturbation::set_ring_pos(const Vector3D& vC)
{
  if((m_vCenter - vC).length() > Const_Almost_Zero)
  {
    m_vCenter = vC;
    m_bReady = false;
  }
}

inline double CChargedRingPerturbation::get_ring_radius() const
{
  return m_fRadius;
}

inline DWORD_PTR CChargedRingPerturbation::get_ring_radius_ptr() const
{
  return (DWORD_PTR)&m_fRadius;
}

inline void CChargedRingPerturbation::set_ring_radius(double fR)
{
  if(m_fRadius != fR)
  {
    m_fRadius = fR;
    m_bReady = false;
  }
}

inline double CChargedRingPerturbation::get_ring_charge() const
{
  return m_fCharge;
}

inline DWORD_PTR CChargedRingPerturbation::get_ring_charge_ptr() const
{
  return (DWORD_PTR)&m_fCharge;
}

inline void CChargedRingPerturbation::set_ring_charge(double fCharge)
{
  if(m_fCharge != fCharge)
  {
    m_fCharge = fCharge;
    m_bReady = false;
  }
}

inline double CChargedRingPerturbation::get_dPhi() const
{
  return m_nDivCount > 1 ? Const_PI / (m_nDivCount - 1) : 0;
}

//-----------------------------------------------------------------------------
// CStackRingPerturbation - a stack of charged dielectric rings
//-----------------------------------------------------------------------------
inline UINT CStackRingPerturbation::get_rings_count() const
{
  return m_nRingsCount;
}

inline DWORD_PTR CStackRingPerturbation::get_rings_count_ptr() const
{
  return (DWORD_PTR)&m_nRingsCount;
}

inline void CStackRingPerturbation::set_rings_count(UINT nCount)
{
  if(m_nRingsCount != nCount)
  {
    m_nRingsCount = nCount;
    m_bReady = false;
  }
}

inline Vector3D CStackRingPerturbation::get_stack_beg_pos() const
{
  return m_vStackBeg;
}

inline DWORD_PTR CStackRingPerturbation::get_stack_beg_pos_ptr() const
{
  return (DWORD_PTR)&m_vStackBeg;
}

inline void CStackRingPerturbation::set_stack_beg_pos(const Vector3D& vBeg)
{
  if((m_vStackBeg - vBeg).length() > Const_Almost_Zero)
  {
    m_vStackBeg = vBeg;
    m_bReady = false;
  }
}

inline Vector3D CStackRingPerturbation::get_stack_end_pos() const
{
  return m_vStackEnd;
}

inline DWORD_PTR CStackRingPerturbation::get_stack_end_pos_ptr() const
{
  return (DWORD_PTR)&m_vStackEnd;
}

inline void CStackRingPerturbation::set_stack_end_pos(const Vector3D& vEnd)
{
  if((m_vStackEnd - vEnd).length() > Const_Almost_Zero)
  {
    m_vStackEnd = vEnd;
    m_bReady = false;
  }
}

inline double CStackRingPerturbation::get_ring_radius() const
{
  return m_fRadius;
}

inline DWORD_PTR CStackRingPerturbation::get_ring_radius_ptr() const
{
  return (DWORD_PTR)&m_fRadius;
}

inline void CStackRingPerturbation::set_ring_radius(double fR)
{
  if(m_fRadius != fR)
  {
    m_fRadius = fR;
    m_bReady = false;
  }
}

inline double CStackRingPerturbation::get_sum_charge() const
{
  return m_fSumCharge;
}

inline DWORD_PTR CStackRingPerturbation::get_sum_charge_ptr() const
{
  return (DWORD_PTR)&m_fSumCharge;
}

inline void CStackRingPerturbation::set_sum_charge(double fCharge)  // CGSE.
{
  if(m_fSumCharge != fCharge)
  {
    m_fSumCharge = fCharge;
    m_bReady = false;
  }
}

inline int CStackRingPerturbation::get_charge_distr_type() const
{
  return m_nChargeDistrType;
}

inline DWORD_PTR CStackRingPerturbation::get_charge_distr_type_ptr() const
{
  return (DWORD_PTR)&m_nChargeDistrType;
}

inline void CStackRingPerturbation::set_charge_distr_type(int nType)
{
  if(m_nChargeDistrType != nType)
  {
    m_nChargeDistrType = nType;
    m_bReady = false;
  }
}

//-----------------------------------------------------------------------------
// CUniformAddField
//-----------------------------------------------------------------------------
inline Vector3D CUniformAddField::get_add_Edc() const
{
  return m_vAddEdc;
}

inline DWORD_PTR CUniformAddField::get_add_Edc_ptr() const
{
  return (DWORD_PTR)&m_vAddEdc;
}

inline void CUniformAddField::set_add_Edc(const Vector3D& vE)
{
  m_vAddEdc = vE;
}

inline double CUniformAddField::get_add_Edc_beg_x() const
{
  return m_fAddEdcBegX;
}

inline DWORD_PTR CUniformAddField::get_add_Edc_beg_x_ptr() const
{
  return (DWORD_PTR)&m_fAddEdcBegX;
}

inline void CUniformAddField::set_add_Edc_beg_x(double fX)
{
  m_fAddEdcBegX = fX;
}

inline double CUniformAddField::get_add_Edc_end_x() const
{
  return m_fAddEdcEndX;
}

inline DWORD_PTR CUniformAddField::get_add_Edc_end_x_ptr() const
{
  return (DWORD_PTR)&m_fAddEdcEndX;
}

inline void CUniformAddField::set_add_Edc_end_x(double fX)
{
  m_fAddEdcEndX = fX;
}

//-----------------------------------------------------------------------------
// CDoubleLayerField - a model of the field from a charged thin dielectric film.
//-----------------------------------------------------------------------------
inline double CDoubleLayerField::get_film_depth() const
{
  return m_fFilmDepth;
}

inline DWORD_PTR CDoubleLayerField::get_film_depth_ptr() const
{
  return (DWORD_PTR)&m_fFilmDepth;
}

inline void CDoubleLayerField::set_film_depth(double fD)
{
  m_fFilmDepth = fD;
  m_fScale = 2 * m_fFilmDepth * m_fChargeSrfDens;
  invalidate_contours();
}

inline double CDoubleLayerField::get_charge_srf_dens() const
{
  return m_fChargeSrfDens;
}

inline DWORD_PTR CDoubleLayerField::get_charge_srf_dens_ptr() const
{
  return (DWORD_PTR)&m_fChargeSrfDens;
}

inline void CDoubleLayerField::set_charge_srf_dens(double fD)
{
  m_fChargeSrfDens = fD;
  m_fScale = 2 * m_fFilmDepth * m_fChargeSrfDens;
  invalidate_contours();
}


};  // namespace EvaporatingParticle
