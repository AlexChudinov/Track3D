
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
    : m_bEnable(true), m_bOnTheFly(false)
  {
  }

  virtual ~CFieldPerturbation()
  {
  }

  enum
  {
    ptbRing           = 0,  // uniformly charged ring of finite radius and negligible width.
    ptbStackOfRings   = 1,  // a stack of equidistant charged rings; the rings can carry different charges.
    ptbUniform        = 2,  // simple uniform addition dE = const to the external DC field.
    ptbDoubleLayer    = 3,  // a model of field produced by a thin charged dielectric film.
    ptbFlatChannelRF  = 4,  // a 2D analytic model of RF field caused by a set of metallic stripes (MEMS), flat substrate.
    ptbCylSubstrateRF = 5,  // the same as ptbFlatChannelRF but for curved (round cylinder) substrate.
    ptbElliptSubstrRF = 6,  // the same as ptbFlatChannelRF but for curved (elliptical) substrate.
    ptbCount          = 7
  };

  static 
  CFieldPerturbation* create(int nType);

  static
  const char*         perturbation_name(int nType);

  virtual bool        prepare() { return true; }

  virtual Vector3D    field_ptb(const Vector3D& vPos) = 0;
  virtual Vector3D    field_on_the_fly(const Vector3D& vPos, double fTime, double fPhase) { return Vector3D(); }

  virtual bool        calc_field() { return true; }
  virtual Vector3D    get_field(size_t nNodeId) { return Vector3D(); }
  virtual float       get_phi(size_t nNodeId) { return 0; }

  virtual bool        sel_region_changed(CStringVector* pRegNames) { return false; }

  int                 type() const;
  const char*         name() const;

  bool                get_enable() const;
  DWORD_PTR           get_enable_ptr() const;

  bool                get_on_the_fly() const { return m_bOnTheFly; }

  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

protected:
  bool                m_bEnable;
  bool                m_bOnTheFly;  // never changed, true for those perturbations, which need to be applied during ion tracking.
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

  bool                sel_region_changed(CStringVector* pRegNames);

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

  virtual bool        prepare();

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

  virtual bool        prepare();

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

  bool                  get_enable_multi_thread() const;
  DWORD_PTR             get_enable_multi_thread_ptr() const;

  virtual void          save(CArchive& ar);
  virtual void          load(CArchive& ar);

protected:
  void                  set_default();

  void                  calc_field_from_face(CFace* pFace, const Vector3D& vPos, Vector3F& vField, float& fPhi) const;
  void                  do_calc_field_from_face(const Vector3D& vPos, const Vector3D& vC, const Vector3D& vN, double fS, Vector3F& vField, float& fPhi) const;

private:
  double                m_fFilmDepth,
                        m_fChargeSrfDens;

  bool                  m_bEnableMultiThreading;

// Run-time:
  std::vector<float>    m_vPhi;   // for field visualization.
  std::vector<Vector3F> m_vField;

  float                 m_fScale;

  bool                  m_bReady;
};

class CPrimitive;
//-----------------------------------------------------------------------------
// CAnalytRFField - a base class for all models of the RF field in a narrow 
//                  channel from a set of narrow stripes (MEMS devices).
//-----------------------------------------------------------------------------
class CAnalytRFField : public CFieldPerturbation
{
public:
  CAnalytRFField();
  virtual ~CAnalytRFField();

  virtual bool        prepare();
  virtual Vector3D    field_ptb(const Vector3D& vPos) { return Vector3D(); }
  virtual Vector3D    field_on_the_fly(const Vector3D& vPos, double fTime, double fPhase);

  virtual bool        sel_region_changed(CStringVector* pRegNames);

// User interface:
  CStringVector&      get_region_names();
  DWORD_PTR           get_region_names_ptr() const;

  bool                get_visibility_flag() const;
  DWORD_PTR           get_visibility_flag_ptr() const;

// The existing regions can be merged with those in selected Named Areas. Three merge options are available: add, substitute and remove.
  int                 get_merge_option() const;
  DWORD_PTR           get_merge_option_ptr() const;
  void                set_merge_option(int nOpt);

  CString             get_last_merged() const;
  void                set_last_merged(const CString& sName);

  double              get_stripe_width() const;     // by default, 24 mcm in CGS.
  DWORD_PTR           get_stripe_width_ptr() const;
  void                set_stripe_width(double fW);

  double              get_gap_width() const;        // by default, 16 mcm in CGS.
  DWORD_PTR           get_gap_width_ptr() const;
  void                set_gap_width(double fW);

  double              get_channel_height() const;
  DWORD_PTR           get_channel_height_ptr() const;
  void                set_channel_height(double fH);

  double              get_rf_freq() const;
  DWORD_PTR           get_rf_freq_ptr() const;
  void                set_rf_freq(double fFreq);

  double              get_rf_ampl() const;
  DWORD_PTR           get_rf_ampl_ptr() const;
  void                set_rf_ampl(double fAmpl);

  UINT                get_decomp_count() const;
  DWORD_PTR           get_decomp_count_ptr() const;
  void                set_decomp_count(UINT nCount);

  enum // Stripes orientation. The virtual stripes can be stretched either orthogonal to the flow or parallel to the flow.
  { 
    strOrtFlowDir   = 0,
    strAlongFlowDir = 1,
    strCount        = 2
  };

  const char*         get_trans_dir_name(int nDir) const;

  int                 get_trans_dir() const;
  DWORD_PTR           get_trans_dir_ptr() const;
  void                set_trans_dir(int nTransDir);

  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

  CPrimitive*         get_bound_shape() const;

protected:
  void                set_default();

  void                allocate_arrays();
  void                delete_arrays();

  bool                calc_coeff();

  virtual bool        calc_bounding_box();

  double              coeff_Fourier(UINT k);

// Transformation from world to model coordinates (x, y), such that 0 < x < L; 0 < y < H. 
  Vector2D            world_to_model_coord(const Vector3D& vPos) const;

  double              phi_in_model_coord(double x, double y) const;
  Vector2D            field_in_model_coord(double x, double y) const;

  Vector2D            model_coord(UINT ix, UINT jy) const;
  bool                cell_indices(const Vector2D& vModPos, UINT& ix, double& ksix, UINT& jy, double& ksiy) const;

// DEBUG
  bool                test_output();
// END DEBUG

  UINT                get_version() const;

// Data:
// Bounding primitive. In the case of the flat substrate (this class) it is CPrimitiveBox; in the case of the
// cylinder (curved) substrate it is CPrimCylSector.
  CPrimitive*         m_pBoundShape;
// A set of region names belonging to one and the same surface - substrate.
  CStringVector       m_vRegNames;
  bool                m_bRegVisible;

  int                 m_nMergeOpt;
  CString             m_sLastMerged;

  double              m_fStripeWidth,   // structure parameters, stripe width and...
                      m_fGapWidth,      // ... gap between stripes.

                      m_fChannelHeight;

  double              m_fFreq,  // radio frequency, Hz.
                      m_fAmpl;  // amplitude, CGSE.

  int                 m_nTransDir;    // stripes orientation, can be either strOrtFlowDir or strAlongFlowDir.

  double              m_fOmega,       // run-time, m_fOmega = 2*pi*m_fFreq.
                      m_fWaveLength;  // run-time, m_fWaveLength = 2 * (m_fStripeWidth + m_fGapWidth).

  UINT                m_nDecompCount; // count of items in the Fourier decomposition.

  double*             m_pA;   // coefficients of Fourier decomposition of the boundary potential.
  double*             m_pB;   // coefficients at sinh(lambda*y) of U(y) = cosh(lambda*y) + B*sinh(lambda*y) decomposition.
  double*             m_pLmb; // m_pLmb[k] = Const_2PI*k/m_fWaveLength.

  UINT                m_nNx,  // dimensions of the result array m_pRes.
                      m_nNy;

  Vector2D**          m_pRes;  // pre-computed RF field's amplitude for a single cell (x = [0, L], y = [0, H]).

  bool                m_bReady; // run-time, false by default.

private:
  UINT                m_nVersion; // backward compatibility.
};

//-----------------------------------------------------------------------------
// CCurvedSubstrateRF
//-----------------------------------------------------------------------------
class CCurvedSubstrateRF : public CAnalytRFField
{
public:
  CCurvedSubstrateRF();

  double              get_x0() const;
  DWORD_PTR           get_x0_ptr() const;
  void                set_x0(double fX0);

  double              get_y0() const;
  DWORD_PTR           get_y0_ptr() const;
  void                set_y0(double fY0);

  double              get_radius() const;
  DWORD_PTR           get_radius_ptr() const;
  void                set_radius(double fR);

  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

protected:
  virtual bool        calc_bounding_box();

private:
  double              m_fx0,    // cylinder axix position.
                      m_fy0,
                      m_fRad;   // cylinder radius.
};

//-----------------------------------------------------------------------------
// CEllipticalSubstrateRF
//-----------------------------------------------------------------------------
class CEllipticalSubstrateRF : public CAnalytRFField
{
public:
  CEllipticalSubstrateRF();

  double              get_x0() const;
  DWORD_PTR           get_x0_ptr() const;
  void                set_x0(double fX0);

  double              get_y0() const;
  DWORD_PTR           get_y0_ptr() const;
  void                set_y0(double fY0);

  double              get_a() const;
  DWORD_PTR           get_a_ptr() const;
  void                set_a(double fA);

  double              get_b() const;
  DWORD_PTR           get_b_ptr() const;
  void                set_b(double fB);

  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

protected:
  virtual bool        calc_bounding_box();

private:
// Ellipse parameters:
  double              m_fx0,
                      m_fy0,
                      m_fa,
                      m_fb;
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

inline bool CDoubleLayerField::get_enable_multi_thread() const
{
  return m_bEnableMultiThreading;
}

inline DWORD_PTR CDoubleLayerField::get_enable_multi_thread_ptr() const
{
  return (DWORD_PTR)&m_bEnableMultiThreading;
}

//-----------------------------------------------------------------------------
// CAnalytRFField - inline implementation.
//-----------------------------------------------------------------------------
inline CStringVector& CAnalytRFField::get_region_names()
{
  return m_vRegNames;
}

inline DWORD_PTR CAnalytRFField::get_region_names_ptr() const
{
  return (DWORD_PTR)&m_vRegNames;
}

inline bool CAnalytRFField::get_visibility_flag() const
{
  return m_bRegVisible;
}

inline DWORD_PTR CAnalytRFField::get_visibility_flag_ptr() const
{
  return (DWORD_PTR)&m_bRegVisible;
}

inline int CAnalytRFField::get_merge_option() const
{
  return m_nMergeOpt;
}

inline DWORD_PTR CAnalytRFField::get_merge_option_ptr() const
{
  return (DWORD_PTR)&m_nMergeOpt;
}

inline void CAnalytRFField::set_merge_option(int nOpt)
{
  m_nMergeOpt = nOpt;
}

inline CString CAnalytRFField::get_last_merged() const
{
  return m_sLastMerged;
}

inline void CAnalytRFField::set_last_merged(const CString& sName)
{
  m_sLastMerged = sName;
}

inline double CAnalytRFField::get_gap_width() const
{
  return m_fGapWidth;
}

inline DWORD_PTR CAnalytRFField::get_gap_width_ptr() const
{
  return (DWORD_PTR)&m_fGapWidth;
}

inline void CAnalytRFField::set_gap_width(double fW)
{
  if(m_fGapWidth != fW)
  {
    m_fGapWidth = fW;
    m_fWaveLength = 2 * (m_fStripeWidth + m_fGapWidth);
    m_bReady = false;
  }
}

inline double CAnalytRFField::get_stripe_width() const
{
  return m_fStripeWidth;
}

inline DWORD_PTR CAnalytRFField::get_stripe_width_ptr() const
{
  return (DWORD_PTR)&m_fStripeWidth;
}

inline void CAnalytRFField::set_stripe_width(double fW)
{
  if(m_fStripeWidth != fW)
  {
    m_fStripeWidth = fW;
    m_fWaveLength = 2 * (m_fStripeWidth + m_fGapWidth);
    m_bReady = false;
  }
}

inline double CAnalytRFField::get_channel_height() const
{
  return m_fChannelHeight;
}

inline DWORD_PTR CAnalytRFField::get_channel_height_ptr() const
{
  return (DWORD_PTR)&m_fChannelHeight;
}

inline void CAnalytRFField::set_channel_height(double fH)
{
  if(m_fChannelHeight != fH)
  {
    m_fChannelHeight = fH;
    m_bReady = false;
  }
}

inline double CAnalytRFField::get_rf_freq() const
{
  return m_fFreq;
}

inline DWORD_PTR CAnalytRFField::get_rf_freq_ptr() const
{
  return (DWORD_PTR)&m_fFreq;
}

inline void CAnalytRFField::set_rf_freq(double fFreq)
{
  m_fFreq = fFreq;
  m_fOmega = Const_2PI * m_fFreq;
}

inline double CAnalytRFField::get_rf_ampl() const
{
  return m_fAmpl;
}

inline DWORD_PTR CAnalytRFField::get_rf_ampl_ptr() const
{
  return (DWORD_PTR)&m_fAmpl;
}

inline void CAnalytRFField::set_rf_ampl(double fAmpl)
{
  m_fAmpl = fAmpl;
}

inline UINT CAnalytRFField::get_decomp_count() const
{
  return m_nDecompCount;
}

inline DWORD_PTR CAnalytRFField::get_decomp_count_ptr() const
{
  return (DWORD_PTR)&m_nDecompCount;
}

inline void CAnalytRFField::set_decomp_count(UINT nCount)
{
  if(m_nDecompCount != nCount)
  {
    m_nDecompCount = nCount;
    m_bReady = false;
  }
}

inline int CAnalytRFField::get_trans_dir() const
{
  return m_nTransDir;
}

inline DWORD_PTR CAnalytRFField::get_trans_dir_ptr() const
{
  return (DWORD_PTR)&m_nTransDir;
}

inline void CAnalytRFField::set_trans_dir(int nTransDir)
{
  if(m_nTransDir != nTransDir)
  {
    m_nTransDir = nTransDir;
    m_bReady = false;
  }
}

inline CPrimitive* CAnalytRFField::get_bound_shape() const
{
  return m_pBoundShape;
}

inline UINT CAnalytRFField::get_version() const
{
  return m_nVersion;
}

//-----------------------------------------------------------------------------
// CCurvedSubstrateRF - in fact, a round cylinder substrate.
//-----------------------------------------------------------------------------
inline double CCurvedSubstrateRF::get_x0() const
{
  return m_fx0;
}

inline DWORD_PTR CCurvedSubstrateRF::get_x0_ptr() const
{
  return (DWORD_PTR)&m_fx0;
}

inline void CCurvedSubstrateRF::set_x0(double fX0)
{
  if(m_fx0 != fX0)
  {
    m_fx0 = fX0;
    m_bReady = false;
  }
}

inline double CCurvedSubstrateRF::get_y0() const
{
  return m_fy0;
}

inline DWORD_PTR CCurvedSubstrateRF::get_y0_ptr() const
{
  return (DWORD_PTR)&m_fy0;
}

inline void CCurvedSubstrateRF::set_y0(double fY0)
{
  if(m_fy0 != fY0)
  {
    m_fy0 = fY0;
    m_bReady = false;
  }
}

inline double CCurvedSubstrateRF::get_radius() const
{
  return m_fRad;
}

inline DWORD_PTR CCurvedSubstrateRF::get_radius_ptr() const
{
  return (DWORD_PTR)&m_fRad;
}

inline void CCurvedSubstrateRF::set_radius(double fR)
{
  if(m_fRad != fR)
  {
    m_fRad = fR;
    m_bReady = false;
  }
}

//-----------------------------------------------------------------------------
// CEllipticalSubstrateRF
//-----------------------------------------------------------------------------
inline double CEllipticalSubstrateRF::get_x0() const
{
  return m_fx0;
}

inline DWORD_PTR CEllipticalSubstrateRF::get_x0_ptr() const
{
  return (DWORD_PTR)&m_fx0;
}

inline void CEllipticalSubstrateRF::set_x0(double fX0)
{
  if(m_fx0 != fX0)
  {
    m_fx0 = fX0;
    m_bReady = false;
  }
}

inline double CEllipticalSubstrateRF::get_y0() const
{
  return m_fy0;
}

inline DWORD_PTR CEllipticalSubstrateRF::get_y0_ptr() const
{
  return (DWORD_PTR)&m_fy0;
}

inline void CEllipticalSubstrateRF::set_y0(double fY0)
{
  if(m_fy0 != fY0)
  {
    m_fy0 = fY0;
    m_bReady = false;
  }
}

inline double CEllipticalSubstrateRF::get_a() const
{
  return m_fa;
}

inline DWORD_PTR CEllipticalSubstrateRF::get_a_ptr() const
{
  return (DWORD_PTR)&m_fa;
}

inline void CEllipticalSubstrateRF::set_a(double fA)
{
  if(m_fa != fA)
  {
    m_fa = fA;
    m_bReady = false;
  }
}

inline double CEllipticalSubstrateRF::get_b() const
{
  return m_fb;
}

inline DWORD_PTR CEllipticalSubstrateRF::get_b_ptr() const
{
  return (DWORD_PTR)&m_fb;
}

inline void CEllipticalSubstrateRF::set_b(double fB)
{
  if(m_fb != fB)
  {
    m_fb = fB;
    m_bReady = false;
  }
}


};  // namespace EvaporatingParticle
