
#pragma once

#include "CObject.h"
#include "constant.hpp"
#include "../field_solver/MeshData.h"

#include <vector>
#include <string>

namespace EvaporatingParticle
{

//---------------------------------------------------------------------------------------
// CPotentialBoundCond
//---------------------------------------------------------------------------------------
struct CPotentialBoundCond
{
  CPotentialBoundCond(BoundaryMesh::BoundaryType type = BoundaryMesh::FIXED_VAL, int val = fvPlusUnity);

  enum  // fixed value types:
  {
    fvPlusUnity   = 0,
    fvMinusUnity  = 1,
    fvStepLike    = 2,
    fvCoulomb     = 3,    // attempt to calculate a mirror Coulomb field.
    fvCount       = 4
  };

  BoundaryMesh::BoundaryType    nType;

  CStringVector                 vRegNames;      // names of the regions with non-trivial boundary conditions.
  bool                          bVisible;       // visibility status of the selected regions.

  int                           nFixedValType;  // the value at the boundary can be +1V, -1V, step-wise potential or Coulomb potential.
  std::string                   sName;

// Step-like potential:
  double                        fStartX,
                                fStepX,
                                fEndX;

  static const char*            get_bc_type_name(BoundaryMesh::BoundaryType nType);
  static const char*            get_fixed_value_name(int nType);

  void                          save(CArchive& ar);
  void                          load(CArchive& ar);
};

struct CNode3D;
struct CRegion;
class CFiniteVolumesSolver;
typedef std::vector<CPotentialBoundCond*> CPotentialBoundCondColl;
//---------------------------------------------------------------------------------------
// CElectricFieldData - an auxiliary class for simulating scalable electric fields.
//---------------------------------------------------------------------------------------
class CElectricFieldData : public CObject
{
public:

  enum
  {
    cmLaplacian3    = 0,
    cmDirTessLap3   = 1,
    cmFinVolJacobi  = 2,
    cmCount         = 3
  };

  enum
  {
    typeFieldDC = 0,
    typeFieldRF = 1,
    typeMirror  = 2,    // attempt to calculate mirror Coulomb field.
    typeCount   = 3
  };

  CElectricFieldData(int nType = typeFieldDC);
  virtual ~CElectricFieldData();

  bool                    get_enable_field() const;
  DWORD_PTR               get_enable_field_ptr() const;

  int                     get_calc_method() const;
  DWORD_PTR               get_calc_method_ptr() const;
  void                    set_calc_method(int nCalcMethod);

  int                     get_type() const;
  DWORD_PTR               get_type_ptr() const;
  void                    set_type(int nType);

  double                  get_scale() const;        // use this function in the UI only, as it returns the scale in volts.
  DWORD_PTR               get_scale_ptr() const;
  void                    set_scale(double fScale);

  double                  get_freq() const;
  DWORD_PTR               get_freq_ptr() const;
  void                    set_freq(double fFreq);

  UINT                    get_iter_count() const;
  DWORD_PTR               get_iter_count_ptr() const;
  void                    set_iter_count(UINT nCount);

  bool                    get_analyt_field() const;
  DWORD_PTR               get_analyt_field_ptr() const;

  double                  get_inscr_radius() const;
  DWORD_PTR               get_inscr_radius_ptr() const;
  void                    set_inscr_radius(double fRadius);

  double                  get_low_analyt_lim() const;
  DWORD_PTR               get_low_analyt_lim_ptr() const;
  void                    set_low_analyt_lim(double fLowLimX);

  CString                 get_field_name() const;
  DWORD_PTR               get_field_name_ptr() const;
  void                    set_field_name(CString sName);

  size_t                  get_bc_count() const;
  CPotentialBoundCond*    get_bc(size_t nId) const;

  Vector3D                get_field(size_t nInd) const;
  double                  get_omega() const;
  double                  get_ampl() const;   // returns the scale in CGS units.

  void                    add_bc();
  void                    remove_bc(size_t nId);

  static const char*      get_field_type_name(int nType);
  static const char*      get_calc_method_name(int nCalcMethod);

  bool                    calc_field(bool bTest = false);
  bool                    need_recalc() const;
  void                    invalidate();

  void                    save(CArchive& ar);
  void                    load(CArchive& ar);

protected:
  void                    set_default();
  void                    clear_bc();

  CRegion*                get_region(const std::string& sName) const; // returns a region pointer by its name or NULL if the name is not found.
  bool                    is_selected(CRegion* pReg) const;           // returns true if the region is selected for boundary conditions.
  CIndexVector            get_reg_nodes(CRegion* pReg) const;         // returns the vector of global indices of the region nodes.

  bool                    set_boundary_conditions(CMeshAdapter& mesh);
  bool                    set_boundary_conditions(CFiniteVolumesSolver& solver);

  bool                    set_default_boundary_conditions(CMeshAdapter& mesh);
  bool                    set_default_boundary_conditions(CFiniteVolumesSolver& solver);

  void                    set_boundary_values(CMeshAdapter& mesh, CRegion* pReg, CPotentialBoundCond* pBC = NULL);

  double                  step_potential(CPotentialBoundCond* pBC, const Vector3D& vPos) const; // step-wise boundary conditions support.
  bool                    coulomb_potential(const CIndexVector& vNodeIds, std::vector<float>& vPhi) const;  // Coulomb boundary conditions support.

  Vector3D                calc_norm(CNode3D* pNode) const;

  bool                    get_result(bool bTest) const;
  void                    notify_scene(); // let the scene objects know that the potential field has been changed.

  void                    apply_analytic_field(const Vector3D& vPos, Vector3F& vField);

  bool                    calc_lap3(bool bTest);
  bool                    calc_dirichlet_lap3(bool bTest);
  bool                    calc_finite_vol_jacobi(bool bTest);

private:
  bool                    m_bEnable;
  int                     m_nCalcMethod,
                          m_nType;

  double                  m_fScale,       // potential scale, in V.
                          m_fOmega;       // circular frequency, for radio-frequency fields only.

  UINT                    m_nIterCount;

  CPotentialBoundCondColl m_vBoundCond;

// An attempt to get analytic field in the flatapole. Alpha version.
  bool                    m_bAnalytField;

  double                  m_fRadius,    // inscribed radius of the flatapole electrodes.
                          m_fLowLimX,   // an analytic formula will be used if m_fLowLimX < x < m_fHighLimX;
                          m_fHighLimX;

  std::vector<Vector3F>   m_vField;
  CString                 m_sName;

// Run-time:
  bool                    m_bNeedRecalc;
};

//-------------------------------------------------------------------------------------------------
// CFieldDataCollection.
//-------------------------------------------------------------------------------------------------
class CFieldDataColl : public std::vector<CElectricFieldData*>
{
public:
  CFieldDataColl();
  virtual ~CFieldDataColl();

// CFieldDataColl::calc_fields(true) is called only from CTracker::create_BH_object(...).
  bool              calc_fields(bool bMirrorClmb = false);
  bool              need_recalc() const;

  void              clear_fields_in_nodes();
  void              clear_fields();

  bool              sel_region_changed(CStringVector* pRegNames);
  bool              remove_bound_cond(CPotentialBoundCond* pBC);

// Visibility of regions support. 
  void              update_visibility_status(); // set visibility flag to all regions selected for boundary conditions.

  void              save(CArchive& ar);
  void              load(CArchive& ar);

  int               get_curr_field_index() const;
  void              set_curr_field_index(int nInd);

private:
  int               m_nCurrField;   // the field which properties can be currently edited by user.
};

inline bool CFieldDataColl::need_recalc() const
{
  for(size_t i = 0; i < size(); i++)
    if(at(i)->need_recalc())
      return true;

  return false;
}

inline int CFieldDataColl::get_curr_field_index() const
{
  return m_nCurrField < size() ? m_nCurrField : -1;
}

inline void CFieldDataColl::set_curr_field_index(int nInd)
{
  m_nCurrField = nInd;
}

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
inline bool CElectricFieldData::need_recalc() const
{
  return m_bNeedRecalc;
}

inline bool CElectricFieldData::get_enable_field() const
{
  return m_bEnable;
}

inline DWORD_PTR CElectricFieldData::get_enable_field_ptr() const
{
  return (DWORD_PTR)&m_bEnable;
}

inline int CElectricFieldData::get_type() const
{
  return m_nType;
}

inline DWORD_PTR CElectricFieldData::get_type_ptr() const
{
  return (DWORD_PTR)&m_nType;
}

inline void CElectricFieldData::set_type(int nType)
{
  if(m_nType != nType)
  {
    m_nType = nType;
    m_bNeedRecalc = true;
  }
}

inline int CElectricFieldData::get_calc_method() const
{
  return m_nCalcMethod;
}

inline DWORD_PTR CElectricFieldData::get_calc_method_ptr() const
{
  return (DWORD_PTR)&m_nCalcMethod;
}

inline void CElectricFieldData::set_calc_method(int nCalcMethod)
{
  if(m_nCalcMethod != nCalcMethod)
  {
    m_nCalcMethod = nCalcMethod;
    m_bNeedRecalc = true;
  }
}

inline double CElectricFieldData::get_scale() const
{
  return m_fScale / SI_to_CGS_Voltage;
}

inline DWORD_PTR CElectricFieldData::get_scale_ptr() const
{
  return (DWORD_PTR)&m_fScale;
}

inline void CElectricFieldData::set_scale(double fScale)
{
  m_fScale = fScale * SI_to_CGS_Voltage;
}

inline double CElectricFieldData::get_freq() const
{
  return m_fOmega / Const_2PI;
}

inline DWORD_PTR CElectricFieldData::get_freq_ptr() const
{
  return (DWORD_PTR)&m_fOmega;
}

inline void CElectricFieldData::set_freq(double fFreq)
{
  m_fOmega = Const_2PI * fFreq;
}

inline UINT CElectricFieldData::get_iter_count() const
{
  return m_nIterCount;
}

inline DWORD_PTR CElectricFieldData::get_iter_count_ptr() const
{
  return (DWORD_PTR)&m_nIterCount;
}

inline void CElectricFieldData::set_iter_count(UINT nCount)
{
  if(m_nIterCount != nCount)
  {
    m_nIterCount = nCount;
    m_bNeedRecalc = true;
  }
}

inline bool CElectricFieldData::get_analyt_field() const
{
  return m_bAnalytField;
}

inline DWORD_PTR CElectricFieldData::get_analyt_field_ptr() const
{
  return (DWORD_PTR)&m_bAnalytField;
}

inline double CElectricFieldData::get_inscr_radius() const
{
  return m_fRadius;
}

inline DWORD_PTR CElectricFieldData::get_inscr_radius_ptr() const
{
  return (DWORD_PTR)&m_fRadius;
}

inline void CElectricFieldData::set_inscr_radius(double fRadius)
{
  if(m_fRadius != fRadius)
  {
    m_fRadius = fRadius;
    m_bNeedRecalc = true;
  }
}

inline double CElectricFieldData::get_low_analyt_lim() const
{
  return m_fLowLimX;
}

inline DWORD_PTR CElectricFieldData::get_low_analyt_lim_ptr() const
{
  return (DWORD_PTR)&m_fLowLimX;
}

inline void CElectricFieldData::set_low_analyt_lim(double fLowLimX)
{
  if(m_fLowLimX != fLowLimX)
  {
    m_fLowLimX = fLowLimX;
    m_bNeedRecalc = true;
  }
}

inline CString CElectricFieldData::get_field_name() const
{
  return m_sName;
}

inline DWORD_PTR CElectricFieldData::get_field_name_ptr() const
{
  return (DWORD_PTR)&m_sName;
}

inline void CElectricFieldData::set_field_name(CString sName)
{
  m_sName = sName;
}

inline size_t CElectricFieldData::get_bc_count() const
{
  return m_vBoundCond.size();
}

inline CPotentialBoundCond* CElectricFieldData::get_bc(size_t nId) const
{
  return nId < m_vBoundCond.size() ? m_vBoundCond.at(nId) : NULL;
}

inline void CElectricFieldData::invalidate()
{
  m_bNeedRecalc = true;
}

inline Vector3D CElectricFieldData::get_field(size_t nInd) const
{
  return m_bNeedRecalc ? Vector3D(0, 0, 0) : Vector3D(m_vField[nInd].x, m_vField[nInd].y, m_vField[nInd].z);
}

inline double CElectricFieldData::get_omega() const
{
  return m_fOmega;
}

inline double CElectricFieldData::get_ampl() const
{
  return m_fScale;
}

};  // namespace EvaporatingParticle.