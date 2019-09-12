
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
class CPotentialBoundCond
{
public:
  CPotentialBoundCond(BoundaryMesh::BoundaryType type = BoundaryMesh::FIXED_VAL, int val = fvPlusUnity);

  enum  // fixed value types:
  {
    fvPlusUnity     = 0,
    fvMinusUnity    = 1,
    fvLinearStepsX  = 2,
    fvLinearStepsY  = 3,
    fvQuadricStepsX = 4,
    fvQuadricStepsY = 5,
    fvCoulomb       = 6,    // attempt to calculate a mirror Coulomb field.
    fvCount         = 7
  };

// All changes in boundary conditions MUST result in the field re-calculation. The following GET-functions return "true" if something changes.
  int                           get_bc_type() const;
  DWORD_PTR                     get_bc_type_ptr() const;
  bool                          set_bc_type(int nNewType);

  int                           get_fixed_val_type() const;
  DWORD_PTR                     get_fixed_val_type_ptr() const;
  bool                          set_fixed_val_type(int nNewType);

  UINT                          get_steps_count() const;
  DWORD_PTR                     get_steps_count_ptr() const;
  bool                          set_steps_count(UINT nNewCount);

  double                        get_center_first_electr() const;
  DWORD_PTR                     get_center_first_electr_ptr() const;
  bool                          set_center_first_electr(double fNewCenter);

  double                        get_step_coord() const;
  DWORD_PTR                     get_step_coord_ptr() const;
  bool                          set_step_coord(double fNewStep);

  double                        get_start_coord() const;
  DWORD_PTR                     get_start_coord_ptr() const;
  bool                          set_start_coord(double fNewCoord);

  double                        get_end_coord() const;
  DWORD_PTR                     get_end_coord_ptr() const;
  bool                          set_end_coord(double fNewCoord);

  double                        get_start_phi() const;
  DWORD_PTR                     get_start_phi_ptr() const;
  bool                          set_start_phi(double fNewPhi);

  double                        get_end_phi() const;
  DWORD_PTR                     get_end_phi_ptr() const;
  bool                          set_end_phi(double fNewPhi);

  double                        get_first_dphi() const;
  DWORD_PTR                     get_first_dphi_ptr() const;
  bool                          set_first_dphi(double fNewdPhi);

  double                        get_last_dphi() const;
  DWORD_PTR                     get_last_dphi_ptr() const;
  bool                          set_last_dphi(double fNewdPhi);

  CString                       get_name() const;
  void                          set_name(const CString& sName);

  CStringVector&                get_region_names();
  const CStringVector&          get_region_names() const;
  DWORD_PTR                     get_region_names_ptr() const;

  bool                          get_visibility_flag() const;
  DWORD_PTR                     get_visibility_flag_ptr() const;

// The existing regions can be merged with those in selected Named Areas. Three merge options are available: add, substitute and remove.
  int                           get_merge_option() const;
  DWORD_PTR                     get_merge_option_ptr() const;
  void                          set_merge_option(int nOpt);

  CString                       get_last_merged() const;
  void                          set_last_merged(const CString& sName);

  double                        linear_potential(double x);

  static const char*            get_bc_type_name(int nType);
  static const char*            get_fixed_value_name(int nType);

  enum // user interface control type.
  {
    uitStartX       = 0,
    uitStepX        = 1,
    uitEndX         = 2,
    uitStartPhi     = 3,
    uitEndPhi       = 4,
    uitFirstStepPhi = 5,
    uitLastStepPhi  = 6,
    uitStepsCount   = 7,
    uitCenterFirst  = 8
  };

  static const char*            get_control_title(int nFixedValType, int nCtrlType);
  static const char*            get_hint(int nFixedValType, int nCtrlType);

  void                          save(CArchive& ar);
  void                          load(CArchive& ar);

private:
  int                           m_nType;

  CStringVector                 m_vRegNames;      // names of the regions with non-trivial boundary conditions.
  bool                          m_bVisible;       // visibility status of the selected regions.
// A flag showing how to merge the existing regions with the ones contained in the selected Named Area. Can be "add", "substitute" and "remove", see SelectedAreas.h.
  int                           m_nMergeOpt;
// The name of the Named Area last merged:
  CString                       m_sLastMerged;

  int                           m_nFixedValType;  // the value at the boundary can be +1V, -1V, step-wise potential or Coulomb potential.
  CString                       m_sName;

// Step-like potential:
  UINT                          m_nStepsCount;
  double                        m_fCenterFirstElectr,
                                m_fStepX;

// For linear potential the following three parameters strongly define the linear function phi(x). The center of the first electrode
// is not always equal to fStartX. Two or more step-like potentials can be linked to one the same linear phi(x).
  double                        m_fStartX,
                                m_fEndX,
// Dimensionless value of potential at the end of the interval. Corresponding "fStartPhi" is assumed always to be unity. Linear potential only.
                                m_fEndPhi;

// Step-like potential (quadric):
  double                        m_fStartPhi,      // absolute value of the potential on the first electrode.
                                m_fFirstStepPhi,  // absolute value of the potential gradient between the centers of the first and second steps.
                                m_fLastStepPhi;   // absolute value of the potential gradient between the centers of the last but one and last steps.
};

inline int CPotentialBoundCond::get_bc_type() const
{
  return m_nType;
}

inline DWORD_PTR CPotentialBoundCond::get_bc_type_ptr() const
{
  return (DWORD_PTR)&m_nType;
}

inline bool CPotentialBoundCond::set_bc_type(int nNewType)
{
  if(m_nType != nNewType)
  {
    m_nType = nNewType;
    return true;
  }
  return false;
}

inline int CPotentialBoundCond::get_fixed_val_type() const
{
  return m_nFixedValType;
}

inline DWORD_PTR CPotentialBoundCond::get_fixed_val_type_ptr() const
{
  return (DWORD_PTR)&m_nFixedValType;
}

inline bool CPotentialBoundCond::set_fixed_val_type(int nNewType)
{
  if(m_nFixedValType != nNewType)
  {
    m_nFixedValType = nNewType;
    return true;
  }
  return false;
}

inline CString CPotentialBoundCond::get_name() const
{
  return m_sName;
}

inline void CPotentialBoundCond::set_name(const CString& sName)
{
  m_sName = sName;
}

inline CStringVector& CPotentialBoundCond::get_region_names()
{
  return m_vRegNames;
}

inline const CStringVector& CPotentialBoundCond::get_region_names() const
{
  return m_vRegNames;
}

inline DWORD_PTR CPotentialBoundCond::get_region_names_ptr() const
{
  return (DWORD_PTR)&m_vRegNames;
}

inline bool CPotentialBoundCond::get_visibility_flag() const
{
  return m_bVisible;
}

inline DWORD_PTR CPotentialBoundCond::get_visibility_flag_ptr() const
{
  return (DWORD_PTR)&m_bVisible;
}

// Step-like potential:
inline UINT CPotentialBoundCond::get_steps_count() const
{
  return m_nStepsCount;
}

inline DWORD_PTR CPotentialBoundCond::get_steps_count_ptr() const
{
  return (DWORD_PTR)&m_nStepsCount;
}

inline bool CPotentialBoundCond::set_steps_count(UINT nNewCount)
{
  if(m_nStepsCount != nNewCount)
  {
    m_nStepsCount = nNewCount;
    return true;
  }
  return false;
}

inline double CPotentialBoundCond::get_center_first_electr() const
{
  return m_fCenterFirstElectr;
}

inline DWORD_PTR CPotentialBoundCond::get_center_first_electr_ptr() const
{
  return (DWORD_PTR)&m_fCenterFirstElectr;
}

inline bool CPotentialBoundCond::set_center_first_electr(double fNewCenter)
{
  if(m_fCenterFirstElectr != fNewCenter)
  {
    m_fCenterFirstElectr = fNewCenter;
    return true;
  }
  return false;
}

inline double CPotentialBoundCond::get_step_coord() const
{
  return m_fStepX;
}

inline DWORD_PTR CPotentialBoundCond::get_step_coord_ptr() const
{
  return (DWORD_PTR)&m_fStepX;
}

inline bool CPotentialBoundCond::set_step_coord(double fNewStep)
{
  if(m_fStepX != fNewStep)
  {
    m_fStepX = fNewStep;
    return true;
  }
  return false;
}

inline double CPotentialBoundCond::get_start_coord() const
{
  return m_fStartX;
}

inline DWORD_PTR CPotentialBoundCond::get_start_coord_ptr() const
{
  return (DWORD_PTR)&m_fStartX;
}

inline bool CPotentialBoundCond::set_start_coord(double fNewCoord)
{
  if(m_fStartX != fNewCoord)
  {
    m_fStartX = fNewCoord;
    return true;
  }
  return false;
}

inline double CPotentialBoundCond::get_end_coord() const
{
  return m_fEndX;
}

inline DWORD_PTR CPotentialBoundCond::get_end_coord_ptr() const
{
  return (DWORD_PTR)&m_fEndX;
}

inline bool CPotentialBoundCond::set_end_coord(double fNewCoord)
{
  if(m_fEndX != fNewCoord)
  {
    m_fEndX = fNewCoord;
    return true;
  }
  return false;
}

inline double CPotentialBoundCond::get_start_phi() const
{
  return m_fStartPhi;
}

inline DWORD_PTR CPotentialBoundCond::get_start_phi_ptr() const
{
  return (DWORD_PTR)&m_fStartPhi;
}

inline bool CPotentialBoundCond::set_start_phi(double fNewPhi)
{
  if(m_fStartPhi != fNewPhi)
  {
    m_fStartPhi = fNewPhi;
    return true;
  }
  return false;
}

inline double CPotentialBoundCond::get_end_phi() const
{
  return m_fEndPhi;
}

inline DWORD_PTR CPotentialBoundCond::get_end_phi_ptr() const
{
  return (DWORD_PTR)&m_fEndPhi;
}

inline bool CPotentialBoundCond::set_end_phi(double fNewPhi)
{
  if(m_fEndPhi != fNewPhi)
  {
    m_fEndPhi = fNewPhi;
    return true;
  }
  return false;
}

inline double CPotentialBoundCond::get_first_dphi() const
{
  return m_fFirstStepPhi;
}

inline DWORD_PTR CPotentialBoundCond::get_first_dphi_ptr() const
{
  return (DWORD_PTR)&m_fFirstStepPhi;
}

inline bool CPotentialBoundCond::set_first_dphi(double fNewdPhi)
{
  if(m_fFirstStepPhi != fNewdPhi)
  {
    m_fFirstStepPhi = fNewdPhi;
    return true;
  }
  return false;
}

inline double CPotentialBoundCond::get_last_dphi() const
{
  return m_fLastStepPhi;
}

inline DWORD_PTR CPotentialBoundCond::get_last_dphi_ptr() const
{
  return (DWORD_PTR)&m_fLastStepPhi;
}

inline bool CPotentialBoundCond::set_last_dphi(double fNewdPhi)
{
  if(m_fLastStepPhi != fNewdPhi)
  {
    m_fLastStepPhi = fNewdPhi;
    return true;
  }
  return false;
}

inline int CPotentialBoundCond::get_merge_option() const
{
  return m_nMergeOpt;
}

inline DWORD_PTR CPotentialBoundCond::get_merge_option_ptr() const
{
  return (DWORD_PTR)&m_nMergeOpt;
}

inline void CPotentialBoundCond::set_merge_option(int nOpt)
{
  m_nMergeOpt = nOpt;
}

inline CString CPotentialBoundCond::get_last_merged() const
{
  return m_sLastMerged;
}

inline void CPotentialBoundCond::set_last_merged(const CString& sName)
{
  m_sLastMerged = sName;
}

struct CNode3D;
struct CRegion;
class CFiniteVolumesSolver;
typedef std::vector<CPotentialBoundCond*> CPotentialBoundCondColl;
//---------------------------------------------------------------------------------------
// CElectricFieldData - a class for simulating scalable electric fields.
//---------------------------------------------------------------------------------------
class CElectricFieldData : public CObject
{
public:

  enum
  {
    cmLaplacian3    = 0,
    cmDirTessLap3   = 1,
    cmFinVolJacobi  = 2,
	cmEigenLibSolver= 3,
    cmCount         = 4
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

  bool                    get_enable_multithread() const;
  DWORD_PTR               get_enable_multithread_ptr() const;

  bool                    get_enable_vis() const;
  DWORD_PTR               get_enable_vis_ptr() const;

  double                  get_tol() const;
  DWORD_PTR               get_tol_ptr() const;
  void                    set_tol(double fTol);

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

  float                   get_phi(size_t nInd) const;

// For Coulomb potential visualization only. These functions must be called only from CTracker::average_clmb_field(UINT nIter).
  void                    init_clmb_phi(size_t nNodeCount);      // here m_vClmbPhi is initialized ...
  void                    set_clmb_phi(size_t nInd, float fPhi); // ... and modified at every iteration.
  float                   get_clmb_phi(size_t nInd) const;

  Vector3D                get_field(size_t nInd) const;

  double                  get_omega() const;
  double                  get_ampl() const;   // returns the scale in CGS units.

  void                    add_bc();
  void                    remove_bc(size_t nId);

// Stability parameter q; the ion parameters are taken from CTracker. All input data are in CGSE.
  static double           get_stab_param(double fAmpl, double fFreq, double fInscR);

  static const char*      get_field_type_name(int nType);
  static const char*      get_calc_method_name(int nCalcMethod);

  bool                    calc_field();

  bool                    need_recalc() const;
  void                    invalidate();

  void                    output_convergence_history(const std::vector<float>& vMaxRelErr);

  void                    save(CArchive& ar);
  void                    load(CArchive& ar);

protected:
  void                    set_default();
  CString                 default_name();
  void                    clear_bc();

  bool                    is_selected(CRegion* pReg) const;           // returns true if the region is selected for boundary conditions.
  CIndexVector            get_reg_nodes(CRegion* pReg) const;         // returns the vector of global indices of the region nodes.

  void                    check_regions();  // run over all boundary conditions and check if all region names correspond to relevant regions.

  bool                    set_boundary_conditions(CMeshAdapter& mesh);
  bool                    set_boundary_conditions(CFiniteVolumesSolver& solver);

  bool                    set_default_boundary_conditions(CMeshAdapter& mesh);
  bool                    set_default_boundary_conditions(CFiniteVolumesSolver& solver);

  void                    set_boundary_values(CMeshAdapter& mesh, CRegion* pReg, CPotentialBoundCond* pBC = NULL);

public:
  double                  linear_step_potential(CPotentialBoundCond* pBC, const Vector3D& vPos) const;  // linear step-wise boundary conditions support.
  double                  quadric_step_potential(CPotentialBoundCond* pBC, const Vector3D& vPos) const; // parabolic step-wise boundary conditions support.
protected:
  bool                    coulomb_potential(const CIndexVector& vNodeIds, std::vector<float>& vPhi) const;  // Coulomb boundary conditions support.

  Vector3D                calc_norm(CNode3D* pNode) const;

  bool                    get_result() const;
  void                    notify_scene(); // let the scene objects know that the potential field has been changed.

  void                    apply_analytic_field(const Vector3D& vPos, Vector3F& vField, float& fPhi);

  bool                    calc_lap3();
  bool                    calc_dirichlet_lap3();
  bool                    calc_finite_vol_jacobi();
  bool					  calc_eigen_lib_lap();

  enum  { jobDfltCond  = 0, jobUserCond  = 1,  jobCalcField = 2 };

  const char*             job_name(int nJobType) const;

private:
  bool                    m_bEnable,
                          m_bMultiThread, // for tests, never saved to the stream.
                          m_bEnableVis;   // for visualization; if true, the potential of this field will be added to node.phi for every node in the mesh. 

  int                     m_nCalcMethod,
                          m_nType;

  double                  m_fScale,       // potential scale, in V.
                          m_fOmega;       // circular frequency, for radio-frequency fields only.

  UINT                    m_nIterCount;
  double                  m_fTol;       // tolerance: if the maxinal relative error becomes less than m_fTol, the solution terminates.

  CPotentialBoundCondColl m_vBoundCond;

// An attempt to get analytic field in the flatapole. Alpha version.
  bool                    m_bAnalytField;

  double                  m_fRadius,    // inscribed radius of the flatapole electrodes.
                          m_fLowLimX,   // an analytic formula will be used if m_fLowLimX < x < m_fHighLimX;
                          m_fHighLimX;

  std::vector<float>      m_vPhi;       // for field visualization and test purposes.
  std::vector<float>      m_vClmbPhi;   // for visualization of the Coulomb potential, which is accumulated at every iteration.
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

inline bool CElectricFieldData::get_enable_multithread() const
{
  return m_bMultiThread;
}

inline DWORD_PTR CElectricFieldData::get_enable_multithread_ptr() const
{
  return (DWORD_PTR)&m_bMultiThread;
}

inline bool CElectricFieldData::get_enable_vis() const
{
  return m_bEnableVis;
}

inline DWORD_PTR CElectricFieldData::get_enable_vis_ptr() const
{
  return (DWORD_PTR)&m_bEnableVis;
}

inline double CElectricFieldData::get_tol() const
{
  return m_fTol;
}

inline DWORD_PTR CElectricFieldData::get_tol_ptr() const
{
  return (DWORD_PTR)&m_fTol;
}

inline void CElectricFieldData::set_tol(double fTol)
{
  if(m_fTol != fTol)
  {
    m_fTol = fTol;
    m_bNeedRecalc = true;
  }
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

inline float CElectricFieldData::get_phi(size_t nInd) const
{
  if(nInd >= m_vPhi.size())
    return 0.0f;

  return m_nType == typeMirror ? m_vPhi[nInd] : (m_bNeedRecalc ? 0.0f : m_vPhi[nInd]);
}

inline void CElectricFieldData::init_clmb_phi(size_t nNodeCount)
{
  if(m_nType != typeMirror)
    return;

  m_vClmbPhi.resize(nNodeCount, 0.0f);
}

inline float CElectricFieldData::get_clmb_phi(size_t nInd) const
{
  if(nInd >= m_vClmbPhi.size())
    return 0.0f;

  return m_nType == typeMirror ? m_vClmbPhi[nInd] : 0.0f;
}

inline void CElectricFieldData::set_clmb_phi(size_t nInd, float fPhi)
{
  if(m_nType != typeMirror)
    return;

  m_vClmbPhi[nInd] = fPhi;
}

inline Vector3D CElectricFieldData::get_field(size_t nInd) const
{
  return m_nType == typeMirror ? 
    Vector3D(m_vField[nInd].x, m_vField[nInd].y, m_vField[nInd].z) : 
    (m_bNeedRecalc ? Vector3D(0, 0, 0) : Vector3D(m_vField[nInd].x, m_vField[nInd].y, m_vField[nInd].z));
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