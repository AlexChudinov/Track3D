
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
  CPotentialBoundCond(BoundaryMesh::BoundaryType type = BoundaryMesh::FIXED_VAL, int val = fvPlusUnity)
    : nType(type), nFixedValType(val)
  {
  }

  enum  // fixed value types:
  {
    fvPlusUnity   = 0,
    fvMinusUnity  = 1,
    fvCount       = 2
  };

  BoundaryMesh::BoundaryType    nType;
  CStringVector                 vRegNames;      // names of the regions with non-trivial boundary conditions.
  int                           nFixedValType;  // the value at the boundary can be either +1V or -1V.
  std::string                   sName;

  static const char*            get_bc_type_name(BoundaryMesh::BoundaryType nType);
  static const char*            get_fixed_value_name(int nType);

  void                          save(CArchive& ar);
  void                          load(CArchive& ar);
};

struct CNode3D;
struct CRegion;
typedef std::vector<CPotentialBoundCond*> CPotentialBoundCondColl;
//---------------------------------------------------------------------------------------
// CElectricFieldData - an auxiliary class for simulating scalable electric fields.
//---------------------------------------------------------------------------------------
class CElectricFieldData : public CObject
{
public:

  enum
  {
    typeFieldDC = 0,
    typeFieldRF = 1,
    typeCount   = 2
  };

  CElectricFieldData(int nType = typeFieldDC);
  virtual ~CElectricFieldData();

  bool                    get_enable_field() const;
  DWORD_PTR               get_enable_field_ptr() const;

  int                     get_type() const;
  DWORD_PTR               get_type_ptr() const;
  void                    set_type(int nType);

  double                  get_scale() const;
  DWORD_PTR               get_scale_ptr() const;
  void                    set_scale(double fScale);

  double                  get_freq() const;
  DWORD_PTR               get_freq_ptr() const;
  void                    set_freq(double fFreq);

  UINT                    get_iter_count() const;
  DWORD_PTR               get_iter_count_ptr() const;
  void                    set_iter_count(UINT nCount);

  size_t                  get_bc_count() const;
  CPotentialBoundCond*    get_bc(size_t nId) const;

  Vector3D                get_field(size_t nInd) const;
  double                  get_omega() const;

  void                    add_bc();
  void                    remove_bc(size_t nId);

  static const char*      get_field_type_name(int nType);

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

  bool                    set_boundary_conditions(CMeshAdapter& mesh);
  bool                    set_default_boundary_conditions(CMeshAdapter& mesh);
  void                    set_boundary_values(CMeshAdapter& mesh, CRegion* pReg, CPotentialBoundCond* pBC = NULL);

  Vector3D                calc_norm(CNode3D* pNode) const;

  bool                    get_result(bool bTest) const;
  void                    notify_scene(); // let the scene objects know that the potential field has been changed.

private:
  bool                    m_bEnable;
  int                     m_nType;

  double                  m_fScale,       // potential scale, in V.
                          m_fOmega;       // circular frequency, for radio-frequency fields only.

  UINT                    m_nIterCount;

  CPotentialBoundCondColl m_vBoundCond;

  std::vector<float>      m_vPotential;
  std::vector<Vector3F>   m_vField;

// Run-time:
  bool                    m_bNeedRecalc;
};

//-------------------------------------------------------------------------------------------------
// CFieldDataCollection.
//-------------------------------------------------------------------------------------------------
class CFieldDataColl : public std::vector<CElectricFieldData*>
{
public:
  virtual ~CFieldDataColl();

  bool              calc_fields();
  bool              need_recalc() const;

  void              clear_fields_in_nodes();
  void              clear_fields();

  bool              sel_region_changed(CStringVector* pRegNames);
  bool              remove_bound_cond(CPotentialBoundCond* pBC);

  void              save(CArchive& ar);
  void              load(CArchive& ar);
};

inline bool CFieldDataColl::need_recalc() const
{
  for(size_t i = 0; i < size(); i++)
    if(at(i)->need_recalc())
      return true;

  return false;
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

inline double CElectricFieldData::get_scale() const
{
  return m_fScale;
}

inline DWORD_PTR CElectricFieldData::get_scale_ptr() const
{
  return (DWORD_PTR)&m_fScale;
}

inline void CElectricFieldData::set_scale(double fScale)
{
  m_fScale = fScale;
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
  m_nIterCount = nCount;
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

};  // namespace EvaporatingParticle.