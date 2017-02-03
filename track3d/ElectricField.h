
#pragma once

#include "CObject.h"
#include "constant.hpp"
#include <LSExport.h>
#include <vector>
#include <string>

namespace EvaporatingParticle
{

//---------------------------------------------------------------------------------------
// CPotentialBoundCond
//---------------------------------------------------------------------------------------
struct CPotentialBoundCond
{
  CPotentialBoundCond(PotentialField::BOUNDARY_TYPE type = PotentialField::FIXED_VAL, int val = fvPlusUnity)
    : nType(type), nFixedValType(val)
  {
  }

  enum  // fixed value types:
  {
    fvPlusUnity   = 0,
    fvMinusUnity  = 1,
    fvCount       = 2
  };

  PotentialField::BOUNDARY_TYPE nType;
  CStringVector                 vRegNames;      // names of the regions with non-trivial boundary conditions.
  int                           nFixedValType;  // the value at the boundary can be either +1V or -1V.
  std::string                   sName;

  static const char*            get_bc_type_name(PotentialField::BOUNDARY_TYPE nType);
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

  void                    add_bc();
  void                    remove_bc(size_t nId);

  static const char*      get_field_type_name(int nType);

  void                    calc_field();
  void                    invalidate();

  void                    save(CArchive& ar);
  void                    load(CArchive& ar);

protected:
  void                    set_default();
  void                    clear_bc();

  Mesh*                   create_mesh();
  CRegion*                get_region(const std::string& sName) const; // returns a region pointer by its name or NULL if the name is not found.
  bool                    is_selected(CRegion* pReg) const;           // returns true if the region is selected for boundary conditions.

  void                    set_boundary_conditions(PotentialField* pField);
  void                    set_default_boundary_conditions(PotentialField* pField);
  void                    set_boundary_values(PotentialField* pField, CRegion* pReg, CPotentialBoundCond* pBC = NULL);

  V3D                     calc_norm(CNode3D* pNode) const;

  bool                    get_result(PotentialField* pField) const;
  void                    notify_scene(); // let the scene objects know that the potential field has been changed.

private:
  int                     m_nType;

  double                  m_fScale,       // potential scale, in V.
                          m_fOmega;       // circular frequency, for radio-frequency fields only.

  UINT                    m_nIterCount;

  CPotentialBoundCondColl m_vBoundCond;

// Temporarily, for DEBUG purposes only, the data contain potentials. Later, in the working
// regime the components of electric field will be stored in this class.
  std::vector<float>      m_vPotential;

// Run-time:
  bool                    m_bNeedRecalc,
                          m_bScaleChanged;
};

//-------------------------------------------------------------------------------------------------
// CFieldDataCollection.
//-------------------------------------------------------------------------------------------------
class CFieldDataColl : public std::vector<CElectricFieldData*>
{
public:
  virtual ~CFieldDataColl();

  void              clear_fields();
  void              calc_fields();

  bool              sel_region_changed(CStringVector* pRegNames);
  bool              remove_bound_cond(CPotentialBoundCond* pBC);

  void              save(CArchive& ar);
  void              load(CArchive& ar);
};


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
  if(m_fScale != fScale)
  {
    m_fScale = fScale;
    m_bScaleChanged = true;
  }
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

};  // namespace EvaporatingParticle.