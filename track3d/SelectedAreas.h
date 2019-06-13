#pragma once

#include "CObject.h"

namespace EvaporatingParticle
{

//---------------------------------------------------------------------------------------
// CSelectedAreas - analog of Named Selections in Ansys.
//---------------------------------------------------------------------------------------
class CSelectedAreas : public CStringVector  // analog of Named Selections in Ansys.
{
public:
  CSelectedAreas();

  CString           get_name() const;
  DWORD_PTR         get_name_ptr() const;
  void              set_name(const CString& sName);

  bool              get_visibility_flag() const;
  DWORD_PTR         get_visibility_flag_ptr() const;

  static CString    default_name();
  static CString    merge_opt_name(int nMergeOpt);

  enum  // Merge options.
  {
    optAdd    = 0,
    optSubst  = 1,
    optRem    = 2,
    optCount  = 3
  };

  void              merge_items(CStringVector& vDest, int nMergeOpt = optAdd) const;

  void              save(CArchive& ar);
  void              load(CArchive& ar);

protected:
  void              add_items(CStringVector& vDest) const;
  void              subst_items(CStringVector& vDest) const;
  void              remove_items(CStringVector& vDest) const;

private:
  CString           m_sName;      // the whole set of regions name, e.g. "RF Plus Electrodes".
  bool              m_bVisible;   // visibility of the whole set of regions.
};

//---------------------------------------------------------------------------------------
// CSelAreasColl - vector of CSelectedAreas.
//---------------------------------------------------------------------------------------
class CSelAreasColl : public std::vector<CSelectedAreas*>
{
public:
  void        save(CArchive& ar);
  void        load(CArchive& ar);

  void        clear_all();
};

inline CString CSelectedAreas::get_name() const
{
  return m_sName;
}

inline DWORD_PTR CSelectedAreas::get_name_ptr() const
{
  return (DWORD_PTR)&m_sName;
}

inline void CSelectedAreas::set_name(const CString& sName)
{
  m_sName = sName;
}

inline bool CSelectedAreas::get_visibility_flag() const
{
  return m_bVisible;
}

inline DWORD_PTR CSelectedAreas::get_visibility_flag_ptr() const
{
  return (DWORD_PTR)&m_bVisible;
}

//typedef std::vector<CSelectedAreas*> CSelAreasColl;

};  // namespace EvaporatingParticle