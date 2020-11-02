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

  COLORREF          get_faces_color() const;
  DWORD_PTR         get_faces_color_ptr() const;
  void              set_faces_color(COLORREF dwClr);

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

  bool              reg_is_selected(const std::string& sRegName) const;

  double            get_square() const;
  CString           get_square_str() const;

protected:
  void              add_items(CStringVector& vDest) const;
  void              subst_items(CStringVector& vDest) const;
  void              remove_items(CStringVector& vDest) const;

private:
  CString           m_sName;      // the whole set of regions name, e.g. "RF Plus Electrodes".
  bool              m_bVisible;   // visibility of the whole set of regions.
  COLORREF          m_dwClr;
};

//---------------------------------------------------------------------------------------
// CSelAreasColl - vector of CSelectedAreas.
//---------------------------------------------------------------------------------------
class CSelAreasColl : public std::vector<CSelectedAreas*>
{
public:
  void            save(CArchive& ar);
  void            load(CArchive& ar);

  void            clear_all();
  void            make_all_visible();

  CSelectedAreas* get_default_area();

protected:
  CSelectedAreas  m_DefaultArea;  // a run-time object to control the color and visibility of the unselected regions.

  void            populate_default_area();
  bool            selected(const std::string& sRegName) const;
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

inline COLORREF CSelectedAreas::get_faces_color() const
{
  return m_dwClr;
}

inline DWORD_PTR CSelectedAreas::get_faces_color_ptr() const
{
  return (DWORD_PTR)&m_dwClr;
}


inline CSelectedAreas* CSelAreasColl::get_default_area()
{
  populate_default_area();
  return &m_DefaultArea;
}

};  // namespace EvaporatingParticle