
#include "stdafx.h"

#include "constant.hpp"
#include "ParticleTracking.h"


namespace EvaporatingParticle
{

//---------------------------------------------------------------------------------------
// CSelectedAreas - analog of Named Selections in Ansys.
//---------------------------------------------------------------------------------------
CSelectedAreas::CSelectedAreas()
{
  m_sName = default_name();
  m_bVisible = true;
  m_dwClr = clLtGray;
}

CString CSelectedAreas::default_name()
{
  CSelAreasColl* pSelAreasColl = CParticleTrackingApp::Get()->GetSelAreas();
  size_t nSelAreasCount = pSelAreasColl->size();

  CString sName;
  char buff[8];
  int i = nSelAreasCount + 1;
  while(true)
  {
    bool bProperName = true;
    sName = CString(_T("Selection ")) + CString(itoa(i, buff, 10));
    for(size_t j = 0; j < nSelAreasCount; j++)
    {
      CSelectedAreas* pSelAreas = pSelAreasColl->at(j);
      if(sName == pSelAreas->get_name())
      {
        bProperName = false;
        break;
      }
    }

    if(bProperName)
      break;

    i++;
  }

  return sName;
}

CString CSelectedAreas::merge_opt_name(int nMergeOpt)
{
  switch(nMergeOpt)
  {
    case optAdd: return CString(_T("Add"));
    case optSubst: return CString(_T("Substitute"));
    case optRem: return CString(_T("Subtract"));
  }

  return CString(_T(""));
}

void CSelectedAreas::merge_items(CStringVector& vDest, int nMergeOpt) const
{
  switch(nMergeOpt)
  {
    case optAdd: add_items(vDest); break;
    case optSubst: subst_items(vDest); break;
    case optRem: remove_items(vDest); break;
  }
}

void CSelectedAreas::add_items(CStringVector& vDest) const
{
  size_t nSize = size();
  for(size_t i = 0; i < nSize; i++)
  {
    if(std::find(vDest.begin(), vDest.end(), at(i)) == vDest.end())
      vDest.push_back(at(i));
  }
}

void CSelectedAreas::subst_items(CStringVector& vDest) const
{
  vDest.clear();
  size_t nSize = size();
  for(size_t i = 0; i < nSize; i++)
    vDest.push_back(at(i));
}

void CSelectedAreas::remove_items(CStringVector& vDest) const
{
  size_t nSize = size();
  for(size_t i = 0; i < nSize; i++)
  {
    std::vector<std::string>::iterator iter = std::find(vDest.begin(), vDest.end(), at(i));
    if(iter != vDest.end())
      vDest.erase(iter);
  }
}

void CSelectedAreas::set_faces_color(COLORREF dwClr)
{
  if(m_dwClr != dwClr)
  {
    m_dwClr = dwClr;
    CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
    pDrawObj->invalidate_faces();
  }
}

bool CSelectedAreas::reg_is_selected(const std::string& sRegName) const
{
  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    if(at(i) == sRegName)
      return true;
  }

  return false;
}

CString CSelectedAreas::get_square_str() const
{
  double fSquare = get_square();
  if(fSquare < Const_Almost_Zero)
    return CString("");

  fSquare *= 100; // mm^2
  double fTrunc = 0.01 * int(0.5 + 100 * fSquare);
  std::ostringstream buffer;
  buffer << fTrunc;
  return buffer.str().c_str();
}

double CSelectedAreas::get_square() const
{
  double fSquare = 0;
  CRegion* pReg = NULL;
  CFace* pFace = NULL;
  size_t nRegCount = size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    pReg = CAnsysMesh::get_region(at(i));
    if(pReg == NULL)
      continue;

    size_t nFaceCount = pReg->vFaces.size();
    for(size_t j = 0; j < nFaceCount; j++)
    {
      pFace = pReg->vFaces.at(j);
      fSquare += pFace->square();
    }
  }

  return fSquare; // CGS, cm^2
}

void CSelectedAreas::save(CArchive& ar)
{
  const UINT nVersion = 1;  // 1 - m_dwClr.
  ar << nVersion;

  ar << m_sName;
  ar << m_bVisible;

  size_t nCount = size();
  ar << nCount;
  for(size_t i = 0; i < nCount; i++)
  {
    CString sName(at(i).c_str());
    ar << sName;
  }

  ar << m_dwClr;
}

void CSelectedAreas::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_sName;
  ar >> m_bVisible;

  clear();
  size_t nCount;
  ar >> nCount;   // count of regions belonging to this area.
  reserve(nCount);
  for(size_t i = 0; i < nCount; i++)
  {
    CString sName;
    ar >> sName;
    push_back(std::string((const char*)sName));

// Set visibility flag for the region:
    CRegion* pReg = CAnsysMesh::get_region(at(i));
    if(pReg != NULL)
      pReg->bEnabled = m_bVisible;
  }

  if(nVersion >= 1)
    ar >> m_dwClr;
}

//---------------------------------------------------------------------------------------
// CSelAreasColl - vector of CSelectedAreas.
//---------------------------------------------------------------------------------------
void CSelAreasColl::save(CArchive& ar)
{
  const UINT nVersion = 1;  // since 1 the default area is saved.
  ar << nVersion;

  size_t nCount = size();
  ar << nCount;
  for(size_t i = 0; i < nCount; i++)
  {
    CSelectedAreas* pSelAreas = at(i);
    pSelAreas->save(ar);
  }

  m_DefaultArea.save(ar);
}

void CSelAreasColl::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  clear();
  size_t nCount;
  ar >> nCount;
  reserve(nCount);
  for(size_t i = 0; i < nCount; i++)
  {
    CSelectedAreas* pSelAreas = new CSelectedAreas();
    pSelAreas->load(ar);
    push_back(pSelAreas);
  }

  if(nVersion >= 1)
    m_DefaultArea.load(ar);
}

void CSelAreasColl::clear_all()
{
  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
    delete at(i);

  clear();
}

void CSelAreasColl::populate_default_area()
{
  m_DefaultArea.clear();
  m_DefaultArea.set_name("Default Area");
  CRegionsCollection& vRegions = CParticleTrackingApp::Get()->GetTracker()->get_regions(false);
  size_t nRegCount = vRegions.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    CRegion* pReg = vRegions.at(i);
    std::string sRegName = pReg->sName;
    if(selected(sRegName))
      continue;

    pReg->bEnabled = m_DefaultArea.get_visibility_flag();
    m_DefaultArea.push_back(sRegName);
  }
}

bool CSelAreasColl::selected(const std::string& sRegName) const
{
  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    if(at(i)->reg_is_selected(sRegName))
      return true;
  }

  return false;
}

void CSelAreasColl::make_all_visible()
{
  bool* pVisFlag = NULL;
  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    pVisFlag = (bool*)(at(i)->get_visibility_flag_ptr());
    *pVisFlag = true;
  }

  pVisFlag = (bool*)(m_DefaultArea.get_visibility_flag_ptr());
  *pVisFlag = true;
}

};  // namespace EvaporatingParticle