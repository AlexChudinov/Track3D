
#include "stdafx.h"

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

void CSelectedAreas::save(CArchive& ar)
{
  const UINT nVersion = 0;
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
}

void CSelectedAreas::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_sName;
  ar >> m_bVisible;

  clear();
  size_t nCount;
  ar >> nCount;
  reserve(nCount);
  for(size_t i = 0; i < nCount; i++)
  {
    CString sName;
    ar >> sName;
    push_back(std::string((const char*)sName));
  }
}

//---------------------------------------------------------------------------------------
// CSelAreasColl - vector of CSelectedAreas.
//---------------------------------------------------------------------------------------
void CSelAreasColl::save(CArchive& ar)
{
  const UINT nVersion = 0;
  ar << nVersion;

  size_t nCount = size();
  ar << nCount;
  for(size_t i = 0; i < nCount; i++)
  {
    CSelectedAreas* pSelAreas = at(i);
    pSelAreas->save(ar);
  }
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
}

void CSelAreasColl::clear_all()
{
  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
    delete at(i);

  clear();
}

};  // namespace EvaporatingParticle