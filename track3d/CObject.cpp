
#include "stdafx.h"
#include "CObject.h"
#include "MainFrm.h"

#include "ParticleTracking.h"


namespace EvaporatingParticle
{

CObject::CObject()
  : m_hProgressBarHandle(NULL), m_hJobNameHandle(NULL)
{
}

CObject::~CObject()
{
}

void CObject::terminate(bool bTermFlag)
{
  CParticleTrackingApp::Get()->terminate(bTermFlag);
}

bool CObject::get_terminate_flag() const
{
  return CParticleTrackingApp::Get()->is_terminated();
}

CString CObject::compile_string(const CStringVector& vNames)
{
  size_t nCount = vNames.size();
  char buff[16];
  std::string cStr(itoa(nCount, buff, 10));
  cStr += "  regions";
  
/*
  std::string cStr("");
  for(size_t j = 0; j < nCount; j++)
  {
    cStr += vNames.at(j);
    if(j < nCount - 1)
      cStr += ", ";
  }
*/
  return CString(cStr.c_str());
}

void CObject::set_status(const char* cAction, int nPercent) const
{
  CMainFrame* pMainWnd = (CMainFrame*)(CParticleTrackingApp::Get()->m_pMainWnd);
  if(pMainWnd == NULL)
    return;

  CMFCStatusBar* pStatusBar = pMainWnd->GetStatusBar();
  if(pStatusBar == NULL)
    return;

  char buff[4];
  std::string cPercentage("");
  if((nPercent >= 0) && (nPercent <= 100) && (_itoa_s(nPercent, buff, 4, 10) == 0))
    cPercentage = std::string("  ") + std::string(buff) + std::string(" %");

  std::string cStatusLine = std::string(cAction) + cPercentage;

  pStatusBar->SetPaneText(0, cStatusLine.c_str());
}

};  // namespace EvaporatingParticle