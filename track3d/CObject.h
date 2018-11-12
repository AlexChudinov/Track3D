
#pragma once

#include <vector>
#include <string>

namespace EvaporatingParticle
{

typedef std::vector<std::string> CStringVector;

class CObject
{
public:
  CObject();
  virtual ~CObject();

  void            get_handlers(HWND& hJobName, HWND& hProgress, HWND& hDlgWnd) const;
  void            set_handlers(HWND hJobName, HWND hProgress, HWND hDlgWnd = NULL);

  bool            get_terminate_flag() const;
  void            terminate(bool bTermFlag = true);

  static CString  compile_string(const CStringVector& vStr);

  void            set_job_name(const char* pJobName);
  void            set_progress(int nPercent);

protected:
  void            show_dlg(int nCmdShow);

  HWND            m_hProgressBarHandle;
  HWND            m_hJobNameHandle;
  HWND            m_hDlgWndHandle;

  void            set_status(const char* cAction, int nPercent) const;
};

inline void CObject::get_handlers(HWND& hJobName, HWND& hProgress, HWND& hDlgWnd) const
{
  hJobName = m_hJobNameHandle;
  hProgress = m_hProgressBarHandle;
  hDlgWnd = m_hDlgWndHandle;
}

inline void CObject::set_handlers(HWND hJobName, HWND hProgress, HWND hDlgWnd)
{
  m_hJobNameHandle = hJobName;
  m_hProgressBarHandle = hProgress;
  m_hDlgWndHandle = hDlgWnd;
}

inline void CObject::set_job_name(const char* pJobName)
{
  CWnd* pEditCtrl = m_hJobNameHandle != NULL ? CWnd::FromHandle(m_hJobNameHandle) : NULL;
  if(pEditCtrl != NULL)
    pEditCtrl->SetWindowTextA(pJobName);
}

inline void CObject::set_progress(int nPercent)
{
  CProgressCtrl* pProgressBar = m_hProgressBarHandle != NULL ? (CProgressCtrl*)CWnd::FromHandle(m_hProgressBarHandle) : NULL;
  if(pProgressBar != NULL)
    pProgressBar->SetPos(nPercent);
}

inline void CObject::show_dlg(int nCmdShow)
{
  CWnd* pDlgWnd = m_hDlgWndHandle != NULL ? CWnd::FromHandle(m_hDlgWndHandle) : NULL;
  if(pDlgWnd != NULL)
    pDlgWnd->ShowWindow(nCmdShow);
}

}; // namespace EvaporatingParticle