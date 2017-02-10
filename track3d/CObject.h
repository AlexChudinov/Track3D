
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

  void            set_handlers(HWND hJobName, HWND hProgress, HWND hDlgWnd = NULL);

  bool            get_terminate_flag() const;
  void            terminate(bool bTermFlag = true);

  static CString  compile_string(const CStringVector& vStr);

  void            set_job_name(const char* pJobName);
  void            set_progress(int nPercent);

protected:
  void            show_dlg(int nCmdShow);

  const bool*     get_terminate_ptr() const;

  HWND            m_hProgressBarHandle;
  HWND            m_hJobNameHandle;
  HWND            m_hDlgWndHandle;
  bool            m_bTerminate;

  void            set_status(const char* cAction, int nPercent) const;
};

inline void CObject::set_handlers(HWND hJobName, HWND hProgress, HWND hDlgWnd)
{
  m_hJobNameHandle = hJobName;
  m_hProgressBarHandle = hProgress;
  m_hDlgWndHandle = hDlgWnd;
}

inline void CObject::set_job_name(const char* pJobName)
{
  CWnd* pEditCtrl = CWnd::FromHandle(m_hJobNameHandle);
  if(pEditCtrl != NULL)
    pEditCtrl->SetWindowTextA(pJobName);
}

inline void CObject::set_progress(int nPercent)
{
  CProgressCtrl* pProgressBar = (CProgressCtrl*)CWnd::FromHandle(m_hProgressBarHandle);
  if(pProgressBar != NULL)
    pProgressBar->SetPos(nPercent);
}

inline void CObject::terminate(bool bTermFlag)
{
  m_bTerminate = bTermFlag;
}

inline bool CObject::get_terminate_flag() const
{
  return m_bTerminate;
}

inline void CObject::show_dlg(int nCmdShow)
{
  CWnd* pDlgWnd = CWnd::FromHandle(m_hDlgWndHandle);
  if(pDlgWnd != NULL)
    pDlgWnd->ShowWindow(nCmdShow);
}

inline const bool* CObject::get_terminate_ptr() const
{
  return &m_bTerminate;
}

}; // namespace EvaporatingParticle