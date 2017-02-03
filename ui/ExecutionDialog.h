#pragma once

#include "resource.h"
#include "afxcmn.h"
#include "afxwin.h"

#include "CObject.h"

//---------------------------------------------------------------------------------------
// CExecutionDialog dialog
//---------------------------------------------------------------------------------------
class CExecutionDialog : public CDialog
{
	DECLARE_DYNAMIC(CExecutionDialog)

public:
  CExecutionDialog(AFX_THREADPROC pThreadFunc, EvaporatingParticle::CObject* pObj, CWnd* pParent = NULL);   // standard constructor
	virtual ~CExecutionDialog();

  virtual BOOL                  OnInitDialog();

  EvaporatingParticle::CObject* GetDialogObject() const;

// Dialog Data
	enum { IDD = IDD_TERMINATE_DIALOG };

private:
  CWinThread*                   m_pThread;
  EvaporatingParticle::CObject* m_pObj;

protected:
	virtual void                  DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()

public:
  afx_msg void                  OnBnClickedOk();

  CProgressCtrl                 m_ProgressBar;
  CEdit                         m_JobNameCtrl;
};

inline EvaporatingParticle::CObject* CExecutionDialog::GetDialogObject() const
{
  return m_pObj;
}
