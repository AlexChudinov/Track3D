// .\ui\CExecutionDialog.cpp : implementation file
//

#include "stdafx.h"
#include "ExecutionDialog.h"
#include "ParticleTracking.h"
#include "Tracker.hpp"

using namespace EvaporatingParticle;

//---------------------------------------------------------------------------------------
// CExecutionDialog dialog
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CExecutionDialog, CDialog)

CExecutionDialog::CExecutionDialog(AFX_THREADPROC pThreadFunc, EvaporatingParticle::CObject* pObj, CWnd* pParent /*=NULL*/)
	: CDialog(CExecutionDialog::IDD, pParent)
{
  m_pThread = AfxBeginThread(pThreadFunc, (LPVOID)this, 0, 0, CREATE_SUSPENDED);
  m_pObj = pObj;
}

CExecutionDialog::~CExecutionDialog()
{
  m_pThread = NULL;
  m_pObj = NULL;
}

void CExecutionDialog::DoDataExchange(CDataExchange* pDX)
{
  CDialog::DoDataExchange(pDX);
  DDX_Control(pDX, IDC_TRACKING_PROGRESS, m_ProgressBar);
  DDX_Control(pDX, IDC_JOB_NAME, m_JobNameCtrl);
}

BOOL CExecutionDialog::OnInitDialog()
{
  CDialog::OnInitDialog();

  if(m_pObj == NULL)
    return TRUE;

  m_pObj->set_handlers(m_JobNameCtrl.m_hWnd, m_ProgressBar.m_hWnd, m_hWnd);
  m_ProgressBar.SetRange(0, 100);

  m_pThread->ResumeThread();

  return TRUE;
}

BEGIN_MESSAGE_MAP(CExecutionDialog, CDialog)
  ON_BN_CLICKED(IDOK, &CExecutionDialog::OnBnClickedOk)
END_MESSAGE_MAP()


// CExecutionDialog message handlers

void CExecutionDialog::OnBnClickedOk()
{
  if(m_pObj != NULL)
  {
    m_pObj->set_handlers(NULL, NULL);
    m_pObj->terminate();
  }

  OnOK();
}
