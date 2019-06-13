
// ParticleTracking.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "afxwinappex.h"
#include "ParticleTracking.h"
#include "MainFrm.h"

#include "ParticleTrackingDoc.h"
#include "ParticleTrackingView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CParticleTrackingApp

BEGIN_MESSAGE_MAP(CParticleTrackingApp, CWinAppEx)
	ON_COMMAND(ID_APP_ABOUT, &CParticleTrackingApp::OnAppAbout)
	// Standard file based document commands
	ON_COMMAND(ID_FILE_NEW, &CWinAppEx::OnFileNew)
	ON_COMMAND(ID_FILE_OPEN, &CWinAppEx::OnFileOpen)
END_MESSAGE_MAP()


// CParticleTrackingApp construction

CParticleTrackingApp::CParticleTrackingApp()
{

	m_bHiColorIcons = TRUE;

	// TODO: add construction code here,
  m_bTerminate = false; // global termination flag [MS] 10-11-2018.

	// Place all significant initialization in InitInstance
}

// The one and only CParticleTrackingApp object

CParticleTrackingApp theApp;


// CParticleTrackingApp initialization

BOOL CParticleTrackingApp::InitInstance()
{
	CWinAppEx::InitInstance();

	// Standard initialization
	// If you are not using these features and wish to reduce the size
	// of your final executable, you should remove from the following
	// the specific initialization routines you do not need
	// Change the registry key under which our settings are stored
	// TODO: You should modify this string to be something appropriate
	// such as the name of your company or organization
	SetRegistryKey(_T("Local AppWizard-Generated Applications"));
	LoadStdProfileSettings(4);  // Load standard INI file options (including MRU)

	InitContextMenuManager();

	InitKeyboardManager();

	InitTooltipManager();
	CMFCToolTipInfo ttParams;
	ttParams.m_bVislManagerTheme = TRUE;
	theApp.GetTooltipManager()->SetTooltipParams(AFX_TOOLTIP_TYPE_ALL,
		RUNTIME_CLASS(CMFCToolTipCtrl), &ttParams);

	// Register the application's document templates.  Document templates
	//  serve as the connection between documents, frame windows and views
	CSingleDocTemplate* pDocTemplate;
	pDocTemplate = new CSingleDocTemplate(
		IDR_MAINFRAME,
		RUNTIME_CLASS(CParticleTrackingDoc),
		RUNTIME_CLASS(CMainFrame),       // main SDI frame window
		RUNTIME_CLASS(CParticleTrackingView));
	if (!pDocTemplate)
		return FALSE;
	AddDocTemplate(pDocTemplate);


	// Enable DDE Execute open
	EnableShellOpen();
	RegisterShellFileTypes(TRUE);

	// Parse command line for standard shell commands, DDE, file open
	CCommandLineInfo cmdInfo;
	ParseCommandLine(cmdInfo);


	// Dispatch commands specified on the command line.  Will return FALSE if
	// app was launched with /RegServer, /Register, /Unregserver or /Unregister.
	if (!ProcessShellCommand(cmdInfo))
		return FALSE;

	// The one and only window has been initialized, so show and update it
	m_pMainWnd->ShowWindow(SW_SHOW);
	m_pMainWnd->UpdateWindow();
	// call DragAcceptFiles only if there's a suffix
	//  In an SDI app, this should occur after ProcessShellCommand
	// Enable drag/drop open
	m_pMainWnd->DragAcceptFiles();
	return TRUE;
}

BOOL CParticleTrackingApp::OnIdle(LONG lCount)  // [MS] 28-05-2013.
{
  if(lCount == 0)
  {
    CMainFrame* pMainWnd = (CMainFrame*)m_pMainWnd;
    pMainWnd->OnIdle();

    if(GetAsyncKeyState(VK_ESCAPE))
      m_Tracker.terminate();
  }

  return CWinAppEx::OnIdle(lCount);
}

BOOL CParticleTrackingApp::PreTranslateMessage(MSG* pMsg)
{
  BOOL bResult = FALSE;

  if(pMsg->message == WM_KEYDOWN && pMsg->wParam  == VK_ESCAPE)
  {
    m_Tracker.terminate();
    m_Exporter.terminate();
    bResult = TRUE;
  }
  else if(pMsg->message == WM_KEYDOWN && pMsg->wParam  == VK_CONTROL)
  {
    m_Drawer.set_ctrl_pressed(true);
  }
  else if(pMsg->message == WM_KEYUP && pMsg->wParam  == VK_CONTROL)
  {
    m_Drawer.set_ctrl_pressed(false);
  }
  else
  {
    bResult = CWinAppEx::PreTranslateMessage(pMsg);
  }

  return bResult;
}

// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()

// App command to run the dialog
void CParticleTrackingApp::OnAppAbout()
{
	CAboutDlg aboutDlg;
	aboutDlg.DoModal();
}

// CParticleTrackingApp customization load/save methods

void CParticleTrackingApp::PreLoadState()
{
	BOOL bNameValid;
	CString strName;
	bNameValid = strName.LoadString(IDS_EDIT_MENU);
	ASSERT(bNameValid);
	GetContextMenuManager()->AddMenu(strName, IDR_POPUP_EDIT);
}

void CParticleTrackingApp::LoadCustomState()
{
}

void CParticleTrackingApp::SaveCustomState()
{
}

CParticleTrackingApp* CParticleTrackingApp::Get()
{
  return &theApp;
}

EvaporatingParticle::CTracker* CParticleTrackingApp::GetTracker()
{
  return &m_Tracker;
}

EvaporatingParticle::CDirichletTesselation* CParticleTrackingApp::GetDirichletTess()
{
  return &m_DirichletTess;
}

EvaporatingParticle::CTrackDraw* CParticleTrackingApp::GetDrawObj()
{
  return &m_Drawer;
}

EvaporatingParticle::CExportOpenFOAM* CParticleTrackingApp::GetExporter()
{
  return &m_Exporter;
}

EvaporatingParticle::CCalcCollection* CParticleTrackingApp::GetCalcs()
{
  return &m_vCalcs;
}

EvaporatingParticle::CFieldDataColl* CParticleTrackingApp::GetFields()
{
  return &m_vFields;
}

EvaporatingParticle::CCrossSectColl* CParticleTrackingApp::GetPlanes()
{
  return &m_vPlanes;
}

EvaporatingParticle::CSelAreasColl* CParticleTrackingApp::GetSelAreas()
{
  return &m_vSelAreas;
}

void CParticleTrackingApp::SelectedRegionChanged(EvaporatingParticle::CNamesVector* pRegNames)
{
  if(pRegNames == (EvaporatingParticle::CNamesVector*)(m_Tracker.get_src()->get_selected_rgn_names_ptr()))
  {
    m_Tracker.get_src()->invalidate();
    return;
  }

  if(m_vCalcs.sel_region_changed(pRegNames))
    return;

  if(m_vFields.sel_region_changed(pRegNames))
    return;
}

// CParticleTrackingApp message handlers



