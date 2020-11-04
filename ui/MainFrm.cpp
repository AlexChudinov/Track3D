
// MainFrm.cpp : implementation of the CMainFrame class
//

#include "stdafx.h"
#include "ParticleTracking.h"

#include "MainFrm.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CMainFrame

IMPLEMENT_DYNCREATE(CMainFrame, CFrameWndEx)

const int  iMaxUserToolbars = 10;
const UINT uiFirstUserToolBarId = AFX_IDW_CONTROLBAR_FIRST + 40;
const UINT uiLastUserToolBarId = uiFirstUserToolBarId + iMaxUserToolbars - 1;

BEGIN_MESSAGE_MAP(CMainFrame, CFrameWndEx)
	ON_WM_CREATE()
	ON_COMMAND(ID_VIEW_CUSTOMIZE, &CMainFrame::OnViewCustomize)
	ON_REGISTERED_MESSAGE(AFX_WM_CREATETOOLBAR, &CMainFrame::OnToolbarCreateNew)
	ON_COMMAND_RANGE(ID_VIEW_APPLOOK_WIN_2000, ID_VIEW_APPLOOK_OFF_2007_AQUA, &CMainFrame::OnApplicationLook)
	ON_UPDATE_COMMAND_UI_RANGE(ID_VIEW_APPLOOK_WIN_2000, ID_VIEW_APPLOOK_OFF_2007_AQUA, &CMainFrame::OnUpdateApplicationLook)
  ON_UPDATE_COMMAND_UI(ID_INDICATOR_REG, &CMainFrame::OnUpdatePane)
  ON_UPDATE_COMMAND_UI(ID_INDICATOR_X, &CMainFrame::OnUpdatePane)  // [MS] 31-05-2013 to support the status bar update.
  ON_UPDATE_COMMAND_UI(ID_INDICATOR_Y, &CMainFrame::OnUpdatePane)  //
  ON_UPDATE_COMMAND_UI(ID_INDICATOR_Z, &CMainFrame::OnUpdatePane)  //
  ON_COMMAND(ID_MOVE, &CMainFrame::OnButtonMove)
  ON_UPDATE_COMMAND_UI(ID_MOVE, &CMainFrame::OnUpdateButtonMove)
  ON_COMMAND(ID_ROTATE_X, &CMainFrame::OnButtonRotateX)
  ON_UPDATE_COMMAND_UI(ID_ROTATE_X, &CMainFrame::OnUpdateButtonRotateX)
  ON_COMMAND(ID_ROTATE_Y, &CMainFrame::OnButtonRotateY)
  ON_UPDATE_COMMAND_UI(ID_ROTATE_Y, &CMainFrame::OnUpdateButtonRotateY)
  ON_COMMAND(ID_ROTATE, &CMainFrame::OnButtonRotate)
  ON_UPDATE_COMMAND_UI(ID_ROTATE, &CMainFrame::OnUpdateButtonRotate)
END_MESSAGE_MAP()

static UINT indicators[] =
{
	ID_SEPARATOR,           // status line indicator
  ID_INDICATOR_REG,
  ID_INDICATOR_X,
  ID_INDICATOR_Y,
  ID_INDICATOR_Z,
	ID_INDICATOR_CAPS,
	ID_INDICATOR_NUM,
	ID_INDICATOR_SCRL,
};

// CMainFrame construction/destruction

CMainFrame::CMainFrame()
{
	// TODO: add member initialization code here
	theApp.m_nAppLook = theApp.GetInt(_T("ApplicationLook"), ID_VIEW_APPLOOK_VS_2005);
}

CMainFrame::~CMainFrame()
{
}

int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CFrameWndEx::OnCreate(lpCreateStruct) == -1)
		return -1;

	BOOL bNameValid;
	// set the visual manager and style based on persisted value
	OnApplicationLook(theApp.m_nAppLook);

	if (!m_wndMenuBar.Create(this))
	{
		TRACE0("Failed to create menubar\n");
		return -1;      // fail to create
	}

	m_wndMenuBar.SetPaneStyle(m_wndMenuBar.GetPaneStyle() | CBRS_SIZE_DYNAMIC | CBRS_TOOLTIPS | CBRS_FLYBY);

	// prevent the menu bar from taking the focus on activation
	CMFCPopupMenu::SetForceMenuFocus(FALSE);

	if (!m_wndToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_TOP | CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC) ||
		!m_wndToolBar.LoadToolBar(/*theApp.m_bHiColorIcons ? IDR_MAINFRAME_256 :*/ IDR_MAINFRAME))
	{
		TRACE0("Failed to create toolbar\n");
		return -1;      // fail to create
	}

	CString strToolBarName;
	bNameValid = strToolBarName.LoadString(IDS_TOOLBAR_STANDARD);
	ASSERT(bNameValid);
	m_wndToolBar.SetWindowText(strToolBarName);

	CString strCustomize;
	bNameValid = strCustomize.LoadString(IDS_TOOLBAR_CUSTOMIZE);
	ASSERT(bNameValid);
	m_wndToolBar.EnableCustomizeButton(TRUE, ID_VIEW_CUSTOMIZE, strCustomize);

  CRect rectClient;
	GetClientRect(rectClient);
	CSize size1 = m_wndToolBar.CalcFixedLayout(FALSE, TRUE);
  m_wndToolBar.SetWindowPos(NULL, rectClient.left, rectClient.top, size1.cx, size1.cy, SWP_NOACTIVATE | SWP_NOZORDER);

// [MS] 13-08-2014
  if(!m_wndMoveRotToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_TOP | CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC) ||
		 !m_wndMoveRotToolBar.LoadToolBar(IDR_TOOLBAR_MOVE_ROTATE))
	{
		TRACE0("Failed to create toolbar\n");
		return -1;      // fail to create
	}

  CSize size2 = m_wndMoveRotToolBar.CalcFixedLayout(FALSE, TRUE);
  m_wndMoveRotToolBar.SetWindowPos(NULL, rectClient.left + size1.cx, rectClient.top, size2.cx, size2.cy, SWP_NOACTIVATE | SWP_NOZORDER);

  CString strMoveRotToolbarName;
	bNameValid = strMoveRotToolbarName.LoadString(IDS_MOVE_ROTATE);
	ASSERT(bNameValid);
	m_wndMoveRotToolBar.SetWindowText(strMoveRotToolbarName);
// [/MS]

	// Allow user-defined toolbars operations:
  InitUserToolbars(NULL, uiFirstUserToolBarId, uiLastUserToolBarId);

	if (!m_wndStatusBar.Create(this))
	{
		TRACE0("Failed to create status bar\n");
		return -1;      // fail to create
	}
	m_wndStatusBar.SetIndicators(indicators, sizeof(indicators)/sizeof(UINT));
  m_wndStatusBar.SetPaneStyle(1, SBPS_NORMAL | SBPS_NOBORDERS);
  m_wndStatusBar.SetPaneWidth(1, 150);  // the first pane is reserved for the region name.

  m_wndStatusBar.SetPaneStyle(2, SBPS_NORMAL | SBPS_NOBORDERS);
  m_wndStatusBar.SetPaneStyle(3, SBPS_NORMAL | SBPS_NOBORDERS);
  m_wndStatusBar.SetPaneStyle(4, SBPS_NORMAL | SBPS_NOBORDERS);

	// TODO: Delete these five lines if you don't want the toolbar and menubar to be dockable
  m_wndMenuBar.EnableDocking(CBRS_ALIGN_ANY);
  m_wndToolBar.EnableDocking(CBRS_ALIGN_ANY);
  m_wndMoveRotToolBar.EnableDocking(CBRS_ALIGN_ANY);  // [MS] 13-08-2014
  EnableDocking(CBRS_ALIGN_ANY);
  DockPane(&m_wndMenuBar);
  DockPane(&m_wndToolBar);
  DockPane(&m_wndMoveRotToolBar); // [MS] 13-08-2014


	// enable Visual Studio 2005 style docking window behavior
  CDockingManager::SetDockingMode(DT_SMART);
	// enable Visual Studio 2005 style docking window auto-hide behavior
  EnableAutoHidePanes(CBRS_ALIGN_ANY);

	// create docking windows
	if (!CreateDockingWindows())
	{
		TRACE0("Failed to create docking windows\n");
		return -1;
	}

	m_wndProperties.EnableDocking(CBRS_ALIGN_ANY);
	DockPane(&m_wndProperties);


	// Enable toolbar and docking window menu replacement
  EnablePaneMenu(TRUE, ID_VIEW_CUSTOMIZE, strCustomize, ID_VIEW_TOOLBAR);

	// enable quick (Alt+drag) toolbar customization
  CMFCToolBar::EnableQuickCustomization();

	if (CMFCToolBar::GetUserImages() == NULL)
	{
		// load user-defined toolbar images
		if (m_UserImages.Load(_T(".\\UserImages.bmp")))
		{
			m_UserImages.SetImageSize(CSize(16, 16), FALSE);
			CMFCToolBar::SetUserImages(&m_UserImages);
		}
	}

	// enable menu personalization (most-recently used commands)
	// TODO: define your own basic commands, ensuring that each pulldown menu has at least one basic command.
	CList<UINT, UINT> lstBasicCommands;

	lstBasicCommands.AddTail(ID_FILE_NEW);
	lstBasicCommands.AddTail(ID_FILE_OPEN);
	lstBasicCommands.AddTail(ID_FILE_SAVE);
	lstBasicCommands.AddTail(ID_FILE_PRINT);
	lstBasicCommands.AddTail(ID_APP_EXIT);
	lstBasicCommands.AddTail(ID_EDIT_CUT);
	lstBasicCommands.AddTail(ID_EDIT_PASTE);
	lstBasicCommands.AddTail(ID_EDIT_UNDO);
	lstBasicCommands.AddTail(ID_APP_ABOUT);
	lstBasicCommands.AddTail(ID_VIEW_STATUS_BAR);
	lstBasicCommands.AddTail(ID_VIEW_TOOLBAR);
	lstBasicCommands.AddTail(ID_VIEW_APPLOOK_OFF_2003);
	lstBasicCommands.AddTail(ID_VIEW_APPLOOK_VS_2005);
	lstBasicCommands.AddTail(ID_VIEW_APPLOOK_OFF_2007_BLUE);
	lstBasicCommands.AddTail(ID_VIEW_APPLOOK_OFF_2007_SILVER);
	lstBasicCommands.AddTail(ID_VIEW_APPLOOK_OFF_2007_BLACK);
	lstBasicCommands.AddTail(ID_VIEW_APPLOOK_OFF_2007_AQUA);

	CMFCToolBar::SetBasicCommands(lstBasicCommands);

	return 0;
}

BOOL CMainFrame::OnCreateClient(LPCREATESTRUCT /*lpcs*/,
	CCreateContext* pContext)
{
	return m_wndSplitter.Create(this,
		2, 2,               // TODO: adjust the number of rows, columns
		CSize(10, 10),      // TODO: adjust the minimum pane size
		pContext);
}

BOOL CMainFrame::PreCreateWindow(CREATESTRUCT& cs)
{
	if( !CFrameWndEx::PreCreateWindow(cs) )
		return FALSE;
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return TRUE;
}

BOOL CMainFrame::CreateDockingWindows()
{
	BOOL bNameValid;
	// Create properties window
	CString strPropertiesWnd;
	bNameValid = strPropertiesWnd.LoadString(IDS_PROPERTIES_WND);
	ASSERT(bNameValid);
	if (!m_wndProperties.Create(strPropertiesWnd, this, CRect(0, 0, 200, 200), TRUE, ID_VIEW_PROPERTIESWND, WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_CLIPCHILDREN | CBRS_RIGHT /*| CBRS_FLOAT_MULTI*/))
	{
		TRACE0("Failed to create Properties window\n");
		return FALSE; // failed to create
	}

	SetDockingWindowIcons(theApp.m_bHiColorIcons);
	return TRUE;
}

void CMainFrame::SetDockingWindowIcons(BOOL bHiColorIcons)
{
	HICON hPropertiesBarIcon = (HICON) ::LoadImage(::AfxGetResourceHandle(), MAKEINTRESOURCE(bHiColorIcons ? IDI_PROPERTIES_WND_HC : IDI_PROPERTIES_WND), IMAGE_ICON, ::GetSystemMetrics(SM_CXSMICON), ::GetSystemMetrics(SM_CYSMICON), 0);
	m_wndProperties.SetIcon(hPropertiesBarIcon, FALSE);

}

// CMainFrame diagnostics

#ifdef _DEBUG
void CMainFrame::AssertValid() const
{
	CFrameWndEx::AssertValid();
}

void CMainFrame::Dump(CDumpContext& dc) const
{
	CFrameWndEx::Dump(dc);
}
#endif //_DEBUG


// CMainFrame message handlers

void CMainFrame::OnViewCustomize()
{
	CMFCToolBarsCustomizeDialog* pDlgCust = new CMFCToolBarsCustomizeDialog(this, TRUE /* scan menus */);
	pDlgCust->EnableUserDefinedToolbars();
	pDlgCust->Create();
}

LRESULT CMainFrame::OnToolbarCreateNew(WPARAM wp,LPARAM lp)
{
	LRESULT lres = CFrameWndEx::OnToolbarCreateNew(wp,lp);
	if (lres == 0)
	{
		return 0;
	}

	CMFCToolBar* pUserToolbar = (CMFCToolBar*)lres;
	ASSERT_VALID(pUserToolbar);

	BOOL bNameValid;
	CString strCustomize;
	bNameValid = strCustomize.LoadString(IDS_TOOLBAR_CUSTOMIZE);
	ASSERT(bNameValid);

	pUserToolbar->EnableCustomizeButton(TRUE, ID_VIEW_CUSTOMIZE, strCustomize);
	return lres;
}

void CMainFrame::OnApplicationLook(UINT id)
{
	CWaitCursor wait;

	theApp.m_nAppLook = id;

	switch (theApp.m_nAppLook)
	{
	case ID_VIEW_APPLOOK_WIN_2000:
		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManager));
		break;

	case ID_VIEW_APPLOOK_OFF_XP:
		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerOfficeXP));
		break;

	case ID_VIEW_APPLOOK_WIN_XP:
		CMFCVisualManagerWindows::m_b3DTabsXPTheme = TRUE;
		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerWindows));
		break;

	case ID_VIEW_APPLOOK_OFF_2003:
		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerOffice2003));
		CDockingManager::SetDockingMode(DT_SMART);
		break;

	case ID_VIEW_APPLOOK_VS_2005:
		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerVS2005));
		CDockingManager::SetDockingMode(DT_SMART);
		break;

	default:
		switch (theApp.m_nAppLook)
		{
		case ID_VIEW_APPLOOK_OFF_2007_BLUE:
			CMFCVisualManagerOffice2007::SetStyle(CMFCVisualManagerOffice2007::Office2007_LunaBlue);
			break;

		case ID_VIEW_APPLOOK_OFF_2007_BLACK:
			CMFCVisualManagerOffice2007::SetStyle(CMFCVisualManagerOffice2007::Office2007_ObsidianBlack);
			break;

		case ID_VIEW_APPLOOK_OFF_2007_SILVER:
			CMFCVisualManagerOffice2007::SetStyle(CMFCVisualManagerOffice2007::Office2007_Silver);
			break;

		case ID_VIEW_APPLOOK_OFF_2007_AQUA:
			CMFCVisualManagerOffice2007::SetStyle(CMFCVisualManagerOffice2007::Office2007_Aqua);
			break;
		}

		CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerOffice2007));
		CDockingManager::SetDockingMode(DT_SMART);
	}

	RedrawWindow(NULL, NULL, RDW_ALLCHILDREN | RDW_INVALIDATE | RDW_UPDATENOW | RDW_FRAME | RDW_ERASE);

	theApp.WriteInt(_T("ApplicationLook"), theApp.m_nAppLook);
}

void CMainFrame::OnUpdateApplicationLook(CCmdUI* pCmdUI)
{
	pCmdUI->SetRadio(theApp.m_nAppLook == pCmdUI->m_nID);
}

BOOL CMainFrame::LoadFrame(UINT nIDResource, DWORD dwDefaultStyle, CWnd* pParentWnd, CCreateContext* pContext) 
{
	// base class does the real work

	if (!CFrameWndEx::LoadFrame(nIDResource, dwDefaultStyle, pParentWnd, pContext))
	{
		return FALSE;
	}


	// enable customization button for all user toolbars
	BOOL bNameValid;
	CString strCustomize;
	bNameValid = strCustomize.LoadString(IDS_TOOLBAR_CUSTOMIZE);
	ASSERT(bNameValid);

	for (int i = 0; i < iMaxUserToolbars; i ++)
	{
		CMFCToolBar* pUserToolbar = GetUserToolBarByIndex(i);
		if (pUserToolbar != NULL)
		{
			pUserToolbar->EnableCustomizeButton(TRUE, ID_VIEW_CUSTOMIZE, strCustomize);
		}
	}

	return TRUE;
}

void CMainFrame::OnIdle() // [MS] 28-05-2013.
{
  m_wndProperties.on_idle();
}

void CMainFrame::OnUpdatePane(CCmdUI* pCmdUI) // [MS] 31-05-2013 to support the status bar update.
{
  pCmdUI->Enable();
}

void CMainFrame::OnButtonMove()
{
  EvaporatingParticle::CTrackDraw* pObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pObj == NULL)
    return;

  pObj->set_context(EvaporatingParticle::CTrackDraw::nContextMove);
}

void CMainFrame::OnUpdateButtonMove(CCmdUI* pCmdUI)
{
  EvaporatingParticle::CTrackDraw* pObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pObj == NULL)
    return;

  pCmdUI->SetCheck(pObj->get_context() == EvaporatingParticle::CTrackDraw::nContextMove);
}

void CMainFrame::OnButtonRotateX()
{
  EvaporatingParticle::CTrackDraw* pObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pObj == NULL)
    return;

  pObj->set_context(EvaporatingParticle::CTrackDraw::nContextRotX);
}

void CMainFrame::OnUpdateButtonRotateX(CCmdUI* pCmdUI)
{
  EvaporatingParticle::CTrackDraw* pObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pObj == NULL)
    return;

  pCmdUI->SetCheck(pObj->get_context() == EvaporatingParticle::CTrackDraw::nContextRotX);
}

void CMainFrame::OnButtonRotateY()
{
  EvaporatingParticle::CTrackDraw* pObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pObj == NULL)
    return;

  pObj->set_context(EvaporatingParticle::CTrackDraw::nContextRotY);
}

void CMainFrame::OnUpdateButtonRotateY(CCmdUI* pCmdUI)
{
  EvaporatingParticle::CTrackDraw* pObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pObj == NULL)
    return;

  pCmdUI->SetCheck(pObj->get_context() == EvaporatingParticle::CTrackDraw::nContextRotY);
}

void CMainFrame::OnButtonRotate()
{
  EvaporatingParticle::CTrackDraw* pObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pObj == NULL)
    return;

  pObj->set_context(EvaporatingParticle::CTrackDraw::nContextRotZ);
}

void CMainFrame::OnUpdateButtonRotate(CCmdUI* pCmdUI)
{
  EvaporatingParticle::CTrackDraw* pObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pObj == NULL)
    return;

  pCmdUI->SetCheck(pObj->get_context() == EvaporatingParticle::CTrackDraw::nContextRotZ);
}

