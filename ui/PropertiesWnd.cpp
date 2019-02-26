
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "Resource.h"
#include "MainFrm.h"
#include "ParticleTracking.h"
#include "ParticleTrackingDoc.h"

#include "EvaporationModel.h"
#include "ExportOpenFOAM.h"
#include "ColorContour.h"

#include "ExecutionDialog.h"
#include "AddCalculatorDialog.h"
#include "AddFieldPtbDialog.h"

#include "ResponseProperty.h"
#include "Button.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

/////////////////////////////////////////////////////////////////////////////
// CResourceViewBar

CPropertiesWnd::CPropertiesWnd()
{
  m_bUpdateAll = false;
  m_bIgnoreIdle = false;
}

CPropertiesWnd::~CPropertiesWnd()
{
}

BEGIN_MESSAGE_MAP(CPropertiesWnd, CDockablePane)
	ON_WM_CREATE()
	ON_WM_SIZE()
	ON_COMMAND(ID_EXPAND_ALL, OnExpandAllProperties)
	ON_UPDATE_COMMAND_UI(ID_EXPAND_ALL, OnUpdateExpandAllProperties)
	ON_COMMAND(ID_SORTPROPERTIES, OnSortProperties)
	ON_UPDATE_COMMAND_UI(ID_SORTPROPERTIES, OnUpdateSortProperties)
	ON_COMMAND(ID_PROPERTIES1, OnDoTracking)
	ON_UPDATE_COMMAND_UI(ID_PROPERTIES1, OnUpdateDoTracking)
	ON_COMMAND(ID_PROPERTIES2, OnSaveImage)
	ON_UPDATE_COMMAND_UI(ID_PROPERTIES2, OnUpdateSaveImage)
  ON_COMMAND(ID_PROPERTIES3, OnExportMesh)
	ON_UPDATE_COMMAND_UI(ID_PROPERTIES3, OnUpdateExportMesh)
  ON_COMMAND(ID_PROPERTIES4, OnCreateBoundaryConditions)
	ON_UPDATE_COMMAND_UI(ID_PROPERTIES4, OnUpdateCreateBoundaryConditions)
  ON_COMMAND(ID_PROPERTIES5, OnImportOpenFOAM)
	ON_UPDATE_COMMAND_UI(ID_PROPERTIES5, OnUpdateImportOpenFOAM)
  ON_COMMAND(ID_PROPERTIES6, OnCreateContour)
	ON_UPDATE_COMMAND_UI(ID_PROPERTIES6, OnUpdateCreateContour)
  ON_COMMAND(ID_PROPERTIES7, OnCreateCalculator)
	ON_UPDATE_COMMAND_UI(ID_PROPERTIES7, OnUpdateCreateCalculator)
  ON_COMMAND(ID_PROPERTIES8, OnAddFieldPeturbation)
	ON_UPDATE_COMMAND_UI(ID_PROPERTIES8, OnUpdateAddFieldPeturbation)
  ON_COMMAND(ID_PROPERTIES9, OnAddField)
	ON_UPDATE_COMMAND_UI(ID_PROPERTIES9, OnUpdateAddField)
  ON_COMMAND(ID_PROPERTIES10, OnAddCrossSection)
	ON_UPDATE_COMMAND_UI(ID_PROPERTIES10, OnUpdateAddCrossSection)
	ON_WM_SETFOCUS()
	ON_WM_SETTINGCHANGE()
  ON_NOTIFY(TCN_SELCHANGING, ID_TAB_HOLDER, OnTabSelChanging)
  ON_NOTIFY(TCN_SELCHANGE, ID_TAB_HOLDER, OnTabSelChange)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CResourceViewBar message handlers

void CPropertiesWnd::AdjustLayout()
{
  if(GetSafeHwnd() == NULL)
    return;

  CRect rectClient, rectTab;
  GetClientRect(rectClient);

	int cyTlb = m_wndToolBar.CalcFixedLayout(FALSE, TRUE).cy;
	m_wndToolBar.SetWindowPos(NULL, rectClient.left, rectClient.top, rectClient.Width(), cyTlb, SWP_NOACTIVATE | SWP_NOZORDER);

  m_wndTabCtrl.SetWindowPos(NULL, rectClient.left, rectClient.top + cyTlb, rectClient.Width(), rectClient.Height() - cyTlb, SWP_NOACTIVATE | SWP_NOZORDER);
  m_wndTabCtrl.AdjustRect(TRUE, &rectTab);
  int cyTab = rectTab.Height();

	m_wndPropList.SetWindowPos(&m_wndTabCtrl, rectClient.left, rectClient.top + cyTlb + cyTab, rectClient.Width(), rectClient.Height() - (cyTab + cyTlb), SWP_NOACTIVATE | SWP_NOZORDER);
}

static const char* scTabItemNames[CPropertiesWnd::nTabCount] = 
{ "Main", "Source", "Tracking", "Droplet Param", "Ion Param", "Drawing", "Export", "Import", "Calculators", "Perturbations", "Field" };

int CPropertiesWnd::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CDockablePane::OnCreate(lpCreateStruct) == -1)
		return -1;

	CRect rectDummy;
	rectDummy.SetRectEmpty();

// Create the tab control:
  if(!m_wndTabCtrl.Create(WS_VISIBLE | WS_CHILD | WS_CLIPSIBLINGS | TCS_TABS  | TCS_MULTILINE  | TCS_FOCUSNEVER  | TCS_HOTTRACK  | TCS_TOOLTIPS, rectDummy, this, ID_TAB_HOLDER))
  {
    TRACE0("Failed to create CTabCtrl window\n");
		return FALSE; // failed to create
  }

  long nRes = -1;
  for(UINT i = 0; i < nTabCount; i++)
    nRes = m_wndTabCtrl.InsertItem(i, scTabItemNames[i]);

  m_wndTabCtrl.SetFont(CFont::FromHandle((HFONT)GetStockObject(DEFAULT_GUI_FONT)));
  m_wndTabCtrl.SetCurSel(0);
  m_wndTabCtrl.HighlightItem(0, TRUE);

// Create the properties list:
	if(!m_wndPropList.Create(WS_VISIBLE | WS_CHILD, rectDummy, this, 1))
	{
		TRACE0("Failed to create Properties Grid \n");
		return -1;      // fail to create
	}

  m_bBusy = false;
	InitPropList();

// Create the property toolbar:
	m_wndToolBar.Create(this, AFX_DEFAULT_TOOLBAR_STYLE, IDR_PROPERTIES);
	m_wndToolBar.LoadToolBar(IDR_PROPERTIES, 0, 0, TRUE /* Is locked */);
	m_wndToolBar.CleanUpLockedImages();
	m_wndToolBar.LoadBitmap(/*theApp.m_bHiColorIcons ? IDB_PROPERTIES_HC :*/ IDR_PROPERTIES, 0, 0, TRUE /* Locked */);

	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() | CBRS_TOOLTIPS | CBRS_FLYBY);
	m_wndToolBar.SetPaneStyle(m_wndToolBar.GetPaneStyle() & ~(CBRS_GRIPPER | CBRS_SIZE_DYNAMIC | CBRS_BORDER_TOP | CBRS_BORDER_BOTTOM | CBRS_BORDER_LEFT | CBRS_BORDER_RIGHT));
	m_wndToolBar.SetOwner(this);

// All commands will be routed via this control , not via the parent frame:
	m_wndToolBar.SetRouteCommandsViaFrame(FALSE);

	AdjustLayout();
	return 0;
}

void CPropertiesWnd::OnSize(UINT nType, int cx, int cy)
{
	CDockablePane::OnSize(nType, cx, cy);
	AdjustLayout();
}

void CPropertiesWnd::OnExpandAllProperties()
{
	m_wndPropList.ExpandAll();
}

void CPropertiesWnd::OnUpdateExpandAllProperties(CCmdUI* pCmdUI)
{
}

void CPropertiesWnd::OnSortProperties()
{
	m_wndPropList.SetAlphabeticMode(!m_wndPropList.IsAlphabeticMode());
}

void CPropertiesWnd::OnUpdateSortProperties(CCmdUI* pCmdUI)
{
	pCmdUI->SetCheck(m_wndPropList.IsAlphabeticMode());
}

static UINT __stdcall main_thread_func(LPVOID pData)
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  if(pObj->get_enable_ansys_field())
  {
    if(pObj->is_ready())
    {
      pObj->read_gasdyn_data(true); // restore only ANSYS-calculated electric fields, the other node data will be unchanged.
      pObj->apply_perturbations();  // field perturbations are applied and stored in the mesh nodes.
    }
  }
  else
  {
    if(!pObj->is_ready() && !pObj->read_data())
    {
      CDialog* pDlg = (CDialog*)pData;
      pDlg->PostMessage(WM_CLOSE);
      return 0;
    }

    EvaporatingParticle::CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
    if(!pFields->calc_fields())
    {
      CDialog* pDlg = (CDialog*)pData;
      pDlg->PostMessage(WM_CLOSE);
      return 0;
    }
  }

  if(pObj->get_use_multi_thread())
    pObj->multi_thread_calculate();
  else
    pObj->single_thread_calculate();

  if(!pObj->get_terminate_flag()) // if the user has not terminated the dialog manually, do it after calculations have ended.
  {
    CDialog* pDlg = (CDialog*)pData;
    pDlg->PostMessage(WM_CLOSE);
  }

  return 0;
}

void CPropertiesWnd::OnDoTracking()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  pObj->terminate(false);

  set_data_to_model();
 
  CExecutionDialog dlg(&main_thread_func, (EvaporatingParticle::CObject*)pObj);
  INT_PTR nRes = dlg.DoModal();

  if(pObj->get_result_flag())
  {
    CParticleTrackingApp::Get()->GetCalcs()->invalidate_calcs();
    EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
    pDrawObj->invalidate_tracks();
    pDrawObj->draw();
  }

  set_update_all();
}

void CPropertiesWnd::OnUpdateDoTracking(CCmdUI* pCmdUI)
{
	pCmdUI->Enable(!m_bBusy);
}

void CPropertiesWnd::OnAddFieldPeturbation()
{
  set_data_to_model();
  
  CAddFieldPtbDialog dlg;
  int nRes = dlg.DoModal();
  if(nRes == IDOK)
    set_update_all();
}

void CPropertiesWnd::OnUpdateAddFieldPeturbation(CCmdUI* pCmdUI)
{
   pCmdUI->Enable(!m_bBusy && CParticleTrackingApp::Get()->GetTracker()->is_ready());
}

void CPropertiesWnd::OnSaveImage()
{
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(!pDrawObj->capture_image())
    return;

	CString cFilter = "Bitmap image|*.bmp|JPEG image|*.jpg|GIF image|*.gif|PNG image|*.png||";
	CFileDialog dlg(FALSE, NULL, NULL, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT | OFN_EXPLORER, cFilter);
  dlg.m_ofn.nFilterIndex = 4;

	HRESULT hResult = (int)dlg.DoModal();
	if(FAILED(hResult))
		return;

	CString cFileName = dlg.m_ofn.lpstrFile, cExt;

// Add the file extension if the user didn't supply one
	if(dlg.m_ofn.nFileExtension == 0) 
		cFileName = cFileName + ".png";

  std::string cName((const char*)cFileName);
  bool bFilled = cName.size() != 0;
  if(bFilled)
    pDrawObj->save_image(cFileName);
}

void CPropertiesWnd::OnUpdateSaveImage(CCmdUI* pCmdUI)
{
	pCmdUI->Enable(!m_bBusy && CParticleTrackingApp::Get()->GetTracker()->is_ready());
}

static UINT __stdcall export_thread_func(LPVOID pData)
{
  EvaporatingParticle::CExportOpenFOAM* pExportObj = CParticleTrackingApp::Get()->GetExporter();

  pExportObj->do_export();

  if(!pExportObj->get_terminate_flag()) // if the user has not terminated the dialog manually, do it after calculations have ended.
  {
    CDialog* pDlg = (CDialog*)pData;
    pDlg->PostMessage(WM_CLOSE);
  }

  return 0;
}

static int CALLBACK BrowseCallbackProc(HWND hWnd, UINT nMsg, LPARAM lParam, LPARAM lpData)
{
  if(nMsg == BFFM_INITIALIZED)
    SendMessage(hWnd, BFFM_SETSELECTION, TRUE, lpData);

  return 0;
}

void CPropertiesWnd::OnExportMesh()
{
  EvaporatingParticle::CExportOpenFOAM* pExportObj = CParticleTrackingApp::Get()->GetExporter();

  BROWSEINFO bi = { 0 };
  bi.lpszTitle = "Select a folder";
  bi.ulFlags = BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE;
  bi.lpfn = BrowseCallbackProc;
  bi.lParam = (LPARAM)(pExportObj->get_path());

  LPMALLOC shMalloc = NULL;
  char folderpath[1024];

  CoInitialize(NULL);
  SHGetMalloc(&shMalloc);

  LPITEMIDLIST pList = SHBrowseForFolder(&bi);

  if(pList != NULL)
  {
    SHGetPathFromIDList(pList, folderpath);
    shMalloc->Free(pList);
    shMalloc->Release();
  
    CWaitCursor wait;
    m_bBusy = true;
    pExportObj->set_path((const char*)folderpath);
    
    set_export_data();
    
    CExecutionDialog dlg(&export_thread_func, (EvaporatingParticle::CObject*)pExportObj);
    INT_PTR nRes = dlg.DoModal();

    m_bBusy = false;
  }
}

void CPropertiesWnd::OnUpdateExportMesh(CCmdUI* pCmdUI)
{
  pCmdUI->Enable(!m_bBusy);
}

static UINT __stdcall import_thread_func(LPVOID pData)
{
  EvaporatingParticle::CImportOpenFOAM& ImpObj = CParticleTrackingApp::Get()->GetTracker()->get_importer();

  ImpObj.do_import();

  if(!ImpObj.get_terminate_flag()) // if the user has not terminated the dialog manually, do it after calculations have ended.
  {
    CDialog* pDlg = (CDialog*)pData;
    pDlg->PostMessage(WM_CLOSE);
  }

  return 0;
}

void CPropertiesWnd::OnImportOpenFOAM()
{
  CWaitCursor wait;
  m_bBusy = true;
  EvaporatingParticle::CImportOpenFOAM& ImpObj = CParticleTrackingApp::Get()->GetTracker()->get_importer();
  
  set_import_data();
  
  CExecutionDialog dlg(&import_thread_func, (EvaporatingParticle::CObject*)&ImpObj);
  INT_PTR nRes = dlg.DoModal();

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(!ImpObj.get_terminate_flag())
  {
    pDrawObj->invalidate_contours();
    pDrawObj->draw();
  }

  m_bBusy = false;
}

void CPropertiesWnd::OnUpdateImportOpenFOAM(CCmdUI* pCmdUI)
{
  pCmdUI->Enable(!m_bBusy && CParticleTrackingApp::Get()->GetTracker()->is_ready());
}

void CPropertiesWnd::OnCreateBoundaryConditions()
{
  set_data_to_model();
  EvaporatingParticle::CExportOpenFOAM* pExportObj = CParticleTrackingApp::Get()->GetExporter();
  pExportObj->add_bound_cond(new EvaporatingParticle::CBoundaryConditions());
  set_update_all();
}

void CPropertiesWnd::OnUpdateCreateBoundaryConditions(CCmdUI* pCmdUI)
{
  int nSelTab = m_wndTabCtrl.GetCurSel();
  bool bProperTab = nSelTab == tabExport;
  pCmdUI->Enable(!m_bBusy && CParticleTrackingApp::Get()->GetTracker()->is_ready() && bProperTab);
}

void CPropertiesWnd::OnCreateContour()
{
  set_data_to_model();
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->add_contour();
  set_update_all();
}

void CPropertiesWnd::OnUpdateCreateContour(CCmdUI* pCmdUI)
{
  int nSelTab = m_wndTabCtrl.GetCurSel();
  bool bProperTab = true; //nSelTab == tabDrawing;
  pCmdUI->Enable(!m_bBusy && CParticleTrackingApp::Get()->GetTracker()->is_ready() && bProperTab);
}

void CPropertiesWnd::OnCreateCalculator()
{
  set_data_to_model();
  
  CAddCalculatorDialog dlg;
  int nRes = dlg.DoModal();
  if(nRes == IDOK)
    set_update_all();
}

void CPropertiesWnd::OnUpdateCreateCalculator(CCmdUI* pCmdUI)
{
  int nSelTab = m_wndTabCtrl.GetCurSel();
  bool bProperTab = true; //nSelTab == tabCalc;
  pCmdUI->Enable(!m_bBusy && CParticleTrackingApp::Get()->GetTracker()->is_ready() && bProperTab);
}

void CPropertiesWnd::OnAddField()
{
  set_data_to_model();
  
  EvaporatingParticle::CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
  pFields->push_back(new EvaporatingParticle::CElectricFieldData());
  pFields->set_curr_field_index(pFields->size() - 1);
  set_update_all();
}

void CPropertiesWnd::OnUpdateAddField(CCmdUI* pCmdUI)
{
  int nSelTab = m_wndTabCtrl.GetCurSel();
  bool bProperTab = true; //nSelTab == tabFields;
  pCmdUI->Enable(!m_bBusy && CParticleTrackingApp::Get()->GetTracker()->is_ready() && bProperTab);
}

void CPropertiesWnd::OnAddCrossSection()
{
  EvaporatingParticle::CCrossSectColl* pPlanes = CParticleTrackingApp::Get()->GetPlanes();
  pPlanes->push_back(new EvaporatingParticle::CDomainCrossSection());
  set_update_all();

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->invalidate_faces();
  pDrawObj->invalidate_aux();
  pDrawObj->draw();
}

void CPropertiesWnd::OnUpdateAddCrossSection(CCmdUI* pCmdUI)
{
  pCmdUI->Enable(!m_bBusy && CParticleTrackingApp::Get()->GetTracker()->is_ready());
}

void CPropertiesWnd::OnSetFocus(CWnd* pOldWnd)
{
	CDockablePane::OnSetFocus(pOldWnd);
	m_wndPropList.SetFocus();
}

void CPropertiesWnd::OnSettingChange(UINT uFlags, LPCTSTR lpszSection)
{
	CDockablePane::OnSettingChange(uFlags, lpszSection);
	SetPropListFont();
}

void CPropertiesWnd::SetPropListFont()
{
	::DeleteObject(m_fntPropList.Detach());

	LOGFONT lf;
	afxGlobalData.fontRegular.GetLogFont(&lf);

	NONCLIENTMETRICS info;
	info.cbSize = sizeof(info);

	afxGlobalData.GetNonClientMetrics(info);

	lf.lfHeight = info.lfMenuFont.lfHeight;
	lf.lfWeight = info.lfMenuFont.lfWeight;
	lf.lfItalic = info.lfMenuFont.lfItalic;

	m_fntPropList.CreateFontIndirect(&lf);

	m_wndPropList.SetFont(&m_fntPropList);
}

void CPropertiesWnd::on_idle()
{
  if(m_bIgnoreIdle)
    return;

  if(m_bUpdateAll)
    update_prop_list();
  else
    update_ctrls();

  RedrawWindow();
}

void CPropertiesWnd::update_prop_list()
{
  int nPos = 0;
  CScrollBar* pScrBar = m_wndPropList.GetScrollBarCtrl(SB_VERT);
  if(pScrBar != NULL)
    nPos = pScrBar->GetScrollPos();

  m_wndPropList.RemoveAll();
  InitPropList();

  if(pScrBar != NULL)
  {
    pScrBar->SetScrollPos(nPos);
    WPARAM wParam = MAKEWPARAM(SB_THUMBPOSITION, nPos);
    LPARAM lParam = MAKELPARAM(pScrBar->m_hWnd, NULL);
    ::SendMessage(m_wndPropList.m_hWnd, WM_VSCROLL, wParam, lParam);
  }

  m_bUpdateAll = false;
}

void CPropertiesWnd::InitPropList()
{
	SetPropListFont();

	m_wndPropList.EnableHeaderCtrl(FALSE);
	m_wndPropList.EnableDescriptionArea();
	m_wndPropList.SetVSDotNetLook();
	m_wndPropList.MarkModifiedProperties();

  int nSelTab = m_wndTabCtrl.GetCurSel();
  switch(nSelTab)
  {
    case tabMain:     add_type_ctrls();     break;
    case tabSource:   add_source_ctrls();   break;
    case tabTrack:    add_tracking_ctrls(); break;
    case tabEvapor:   add_evapor_ctrls();   break;
    case tabIonPar:   add_ion_ctrls();      break;
    case tabDrawing:  add_draw_ctrls();     break;
    case tabExport:   add_export_ctrls();   break;
    case tabImport:   add_import_ctrls();   break;
    case tabCalc:     add_calc_ctrls();     break;
    case tabPtb:      add_ptb_ctrls();      break;
    case tabFields:   add_field_ctrls();    break;
  }
}

void CPropertiesWnd::set_data_to_model()
{
  int nSelTab = m_wndTabCtrl.GetCurSel();
  switch(nSelTab)
  {
    case tabMain:     set_type_data();      break;
    case tabSource:   set_source_data();    break;
    case tabTrack:    set_tracking_data();  break;
    case tabEvapor:   set_evapor_data();    break;
    case tabIonPar:   set_ion_data();       break;
    case tabDrawing:  set_draw_data();      break;
    case tabExport:   set_export_data();    break;
    case tabImport:   set_import_data();    break;
    case tabCalc:     set_calc_data();      break;
    case tabPtb:      set_ptb_data();       break;
    case tabFields:   set_field_data();     break;
  }
}

void CPropertiesWnd::update_ctrls()
{
  int nSelTab = m_wndTabCtrl.GetCurSel();
  switch(nSelTab)
  {
    case tabMain:     update_type_ctrls();      break;
    case tabSource:   update_source_ctrls();    break;
    case tabTrack:    update_tracking_ctrls();  break;
    case tabEvapor:   update_evapor_ctrls();    break;
    case tabIonPar:   update_ion_ctrls();       break;
    case tabDrawing:  update_draw_ctrls();      break;
    case tabExport:   update_export_ctrls();    break;
    case tabImport:   update_import_ctrls();    break;
    case tabCalc:     update_calc_ctrls();      break;
    case tabPtb:      update_ptb_ctrls();       break;
    case tabFields:   update_field_ctrls();     break;
  }
}

void CPropertiesWnd::disable_all_but_one(CMFCPropertyGridProperty* pBtn, CMFCPropertyGridProperty* pParent)
{
  if(pParent != NULL)
  {
    if(pParent != pBtn)
      pParent->Enable(FALSE);

    for(int j = 0; j < pParent->GetSubItemsCount(); j++)
      disable_all_but_one(pBtn, pParent->GetSubItem(j));
  }
  else
  {
    int nPropCount = m_wndPropList.GetPropertyCount();
    for(int i = 0; i < nPropCount; i++)
    {
      CMFCPropertyGridProperty* pProp = m_wndPropList.GetProperty(i);
      if(pProp == pBtn)
        continue;

      disable_all_but_one(pBtn, pProp);
    }
  }
}

//---------------------------------------------------------------------------------------
//  Notifications handlers
//---------------------------------------------------------------------------------------
LRESULT CPropertiesWnd::WindowProc(UINT message, WPARAM wParam, LPARAM lParam)
{
  if(message == WM_NOTIFY)
  {
    if(((LPNMHDR)lParam)->code == TTN_GETDISPINFO)
    {
      LPTOOLTIPTEXT lpttt = (LPTOOLTIPTEXT)lParam; 
// Set the instance of the module that contains the resource.
      lpttt->hinst = CParticleTrackingApp::Get()->m_hInstance;
      UINT_PTR idButton = lpttt->hdr.idFrom;      
      switch(idButton)
      {
        case 1: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_0); break;
        case 2: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_1); break;
        case 3: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_2); break;
        case 4: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_3); break;
        case 5: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_7); break;
        case 6: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_8); break;
        case 7: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_9); break;
        case 8: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_10); break;
        case 9: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_11); break;
        case 10: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_12); break;
        case 11: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_13); break;
        case 12: lpttt->lpszText = MAKEINTRESOURCE(IDS_TOOL_TIP_14); break;
      }
    }
  }

  return CDockablePane::WindowProc(message, wParam, lParam);
}

void CPropertiesWnd::OnTabSelChanging(NMHDR* pNMHDR, LRESULT* pLResult)
{
  if((pNMHDR->hwndFrom == m_wndTabCtrl.m_hWnd) && (pNMHDR->code == TCN_SELCHANGING))
  {
    set_data_to_model();
    int nSel = m_wndTabCtrl.GetCurSel();
    m_wndTabCtrl.HighlightItem(nSel, FALSE);
  }
}

void CPropertiesWnd::OnTabSelChange(NMHDR* pNMHDR, LRESULT* pLResult)
{
  if((pNMHDR->hwndFrom == m_wndTabCtrl.m_hWnd) && (pNMHDR->code == TCN_SELCHANGE))
  {
    int nSel = m_wndTabCtrl.GetCurSel();
    m_wndTabCtrl.HighlightItem(nSel, TRUE);
    set_update_all();

    EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
    pDrawObj->hide_selected();
    pDrawObj->clear_selected_regions();
    pDrawObj->draw();
  }
}

