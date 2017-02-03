
// MainFrm.h : interface of the CMainFrame class
//

#pragma once
#include "PropertiesWnd.h"

class CMainFrame : public CFrameWndEx
{
	
protected: // create from serialization only
	CMainFrame();
	DECLARE_DYNCREATE(CMainFrame)

// Attributes
protected:
	CSplitterWnd m_wndSplitter;
public:

// Operations
public:
  void              OnIdle();             // [MS] 28-05-2013.
  CPropertiesWnd*   GetWndProperties();   // [MS] 30-05-2013.
  CMFCStatusBar*    GetStatusBar();       // [MS] 31-05-2013.

// Overrides
public:
	virtual BOOL OnCreateClient(LPCREATESTRUCT lpcs, CCreateContext* pContext);
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	virtual BOOL LoadFrame(UINT nIDResource, DWORD dwDefaultStyle = WS_OVERLAPPEDWINDOW | FWS_ADDTOTITLE, CWnd* pParentWnd = NULL, CCreateContext* pContext = NULL);

// Implementation
public:
	virtual ~CMainFrame();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:  // control bar embedded members
	CMFCMenuBar       m_wndMenuBar;
	CMFCToolBar       m_wndToolBar;
	CMFCStatusBar     m_wndStatusBar;
	CMFCToolBarImages m_UserImages;
	CPropertiesWnd    m_wndProperties;
  CMFCToolBar       m_wndMoveRotToolBar;  // [MS] 13-08-2014.

// Generated message map functions
protected:
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnViewCustomize();
	afx_msg LRESULT OnToolbarCreateNew(WPARAM wp, LPARAM lp);
	afx_msg void OnApplicationLook(UINT id);
	afx_msg void OnUpdateApplicationLook(CCmdUI* pCmdUI);
  afx_msg void OnUpdatePane(CCmdUI* pCmdUI);
  afx_msg void OnButtonMove();
  afx_msg void OnUpdateButtonMove(CCmdUI* pCmdUI);

  afx_msg void OnButtonRotateX();
  afx_msg void OnUpdateButtonRotateX(CCmdUI* pCmdUI);

  afx_msg void OnButtonRotateY();
  afx_msg void OnUpdateButtonRotateY(CCmdUI* pCmdUI);

  afx_msg void OnButtonRotate();
  afx_msg void OnUpdateButtonRotate(CCmdUI* pCmdUI);
	DECLARE_MESSAGE_MAP()

	BOOL CreateDockingWindows();
	void SetDockingWindowIcons(BOOL bHiColorIcons);
};

inline CPropertiesWnd* CMainFrame::GetWndProperties()
{
  return &m_wndProperties;
}

inline CMFCStatusBar* CMainFrame::GetStatusBar()
{
  return &m_wndStatusBar;
}
