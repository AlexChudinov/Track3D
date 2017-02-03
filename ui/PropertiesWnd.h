
#pragma once

class CPropertiesToolBar : public CMFCToolBar
{
public:
	virtual void OnUpdateCmdUI(CFrameWnd* /*pTarget*/, BOOL bDisableIfNoHndler)
	{
		CMFCToolBar::OnUpdateCmdUI((CFrameWnd*) GetOwner(), bDisableIfNoHndler);
	}

	virtual BOOL AllowShowOnList() const { return FALSE; }
};

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
class CPropertiesWnd : public CDockablePane
{
// Construction
public:
	CPropertiesWnd();

	void AdjustLayout();

  enum
  {
    tabMain     = 0,
    tabSource   = 1,
    tabTrack    = 2,
    tabEvapor   = 3,
    tabIonPar   = 4,
    tabDrawing  = 5,
    tabExport   = 6,
    tabImport   = 7,
    tabCalc     = 8,
    tabPtb      = 9,
    tabFields   = 10,
    nTabCount   = 11
  };

// Attributes
public:
	void SetVSDotNetLook(BOOL bSet)
	{
		m_wndPropList.SetVSDotNetLook(bSet);
		m_wndPropList.SetGroupNameFullWidth(bSet);
	}

protected:
	CFont                 m_fntPropList;
	CComboBox             m_wndObjectCombo;
	CPropertiesToolBar    m_wndToolBar;
	CMFCPropertyGridCtrl  m_wndPropList;
  CTabCtrl              m_wndTabCtrl;

  bool                  m_bBusy;
  bool                  m_bUpdateAll;

// Implementation
public:
	virtual ~CPropertiesWnd();

  void          on_idle();
  void          set_update_ctrls();
  void          set_update_all();

  void          set_data_to_model();
  void          update_ctrls();

protected:
  void          update_prop_list();

  void          add_type_ctrls();
  void          add_source_ctrls();
  void          add_tracking_ctrls();
  void          add_evapor_ctrls();
  void          add_ion_ctrls();
  void          add_draw_ctrls();
  void          add_cs_plane_ctrls(CMFCPropertyGridProperty* pDrawGroup);
  void          add_contour_ctrls(CMFCPropertyGridProperty* pDrawGroup);
  void          add_export_ctrls();
  void          add_bc_ctrls(CMFCPropertyGridProperty* pOpenFOAMGroup);
  void          add_import_ctrls();
  void          add_calc_ctrls();
  void          add_ptb_ctrls();
  void          add_field_ctrls();

  void          set_type_data();
  void          set_source_data();
  void          set_tracking_data();
  void          set_evapor_data();
  void          set_ion_data();
  void          set_draw_data();
  void          set_export_data();
  void          set_bc_data();
  void          set_import_data();
  void          set_calc_data();
  void          set_ptb_data();
  void          set_field_data();

  void          update_type_ctrls();
  void          update_source_ctrls();
  void          update_tracking_ctrls();
  void          update_evapor_ctrls();
  void          update_ion_ctrls();
  void          update_draw_ctrls();
  void          update_cs_plane_ctrls();
  void          update_contour_ctrls();
  void          update_export_ctrls();
  void          update_bc_ctrls(bool bEnable);
  void          update_import_ctrls();
  void          update_calc_ctrls();
  void          update_ptb_ctrls();
  void          update_field_ctrls();

	afx_msg int   OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void  OnSize(UINT nType, int cx, int cy);
	afx_msg void  OnExpandAllProperties();
	afx_msg void  OnUpdateExpandAllProperties(CCmdUI* pCmdUI);
	afx_msg void  OnSortProperties();
	afx_msg void  OnUpdateSortProperties(CCmdUI* pCmdUI);
	afx_msg void  OnDoTracking();
	afx_msg void  OnUpdateDoTracking(CCmdUI* pCmdUI);
  afx_msg void  OnAddFieldPeturbation();
	afx_msg void  OnUpdateAddFieldPeturbation(CCmdUI* pCmdUI);
	afx_msg void  OnSaveImage();
	afx_msg void  OnUpdateSaveImage(CCmdUI* pCmdUI);
  afx_msg void  OnExportMesh();
	afx_msg void  OnUpdateExportMesh(CCmdUI* pCmdUI);
  afx_msg void  OnCreateBoundaryConditions();
	afx_msg void  OnUpdateCreateBoundaryConditions(CCmdUI* pCmdUI);
  afx_msg void  OnCreateContour();
	afx_msg void  OnUpdateCreateContour(CCmdUI* pCmdUI);
  afx_msg void  OnCreateCalculator();
  afx_msg void  OnUpdateCreateCalculator(CCmdUI* pCmdUI);
  afx_msg void  OnImportOpenFOAM();
	afx_msg void  OnUpdateImportOpenFOAM(CCmdUI* pCmdUI);
  afx_msg void  OnAddField();
	afx_msg void  OnUpdateAddField(CCmdUI* pCmdUI);
  afx_msg void  OnAddCrossSection();
	afx_msg void  OnUpdateAddCrossSection(CCmdUI* pCmdUI);
	afx_msg void  OnSetFocus(CWnd* pOldWnd);
	afx_msg void  OnSettingChange(UINT uFlags, LPCTSTR lpszSection);
  afx_msg void  OnTabSelChanging(NMHDR* pNMHDR, LRESULT* pLResult);
  afx_msg void  OnTabSelChange(NMHDR* pNMHDR, LRESULT* pLResult);

	DECLARE_MESSAGE_MAP()

	void          InitPropList();
	void          SetPropListFont();

  virtual LRESULT WindowProc(UINT message, WPARAM wParam, LPARAM lParam);
};

//---------------------------------------------------------------------------------------
// Inline implementation.
//---------------------------------------------------------------------------------------
inline void CPropertiesWnd::set_update_all()
{
  m_bUpdateAll = true;
}
