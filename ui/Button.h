
#pragma once

class CPropertiesWnd;
//---------------------------------------------------------------------------------------
// CSelectRegionButton.
//---------------------------------------------------------------------------------------
class CSelectRegionButton : public CMFCPropertyGridProperty
{
  DECLARE_DYNAMIC(CSelectRegionButton)

  BOOL          m_bNew,
                m_bPressed;

public:
  CSelectRegionButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CMFCPropertyGridProperty(strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
    m_pWndProp = pWndProp;

    m_bNew = TRUE;
    m_bAllowEdit = FALSE;
    m_bPressed = FALSE;
  }

  virtual void    OnClickButton(CPoint point);
  virtual void    OnDrawValue(CDC* pDC, CRect rect);
  virtual void    OnDrawButton(CDC* pDC, CRect rectButton);
  virtual void    AdjustButtonRect();
  virtual BOOL    OnUpdateValue();

  virtual void    OnSetSelection(CMFCPropertyGridProperty* /*pOldSel*/);
	virtual void    OnKillSelection(CMFCPropertyGridProperty* /*pNewSel*/);

  virtual CString ButtonValue();

protected:
  virtual BOOL    HasButton() const;

  void            show_all_regions();
  void            process_click();

  CPropertiesWnd* m_pWndProp;
};

//---------------------------------------------------------------------------------------
// CSelectTrajectButton.
//---------------------------------------------------------------------------------------
class CSelectTrajectButton : public CSelectRegionButton
{
  DECLARE_DYNAMIC(CSelectTrajectButton)

public:
  CSelectTrajectButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CSelectRegionButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
  virtual void    OnSetSelection(CMFCPropertyGridProperty* /*pOldSel*/);
	virtual void    OnKillSelection(CMFCPropertyGridProperty* /*pNewSel*/);

  virtual CString ButtonValue();
};

//---------------------------------------------------------------------------------------
// CSelectFacesButton.
//---------------------------------------------------------------------------------------
class CSelectFacesButton : public CSelectRegionButton
{
  DECLARE_DYNAMIC(CSelectFacesButton)

public:
  CSelectFacesButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CSelectRegionButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
  virtual void    OnSetSelection(CMFCPropertyGridProperty* /*pOldSel*/);
	virtual void    OnKillSelection(CMFCPropertyGridProperty* /*pNewSel*/);

  virtual CString ButtonValue();
};

//---------------------------------------------------------------------------------------
// CProprtyListButton - a base class for simple buttons attached to the properties list.
//---------------------------------------------------------------------------------------
class CProprtyListButton : public CMFCPropertyGridProperty
{
  DECLARE_DYNAMIC(CProprtyListButton)

public:
  CProprtyListButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CMFCPropertyGridProperty(strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
    m_pWndProp = pWndProp;
    m_bAllowEdit = FALSE;
  }

  virtual void    OnClickButton(CPoint point) {}

  virtual void    OnDrawValue(CDC* pDC, CRect rect);
  virtual void    OnDrawButton(CDC* pDC, CRect rectButton);
  virtual void    AdjustButtonRect();

protected:
  virtual BOOL    HasButton() const;

  CPropertiesWnd* m_pWndProp;
};

//---------------------------------------------------------------------------------------
// CSelectFolderButton.
//---------------------------------------------------------------------------------------
class CSelectFolderButton : public CProprtyListButton
{
  DECLARE_DYNAMIC(CSelectFolderButton)

public:
  CSelectFolderButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CProprtyListButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CRemovePropertyButton.
//---------------------------------------------------------------------------------------
class CRemovePropertyButton : public CProprtyListButton
{
  DECLARE_DYNAMIC(CRemovePropertyButton)

public:
  CRemovePropertyButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CProprtyListButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnDrawButton(CDC* pDC, CRect rectButton);
  virtual bool    ConfirmRemove() const;
};

//---------------------------------------------------------------------------------------
// CRemoveBoundCondButton.
//---------------------------------------------------------------------------------------
class CRemoveBoundCondButton : public CRemovePropertyButton
{
  DECLARE_DYNAMIC(CRemoveBoundCondButton)

public:
  CRemoveBoundCondButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CRemovePropertyButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
  virtual bool    ConfirmRemove() const;
};

//---------------------------------------------------------------------------------------
// CRemoveContourButton.
//---------------------------------------------------------------------------------------
class CRemoveContourButton : public CRemovePropertyButton
{
  DECLARE_DYNAMIC(CRemoveContourButton)

public:
  CRemoveContourButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CRemovePropertyButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
  virtual bool    ConfirmRemove() const;
};

//---------------------------------------------------------------------------------------
// CClearLocationsButton.
//---------------------------------------------------------------------------------------
class CClearLocationsButton : public CRemovePropertyButton
{
  DECLARE_DYNAMIC(CClearLocationsButton)

public:
  CClearLocationsButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CRemovePropertyButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CClearSelectedFacesButton.
//---------------------------------------------------------------------------------------
class CClearSelectedFacesButton : public CRemovePropertyButton
{
  DECLARE_DYNAMIC(CClearSelectedFacesButton)

public:
  CClearSelectedFacesButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CRemovePropertyButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
  virtual bool    ConfirmRemove() const;
};

//---------------------------------------------------------------------------------------
// CRemoveCalcButton.
//---------------------------------------------------------------------------------------
class CRemoveCalcButton : public CRemovePropertyButton
{
  DECLARE_DYNAMIC(CRemoveCalcButton)

public:
  CRemoveCalcButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CRemovePropertyButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
  virtual bool    ConfirmRemove() const;
};

//---------------------------------------------------------------------------------------
// CStartCalcButton.
//---------------------------------------------------------------------------------------
class CStartCalcButton : public CProprtyListButton
{
  DECLARE_DYNAMIC(CStartCalcButton)

public:
  CStartCalcButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CProprtyListButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
  virtual void    OnDrawButton(CDC* pDC, CRect rectButton);
};

//---------------------------------------------------------------------------------------
// CRemovePerturbationButton.
//---------------------------------------------------------------------------------------
class CRemovePerturbationButton : public CRemovePropertyButton
{
  DECLARE_DYNAMIC(CRemovePerturbationButton)

public:
  CRemovePerturbationButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CRemovePropertyButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CRemoveFieldButton.
//---------------------------------------------------------------------------------------
class CRemoveFieldButton : public CRemovePropertyButton
{
  DECLARE_DYNAMIC(CRemoveFieldButton)

public:
  CRemoveFieldButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CRemovePropertyButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CAddFieldBoundCondButton - a class for adding boundary condition for electric potential.
//---------------------------------------------------------------------------------------
class CAddFieldBoundCondButton : public CProprtyListButton
{
  DECLARE_DYNAMIC(CAddFieldBoundCondButton)

public:
  CAddFieldBoundCondButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CProprtyListButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
  virtual void    OnDrawButton(CDC* pDC, CRect rectButton);
};

//---------------------------------------------------------------------------------------
// CRemoveFieldBoundCondButton.
//---------------------------------------------------------------------------------------
class CRemoveFieldBoundCondButton : public CRemovePropertyButton
{
  DECLARE_DYNAMIC(CRemoveFieldBoundCondButton)

public:
  CRemoveFieldBoundCondButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CRemovePropertyButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CRemoveCrossSectionButton.
//---------------------------------------------------------------------------------------
class CRemoveCrossSectionButton : public CRemovePropertyButton
{
  DECLARE_DYNAMIC(CRemoveCrossSectionButton)

public:
  CRemoveCrossSectionButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CRemovePropertyButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CCalcFieldButton.
//---------------------------------------------------------------------------------------
class CCalcFieldButton : public CProprtyListButton
{
  DECLARE_DYNAMIC(CCalcFieldButton)

public:
  CCalcFieldButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CProprtyListButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
  virtual void    OnDrawButton(CDC* pDC, CRect rectButton);
};

//---------------------------------------------------------------------------------------
// CCalcFieldPtbButton.
//---------------------------------------------------------------------------------------
class CCalcFieldPtbButton : public CCalcFieldButton
{
  DECLARE_DYNAMIC(CCalcFieldPtbButton)

public:
  CCalcFieldPtbButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CCalcFieldButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CSaveFieldButton.
//---------------------------------------------------------------------------------------
class CSaveFieldButton : public CProprtyListButton
{
  DECLARE_DYNAMIC(CSaveFieldButton)

public:
  CSaveFieldButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CProprtyListButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
  }

  virtual void    OnClickButton(CPoint point);
  virtual void    OnDrawButton(CDC* pDC, CRect rectButton);

protected:
  CString         get_file_name(bool bLoadFromFile = true);
};

//---------------------------------------------------------------------------------------
// CSaveLoadSelFacesButton.
//---------------------------------------------------------------------------------------
class CSaveLoadSelFacesButton : public CSaveFieldButton
{
  DECLARE_DYNAMIC(CSaveLoadSelFacesButton)

public:
  CSaveLoadSelFacesButton(bool bSave, CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CSaveFieldButton(pWndProp, strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
    m_bSave = bSave;
  }

  virtual void    OnClickButton(CPoint point);

protected:
  bool            m_bSave;
};

//---------------------------------------------------------------------------------------
// CCheckBoxButton - a base class for processing boolean variables. Use this class if no
// immediate response is expected. However, if some specific actions are required just 
// after changing the flag (redrawing, etc), you have to use the descendants with the 
// OnClickButton function overridden.
//---------------------------------------------------------------------------------------
class CCheckBoxButton : public CProprtyListButton
{
  DECLARE_DYNAMIC(CCheckBoxButton)

public:
  CCheckBoxButton(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr, DWORD_PTR dwData)
    : CProprtyListButton(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual void    OnClickButton(CPoint point);
  virtual void    OnDrawButton(CDC* pDC, CRect rectButton);
};

//---------------------------------------------------------------------------------------
// CRedrawCheckBox - use this class if only an immediate redrwing is needed after the switch.
//---------------------------------------------------------------------------------------
class CRedrawCheckBox : public CCheckBoxButton
{
  DECLARE_DYNAMIC(CRedrawCheckBox)

public:
  CRedrawCheckBox(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr, DWORD_PTR dwData)
    : CCheckBoxButton(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CContourRangeCheckBox
//---------------------------------------------------------------------------------------
class CContourRangeCheckBox : public CCheckBoxButton
{
  DECLARE_DYNAMIC(CContourRangeCheckBox)

public:
  CContourRangeCheckBox(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr, DWORD_PTR dwData)
    : CCheckBoxButton(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CTrackRangeCheckBox
//---------------------------------------------------------------------------------------
class CTrackRangeCheckBox : public CCheckBoxButton
{
  DECLARE_DYNAMIC(CTrackRangeCheckBox)

public:
  CTrackRangeCheckBox(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr, DWORD_PTR dwData)
    : CCheckBoxButton(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CCrossSectCheckBox
//---------------------------------------------------------------------------------------
class CCrossSectCheckBox : public CCheckBoxButton
{
  DECLARE_DYNAMIC(CCrossSectCheckBox)

public:
  CCrossSectCheckBox(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr, DWORD_PTR dwData)
    : CCheckBoxButton(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CSourceCheckBox
//---------------------------------------------------------------------------------------
class CSourceCheckBox : public CCheckBoxButton
{
  DECLARE_DYNAMIC(CSourceCheckBox)

public:
  CSourceCheckBox(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr, DWORD_PTR dwData)
    : CCheckBoxButton(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CInvalidateFieldCheckBox
//---------------------------------------------------------------------------------------
class CInvalidateFieldCheckBox : public CCheckBoxButton
{
  DECLARE_DYNAMIC(CInvalidateFieldCheckBox)

public:
  CInvalidateFieldCheckBox(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr, DWORD_PTR dwData)
    : CCheckBoxButton(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CPlaneYZCalcCheckBox
//---------------------------------------------------------------------------------------
class CPlaneYZCalcCheckBox : public CCheckBoxButton
{
  DECLARE_DYNAMIC(CPlaneYZCalcCheckBox)

public:
  CPlaneYZCalcCheckBox(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr, DWORD_PTR dwData)
    : CCheckBoxButton(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CHideShowRegsCheckBox
//---------------------------------------------------------------------------------------
class CHideShowRegsCheckBox : public CCheckBoxButton
{
  DECLARE_DYNAMIC(CHideShowRegsCheckBox)

public:
  CHideShowRegsCheckBox(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr, DWORD_PTR dwData)
    : CCheckBoxButton(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// CUserDefCSCheckBox - the check-box for immediate cross-section calculation.
//---------------------------------------------------------------------------------------
class CUserDefCSCheckBox : public CCheckBoxButton
{
  DECLARE_DYNAMIC(CUserDefCSCheckBox)

public:
  CUserDefCSCheckBox(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr, DWORD_PTR dwData)
    : CCheckBoxButton(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual void    OnClickButton(CPoint point);
};

//---------------------------------------------------------------------------------------
// Inline implementation.
//---------------------------------------------------------------------------------------
inline BOOL CSelectRegionButton::HasButton() const
{
  return TRUE;
}

inline BOOL CProprtyListButton::HasButton() const
{
  return TRUE;
}