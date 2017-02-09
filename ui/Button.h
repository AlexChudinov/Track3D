
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

protected:
  virtual BOOL    HasButton() const;

private:
  CPropertiesWnd* m_pWndProp;
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