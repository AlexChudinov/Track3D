
#pragma once

class CPropertiesWnd;
//---------------------------------------------------------------------------------------
//  The "Response" classes have the OnUpdateValue() method overridden which allows to use
//  this kind of controls to change the rendering parameters (color, opacity) on the fly.
//---------------------------------------------------------------------------------------
class CResponseProperty : public CMFCPropertyGridProperty
{
public:
  CResponseProperty(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CMFCPropertyGridProperty(strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
    m_pWndProp = pWndProp;
  }

  virtual BOOL    OnUpdateValue();

protected:
  void            on_update_contour_ctrls();

private:
  CPropertiesWnd* m_pWndProp;
};

//---------------------------------------------------------------------------------------
//  CCrossSectionResponer is a response type of property for cross-sections.
//---------------------------------------------------------------------------------------
class CCrossSectionResponer : public CMFCPropertyGridProperty
{
public:
  CCrossSectionResponer(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CMFCPropertyGridProperty(strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
    m_pWndProp = pWndProp;
  }

  virtual BOOL    OnUpdateValue();

private:
  CPropertiesWnd* m_pWndProp;
};

//---------------------------------------------------------------------------------------
//  CCalcResponseProperty is a response type of property specially used for calculators.
//---------------------------------------------------------------------------------------
class CCalcResponseProperty : public CMFCPropertyGridProperty
{
public:
  CCalcResponseProperty(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CMFCPropertyGridProperty(strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
    m_pWndProp = pWndProp;
  }

  virtual BOOL    OnUpdateValue();

private:
  CPropertiesWnd* m_pWndProp;
};

//---------------------------------------------------------------------------------------
//  CColorResponseProperty
//---------------------------------------------------------------------------------------
class CColorResponseProperty : public CMFCPropertyGridColorProperty
{
public:
  CColorResponseProperty(const CString& strName, const COLORREF& color, CPalette* pPalette = NULL, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0)
    : CMFCPropertyGridColorProperty(strName, color, pPalette, lpszDescr, dwData)
  {
  }

  virtual BOOL  OnUpdateValue();
};
