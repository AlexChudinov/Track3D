
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

//---------------------------------------------------------------------------------------
//  CAreaFacesColorResponder
//---------------------------------------------------------------------------------------
class CAreaFacesColorResponder : public CMFCPropertyGridColorProperty
{
public:
  CAreaFacesColorResponder(const CString& strName, const COLORREF& color, CPalette* pPalette = NULL, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0)
    : CMFCPropertyGridColorProperty(strName, color, pPalette, lpszDescr, dwData)
  {
  }

  virtual BOOL  OnUpdateValue();
};

//---------------------------------------------------------------------------------------
//  CGeneralResponseProperty
//---------------------------------------------------------------------------------------
class CGeneralResponseProperty : public CMFCPropertyGridProperty
{
public:
  CGeneralResponseProperty(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0,
		LPCTSTR lpszEditMask = NULL, LPCTSTR lpszEditTemplate = NULL, LPCTSTR lpszValidChars = NULL)
    : CMFCPropertyGridProperty(strName, varValue, lpszDescr, dwData, lpszEditMask, lpszEditTemplate, lpszValidChars)
  {
    m_pWndProp = pWndProp;
  }

  virtual BOOL    OnUpdateValue();

protected:
  CPropertiesWnd* m_pWndProp;
};

//---------------------------------------------------------------------------------------
//  CElectricFieldResponder
//---------------------------------------------------------------------------------------
class CElectricFieldResponder : public CGeneralResponseProperty
{
public:
  CElectricFieldResponder(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0)
    : CGeneralResponseProperty(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual BOOL    OnUpdateValue();
};

//---------------------------------------------------------------------------------------
//  CNamedAreasSelResponder
//---------------------------------------------------------------------------------------
class CNamedAreasSelResponder : public CGeneralResponseProperty
{
public:
  CNamedAreasSelResponder(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0)
    : CGeneralResponseProperty(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

  virtual BOOL    OnUpdateValue();
};

//---------------------------------------------------------------------------------------
//  CSetAndRedrawResponder. 
//---------------------------------------------------------------------------------------
class CSetAndRedrawResponder : public CGeneralResponseProperty
{
public:
  CSetAndRedrawResponder(CPropertiesWnd* pWndProp, const CString& strName, const COleVariant& varValue, LPCTSTR lpszDescr = NULL, DWORD_PTR dwData = 0)
    : CGeneralResponseProperty(pWndProp, strName, varValue, lpszDescr, dwData)
  {
  }

// Note: this function does NOT set the value from the control to the model. It just calls m_pWndProp->set_data_to_model().
// Thus: make sure a proper set-function is written and called from CPropertiesWnd::set_some_data().
  virtual BOOL    OnUpdateValue();
};

