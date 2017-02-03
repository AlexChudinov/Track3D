
#pragma once

#include "afxwin.h"
#include "resource.h"

//-------------------------------------------------------------------------------------------------
// CAddCalculatorDialog dialog
//-------------------------------------------------------------------------------------------------
class CAddCalculatorDialog : public CDialog
{
	DECLARE_DYNAMIC(CAddCalculatorDialog)

public:
	CAddCalculatorDialog(CWnd* pParent = NULL);   // standard constructor
	virtual ~CAddCalculatorDialog();

// Dialog Data
	enum { IDD = IDD_ADD_CALCULATOR_DIALOG };

protected:
	virtual void    DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
  virtual BOOL    OnInitDialog();


	DECLARE_MESSAGE_MAP()

public:
  afx_msg void    OnBnClickedOk();

  CComboBox       m_CalcSelector;
  
};
