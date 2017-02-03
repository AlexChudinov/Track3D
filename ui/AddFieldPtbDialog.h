#pragma once

#include "afxwin.h"
#include "resource.h"

// CAddFieldPtbDialog dialog

class CAddFieldPtbDialog : public CDialog
{
	DECLARE_DYNAMIC(CAddFieldPtbDialog)

public:
	CAddFieldPtbDialog(CWnd* pParent = NULL);   // standard constructor
	virtual ~CAddFieldPtbDialog();

// Dialog Data
	enum { IDD = IDD_ADD_PERTURB_DIALOG };

protected:
	virtual void    DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
  virtual BOOL    OnInitDialog();

	DECLARE_MESSAGE_MAP()

public:
  CComboBox       m_PtbSelector;

  afx_msg void    OnBnClickedOk();
};
