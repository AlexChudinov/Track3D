// AddCalculatorDialog.cpp : implementation file
//

#include "stdafx.h"

#include "AddCalculatorDialog.h"
#include "ParticleTracking.h"
#include "MainFrm.h"

//-------------------------------------------------------------------------------------------------
// CAddCalculatorDialog dialog
//-------------------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CAddCalculatorDialog, CDialog)

CAddCalculatorDialog::CAddCalculatorDialog(CWnd* pParent /*=NULL*/)
	: CDialog(CAddCalculatorDialog::IDD, pParent)
{
}

CAddCalculatorDialog::~CAddCalculatorDialog()
{
}

BOOL CAddCalculatorDialog::OnInitDialog()
{
  CDialog::OnInitDialog();

  m_CalcSelector.Clear();
  for(int i = 0; i < EvaporatingParticle::CCalculator::ctCount; i++)
    m_CalcSelector.AddString(EvaporatingParticle::CCalculator::calc_name(i));

  m_CalcSelector.SetCurSel(EvaporatingParticle::CCalculator::ctAlongLine);
  return TRUE;
}

void CAddCalculatorDialog::DoDataExchange(CDataExchange* pDX)
{
  CDialog::DoDataExchange(pDX);
  DDX_Control(pDX, IDC_SELECT_CALC_COMBO, m_CalcSelector);
}


BEGIN_MESSAGE_MAP(CAddCalculatorDialog, CDialog)
  ON_BN_CLICKED(IDOK, &CAddCalculatorDialog::OnBnClickedOk)
END_MESSAGE_MAP()


// CAddCalculatorDialog message handlers

void CAddCalculatorDialog::OnBnClickedOk()
{
  if(m_CalcSelector.GetCount() != 0)
  {
    int nSel = m_CalcSelector.GetCurSel();
    EvaporatingParticle::CCalculator* pCalc = EvaporatingParticle::CCalculator::create(nSel);
    if(pCalc != NULL)
    {
      EvaporatingParticle::CCalcCollection* pCalcColl = CParticleTrackingApp::Get()->GetCalcs();
      pCalcColl->push_back(pCalc);
    }
  }

  OnOK();
}
