// AddFieldPtbDialog.cpp : implementation file
//

#include "stdafx.h"
#include "AddFieldPtbDialog.h"
#include "ParticleTracking.h"

// CAddFieldPtbDialog dialog

IMPLEMENT_DYNAMIC(CAddFieldPtbDialog, CDialog)

CAddFieldPtbDialog::CAddFieldPtbDialog(CWnd* pParent /*=NULL*/)
	: CDialog(CAddFieldPtbDialog::IDD, pParent)
{
}

CAddFieldPtbDialog::~CAddFieldPtbDialog()
{
}

BOOL CAddFieldPtbDialog::OnInitDialog()
{
  CDialog::OnInitDialog();

  m_PtbSelector.Clear();
// I use the InsertString function instead of AddString because this does not result in sorting strings even though the 
// combo-box has been created with the CBS_SORT style.
  for(int i = 0; i < EvaporatingParticle::CFieldPerturbation::ptbCount; i++)
    m_PtbSelector.InsertString(m_PtbSelector.GetCount(), EvaporatingParticle::CFieldPerturbation::perturbation_name(i));

  m_PtbSelector.SetCurSel(EvaporatingParticle::CFieldPerturbation::ptbRing);
  return TRUE;
}

void CAddFieldPtbDialog::DoDataExchange(CDataExchange* pDX)
{
  CDialog::DoDataExchange(pDX);
  DDX_Control(pDX, IDC_SELECT_PTB_COMBO, m_PtbSelector);
}


BEGIN_MESSAGE_MAP(CAddFieldPtbDialog, CDialog)
  ON_BN_CLICKED(IDOK, &CAddFieldPtbDialog::OnBnClickedOk)
END_MESSAGE_MAP()


// CAddFieldPtbDialog message handlers

void CAddFieldPtbDialog::OnBnClickedOk()
{
  if(m_PtbSelector.GetCount() != 0)
  {
    int nSel = m_PtbSelector.GetCurSel();
    EvaporatingParticle::CFieldPerturbation* pFieldPtb = EvaporatingParticle::CFieldPerturbation::create(nSel);
    if(pFieldPtb != NULL)
    {
      EvaporatingParticle::CFieldPtbCollection& coll = CParticleTrackingApp::Get()->GetTracker()->get_field_ptb();
      coll.push_back(pFieldPtb);
    }
  }

  OnOK();
}
