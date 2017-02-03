
#include "stdafx.h"

#include "Button.h"
#include "ParticleTracking.h"
#include "PropertiesWnd.h"

#include "ColorContour.h"
#include "OutputEngine.h"
#include "ExecutionDialog.h"
#include "Perturbation.h"

//---------------------------------------------------------------------------------------
//  CSelectRegionButton
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CSelectRegionButton, CMFCPropertyGridProperty)

void CSelectRegionButton::OnClickButton(CPoint point)
{
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(m_dwData == NULL)
  {
    pDrawObj->show_all_regions();
    CSelectRegionButton* pBtn = (CSelectRegionButton*)(m_pWndList->FindItemByData(pDrawObj->get_hidden_reg_names_ptr()));
    if(pBtn != NULL)
    {
      pBtn->SetValue(CString(""));
      pBtn->Redraw();
    }
  }
  else
  {
    EvaporatingParticle::CStringVector* pRegNames = (EvaporatingParticle::CStringVector*)m_dwData;
    if(pDrawObj->get_sel_flag())
    {
      pDrawObj->exit_sel_context(pRegNames);  // here pRegNames is updated using the information stored in the Draw Object.
      SetValue(EvaporatingParticle::CObject::compile_string(*pRegNames));
      pDrawObj->invalidate_contour(m_dwData);
      m_bPressed = FALSE;
      Redraw();
    }
    else
    {
      pDrawObj->enter_sel_context(pRegNames);
      m_bPressed = TRUE;
      Redraw();
    }
  }

  pDrawObj->draw();
}

void CSelectRegionButton::OnSetSelection(CMFCPropertyGridProperty* /*pOldSel*/)
{
  if(m_dwData == NULL)
    return;

  EvaporatingParticle::CStringVector* pRegNames = (EvaporatingParticle::CStringVector*)m_dwData;
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->enter_sel_context(pRegNames, false);
  pDrawObj->draw();
}

void CSelectRegionButton::OnKillSelection(CMFCPropertyGridProperty* /*pNewSel*/)
{
  if(m_dwData == NULL)
    return;

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->hide_selected();
  pDrawObj->draw();
}

void CSelectRegionButton::OnDrawValue(CDC* pDC, CRect rect)
{
  AdjustButtonRect();

  rect.left += m_rectButton.Width();
	CMFCPropertyGridProperty::OnDrawValue(pDC, rect);

  OnDrawButton(pDC, m_rectButton);
}

void CSelectRegionButton::OnDrawButton(CDC* pDC, CRect rect)
{
  COLORREF crBkColor = pDC->GetBkColor();
  COLORREF crTextColor = pDC->GetTextColor();
  
  CRect rectButton = rect;
	rectButton.DeflateRect(1, 1);
	rectButton.top++;
	rectButton.left++;

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  COLORREF crBrushColor = m_bPressed ? clAqua : clLtGray;

	CBrush br(crBrushColor);
	pDC->FillRect(rectButton, &br);
	pDC->Draw3dRect(rectButton, 0, 0);

  pDC->SetTextColor(crTextColor);
  pDC->SetBkColor(crBkColor);
}

void CSelectRegionButton::AdjustButtonRect()
{
  m_rectButton = m_Rect;
  m_rectButton.left = m_pWndList->GetListRect().left + m_pWndList->GetLeftColumnWidth() + 2;
  m_rectButton.right = m_rectButton.left + m_rectButton.Height();
  m_rectButton.top ++;
}

BOOL CSelectRegionButton::OnUpdateValue()
{
  if(m_bNew)
  {
    m_pWndList->OnPropertyChanged(this);
    m_bNew = FALSE;
  }

  return TRUE;
}

//---------------------------------------------------------------------------------------
// CProprtyListButton - a base class for simple buttons attached to the properties list.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CProprtyListButton, CMFCPropertyGridProperty)

void CProprtyListButton::OnDrawValue(CDC* pDC, CRect rect)
{
  AdjustButtonRect();

  rect.left += m_rectButton.Width();
	CMFCPropertyGridProperty::OnDrawValue(pDC, rect);

  OnDrawButton(pDC, m_rectButton);
}

void CProprtyListButton::OnDrawButton(CDC* pDC, CRect rect)
{
  COLORREF crBkColor = pDC->GetBkColor();
  COLORREF crTextColor = pDC->GetTextColor();
  
  CRect rectButton = rect;
	rectButton.DeflateRect(1, 1);
	rectButton.top++;
	rectButton.left++;

  COLORREF crBrushColor = clLtGray;

	CBrush br(crBrushColor);
	pDC->FillRect(rectButton, &br);
	pDC->Draw3dRect(rectButton, 0, 0);

  pDC->SetTextColor(crTextColor);
  pDC->SetBkColor(crBkColor);
}

void CProprtyListButton::AdjustButtonRect()
{
  m_rectButton = m_Rect;
  m_rectButton.left = m_pWndList->GetListRect().left + m_pWndList->GetLeftColumnWidth() + 2;
  m_rectButton.right = m_rectButton.left + m_rectButton.Height();
  m_rectButton.top++;
  m_rectButton.right--;
}

//---------------------------------------------------------------------------------------
// CRemovePropertyButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CRemovePropertyButton, CProprtyListButton)

void CRemovePropertyButton::OnDrawButton(CDC* pDC, CRect rectButton)
{
  COLORREF crBkColor = pDC->GetBkColor();
  COLORREF crTextColor = pDC->GetTextColor();

  COLORREF crBrushColor = RGB(180, 0, 0);
  CBrush brush1(crBrushColor);
  pDC->FillRect(rectButton, &brush1);

  int nx0 = rectButton.left + 5, ny0 = rectButton.top + 5;
  int nx1 = rectButton.right - 5, ny1 = rectButton.bottom - 5;
  CBrush brush2(clWhite);
  pDC->MoveTo(nx0, ny0);
  pDC->LineTo(nx1, ny1);

  pDC->MoveTo(nx0 - 1, ny0);
  pDC->LineTo(nx1 - 1, ny1);

  pDC->MoveTo(nx0 + 1, ny0);
  pDC->LineTo(nx1 + 1, ny1);

  nx0 = rectButton.right - 6;
  nx1 = rectButton.left + 4;

  pDC->MoveTo(nx0, ny0);
  pDC->LineTo(nx1, ny1);

  pDC->MoveTo(nx0 - 1, ny0);
  pDC->LineTo(nx1 - 1, ny1);

  pDC->MoveTo(nx0 + 1, ny0);
  pDC->LineTo(nx1 + 1, ny1);

  pDC->SetTextColor(crTextColor);
  pDC->SetBkColor(crBkColor);
}

//---------------------------------------------------------------------------------------
// CRemoveBoundCondButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CRemoveBoundCondButton, CRemovePropertyButton)

void CRemoveBoundCondButton::OnClickButton(CPoint point)
{
  EvaporatingParticle::CBoundaryConditions* pBC = (EvaporatingParticle::CBoundaryConditions*)m_dwData;
  EvaporatingParticle::CExportOpenFOAM* pExpObj = CParticleTrackingApp::Get()->GetExporter();
  pExpObj->remove_bound_cond(pBC);
  m_pWndProp->set_update_all();
}

//---------------------------------------------------------------------------------------
//  CRemoveContourButton
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CRemoveContourButton, CRemovePropertyButton)

void CRemoveContourButton::OnClickButton(CPoint point)
{
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  EvaporatingParticle::CColorContour* pCont = (EvaporatingParticle::CColorContour*)m_dwData;
  pDrawObj->remove_contour(pCont);

  m_pWndProp->set_update_all();
  pDrawObj->draw();
}

//---------------------------------------------------------------------------------------
// CRemoveCalcButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CRemoveCalcButton, CRemovePropertyButton)

void CRemoveCalcButton::OnClickButton(CPoint point)
{
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  EvaporatingParticle::CCalcCollection* pCalcColl = CParticleTrackingApp::Get()->GetCalcs();
  EvaporatingParticle::CCalculator* pCalc = (EvaporatingParticle::CCalculator*)m_dwData;
  size_t nCalcCount = pCalcColl->size();
  for(size_t i = 0; i < nCalcCount; i++)
  {
    if(pCalcColl->at(i) == pCalc)
    {
      pCalcColl->erase(pCalcColl->begin() + i);
      delete pCalc;
      break;
    }
  }

  m_pWndProp->set_update_all();
  pDrawObj->draw();
}

//---------------------------------------------------------------------------------------
// CStartCalcButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CStartCalcButton, CProprtyListButton)

static UINT __stdcall run_calc_func(LPVOID pData)
{
  CExecutionDialog* pDlg = (CExecutionDialog*)pData;
  EvaporatingParticle::CCalculator* pCalc = (EvaporatingParticle::CCalculator*)(pDlg->GetDialogObject());
  pCalc->run();

  if(!pCalc->get_terminate_flag()) // if the user has not terminated the dialog manually, do it after calculations are over.
    pDlg->PostMessage(WM_CLOSE);

  return 0;
}

void CStartCalcButton::OnClickButton(CPoint point)
{
  if(m_pWndProp == NULL)
    return;

  m_pWndProp->set_data_to_model();

  EvaporatingParticle::CCalculator* pCalc = (EvaporatingParticle::CCalculator*)m_dwData;
  if(pCalc == NULL)
    return;

  CExecutionDialog dlg(&run_calc_func, (EvaporatingParticle::CObject*)pCalc);
  INT_PTR nRes = dlg.DoModal();

  m_pWndProp->set_update_all();
}

void CStartCalcButton::OnDrawButton(CDC* pDC, CRect rectButton)
{
  COLORREF crBkColor = pDC->GetBkColor();
  COLORREF crTextColor = pDC->GetTextColor();

  COLORREF crBrushColor = RGB(0, 120, 0);
  CBrush brush1(crBrushColor);
  pDC->FillRect(rectButton, &brush1);

  int nx0 = rectButton.left + 2;
  int ny0 = rectButton.top + (rectButton.bottom - rectButton.top) / 2;
  int nx1 = rectButton.right - 3;
  CBrush brush2(clWhite);
  pDC->MoveTo(nx0, ny0);
  pDC->LineTo(nx1, ny0);

  pDC->MoveTo(nx0, ny0 - 1);
  pDC->LineTo(nx1 - 1, ny0 - 1);

  pDC->MoveTo(nx0, ny0 + 1);
  pDC->LineTo(nx1 - 1, ny0 + 1);

  int nL = (nx1 - nx0) / 2;   // arrow length.
  int nW = 1 + nL / 2;        // arrow half width.
  int nx2 = nx1 - nL + 1;
  int ny1 = ny0 + nW;
  int ny2 = ny0 - nW - 1;

  while(nx2 < nx1)
  {
    pDC->MoveTo(nx2, ny1);
    pDC->LineTo(nx2, ny2);
    nx2++;
    ny2++;
    ny1--;
  }

  pDC->SetTextColor(crTextColor);
  pDC->SetBkColor(crBkColor);
}

//---------------------------------------------------------------------------------------
// CRemovePerturbationButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CRemovePerturbationButton, CRemovePropertyButton)

void CRemovePerturbationButton::OnClickButton(CPoint point)
{
  EvaporatingParticle::CFieldPtbCollection& coll = CParticleTrackingApp::Get()->GetTracker()->get_field_ptb();
  EvaporatingParticle::CFieldPerturbation* pPtb = (EvaporatingParticle::CFieldPerturbation*)m_dwData;
  size_t nPtbCount = coll.size();
  for(size_t i = 0; i < nPtbCount; i++)
  {
    if(coll.at(i) == pPtb)
    {
      coll.erase(coll.begin() + i);
      delete pPtb;
      break;
    }
  }

  m_pWndProp->set_update_all();
}

//---------------------------------------------------------------------------------------
// CRemoveFieldButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CRemoveFieldButton, CRemovePropertyButton)

void CRemoveFieldButton::OnClickButton(CPoint point)
{
  EvaporatingParticle::CFieldDataColl* pColl = CParticleTrackingApp::Get()->GetFields();
  EvaporatingParticle::CElectricFieldData* pData = (EvaporatingParticle::CElectricFieldData*)m_dwData;
  size_t nCount = pColl->size();
  for(size_t i = 0; i < nCount; i++)
  {
    if(pColl->at(i) == pData)
    {
      pColl->erase(pColl->begin() + i);
      delete pData;
      break;
    }
  }

  m_pWndProp->set_update_all();
}

//---------------------------------------------------------------------------------------
// CAddFieldBoundCondButton - a class for adding boundary condition for electric potential.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CAddFieldBoundCondButton, CProprtyListButton)

void CAddFieldBoundCondButton::OnClickButton(CPoint point)
{
  EvaporatingParticle::CElectricFieldData* pData = (EvaporatingParticle::CElectricFieldData*)m_dwData;
  pData->add_bc();

  m_pWndProp->set_update_all();
}

void CAddFieldBoundCondButton::OnDrawButton(CDC* pDC, CRect rectButton)
{
  COLORREF crBkColor = pDC->GetBkColor();
  COLORREF crTextColor = pDC->GetTextColor();

  COLORREF crBrushColor = RGB(0, 120, 0);
  CBrush brush1(crBrushColor);
  pDC->FillRect(rectButton, &brush1);

  int nx0 = rectButton.left + 3;
  int ny0 = rectButton.top + (rectButton.bottom - rectButton.top) / 2;
  int nx1 = rectButton.right - 3;
  CBrush brush2(clWhite);
  pDC->MoveTo(nx0, ny0);
  pDC->LineTo(nx1, ny0);

  pDC->MoveTo(nx0, ny0 - 1);
  pDC->LineTo(nx1, ny0 - 1);

  pDC->MoveTo(nx0, ny0 + 1);
  pDC->LineTo(nx1, ny0 + 1);

  nx0 = rectButton.left + (rectButton.right - rectButton.left) / 2;
  ny0 = rectButton.top + 4;
  int ny1 = rectButton.bottom - 3;

  pDC->MoveTo(nx0, ny0);
  pDC->LineTo(nx0, ny1);

  pDC->MoveTo(nx0 - 1, ny0);
  pDC->LineTo(nx0 - 1, ny1);

  pDC->MoveTo(nx0 + 1, ny0);
  pDC->LineTo(nx0 + 1, ny1);

  pDC->SetTextColor(crTextColor);
  pDC->SetBkColor(crBkColor);
}

//---------------------------------------------------------------------------------------
// CRemoveFieldBoundCondButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CRemoveFieldBoundCondButton, CRemovePropertyButton)

void CRemoveFieldBoundCondButton::OnClickButton(CPoint point)
{
  EvaporatingParticle::CPotentialBoundCond* pBC = (EvaporatingParticle::CPotentialBoundCond*)m_dwData;
  EvaporatingParticle::CFieldDataColl* pColl = CParticleTrackingApp::Get()->GetFields();
  if(pColl->remove_bound_cond(pBC))
    m_pWndProp->set_update_all();
}

//---------------------------------------------------------------------------------------
// CRemoveCrossSectionButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CRemoveCrossSectionButton, CRemovePropertyButton)

void CRemoveCrossSectionButton::OnClickButton(CPoint point)
{
  EvaporatingParticle::CDomainCrossSection* pPlane = (EvaporatingParticle::CDomainCrossSection*)m_dwData;
  EvaporatingParticle::CCrossSectColl* pColl = CParticleTrackingApp::Get()->GetPlanes();
  size_t nPlanesCount = pColl->size();
  for(size_t i = 0; i < nPlanesCount; i++)
  {
    if(pPlane == pColl->at(i))
    {
      pColl->erase(pColl->begin() + i);
// TO DO: notify the contours, which could be using this plane.
      delete pPlane;
      break;
    }
  }

  m_pWndProp->set_update_all();
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->invalidate_hidden();
  pDrawObj->draw();
}

//---------------------------------------------------------------------------------------
// CCalcFieldButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CCalcFieldButton, CProprtyListButton)

static UINT __stdcall calc_field_func(LPVOID pData)
{
  CExecutionDialog* pDlg = (CExecutionDialog*)pData;
  EvaporatingParticle::CElectricFieldData* pField = (EvaporatingParticle::CElectricFieldData*)(pDlg->GetDialogObject());
  pField->calc_field();

  if(!pField->get_terminate_flag()) // if the user has not terminated the dialog manually, do it after calculations are over.
    pDlg->PostMessage(WM_CLOSE);

  return 0;
}

void CCalcFieldButton::OnClickButton(CPoint point)
{
  if(m_pWndProp == NULL)
    return;

  m_pWndProp->set_data_to_model();

  EvaporatingParticle::CElectricFieldData* pData = (EvaporatingParticle::CElectricFieldData*)m_dwData;
  if(pData == NULL)
    return;

  CExecutionDialog dlg(&calc_field_func, (EvaporatingParticle::CObject*)pData);
  INT_PTR nRes = dlg.DoModal();

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->invalidate_contours();
  pDrawObj->draw();
}

void CCalcFieldButton::OnDrawButton(CDC* pDC, CRect rectButton)
{
  COLORREF crBkColor = pDC->GetBkColor();
  COLORREF crTextColor = pDC->GetTextColor();

  COLORREF crBrushColor = RGB(0, 120, 0);
  CBrush brush1(crBrushColor);
  pDC->FillRect(rectButton, &brush1);

  int nx0 = rectButton.left + 2;
  int ny0 = rectButton.top + (rectButton.bottom - rectButton.top) / 2;
  int nx1 = rectButton.right - 3;
  CBrush brush2(clWhite);
  pDC->MoveTo(nx0, ny0);
  pDC->LineTo(nx1, ny0);

  pDC->MoveTo(nx0, ny0 - 1);
  pDC->LineTo(nx1 - 1, ny0 - 1);

  pDC->MoveTo(nx0, ny0 + 1);
  pDC->LineTo(nx1 - 1, ny0 + 1);

  int nL = (nx1 - nx0) / 2;   // arrow length.
  int nW = 1 + nL / 2;        // arrow half width.
  int nx2 = nx1 - nL + 1;
  int ny1 = ny0 + nW;
  int ny2 = ny0 - nW - 1;

  while(nx2 < nx1)
  {
    pDC->MoveTo(nx2, ny1);
    pDC->LineTo(nx2, ny2);
    nx2++;
    ny2++;
    ny1--;
  }

  pDC->SetTextColor(crTextColor);
  pDC->SetBkColor(crBkColor);
}
