
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
    m_pWndProp->set_data_to_model();
    EvaporatingParticle::CStringVector* pRegNames = (EvaporatingParticle::CStringVector*)m_dwData;
    if(pDrawObj->get_sel_flag())
    {
      pDrawObj->exit_sel_context(pRegNames);  // here pRegNames is updated using the information stored in the Draw Object.
      SetValue(EvaporatingParticle::CObject::compile_string(*pRegNames));
      pDrawObj->invalidate_contour(m_dwData);
      m_bPressed = FALSE;

      m_pWndProp->set_ignore_idle(false);
      m_pWndProp->set_update_all();
      m_pWndProp->enable_tab_ctrl();
    }
    else
    {
      m_pWndProp->set_ignore_idle(true);
      m_pWndProp->disable_all_but_one(this);
      m_pWndProp->enable_tab_ctrl(FALSE);

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

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pDrawObj->get_sel_flag())
    return;

  EvaporatingParticle::CStringVector* pRegNames = (EvaporatingParticle::CStringVector*)m_dwData;
  pDrawObj->enter_sel_context(pRegNames, false);
  pDrawObj->draw();
}

void CSelectRegionButton::OnKillSelection(CMFCPropertyGridProperty* /*pNewSel*/)
{
  if(m_dwData == NULL)
    return;

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pDrawObj->get_sel_flag())
    return;

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

  rect.left += m_rectButton.Width() + 2;
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

  rectButton.DeflateRect(1, 1);
  COLORREF crBrushColor = RGB(180, 0, 0);
  CBrush brush1(crBrushColor);
  pDC->FillRect(rectButton, &brush1);

  int nx0 = rectButton.left + 4, ny0 = rectButton.top + 4;
  int nx1 = rectButton.right - 4, ny1 = rectButton.bottom - 4;
  CBrush brush2(clWhite);
  pDC->MoveTo(nx0, ny0);
  pDC->LineTo(nx1, ny1);

  pDC->MoveTo(nx0 - 1, ny0);
  pDC->LineTo(nx1 - 1, ny1);

  pDC->MoveTo(nx0 + 1, ny0);
  pDC->LineTo(nx1 + 1, ny1);

  nx0 = rectButton.right - 5;
  nx1 = rectButton.left + 3;

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
// CClearLocationsButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CClearLocationsButton, CRemovePropertyButton)

void CClearLocationsButton::OnClickButton(CPoint point)
{
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  EvaporatingParticle::CColorContour* pCont = (EvaporatingParticle::CColorContour*)m_dwData;
  pCont->clear_reg_names();
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
  pField->calc_field(true);

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

//---------------------------------------------------------------------------------------
// CSaveFieldButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CSaveFieldButton, CProprtyListButton)

void CSaveFieldButton::OnClickButton(CPoint point)
{
  CString cFilter = "CSV Data File|*.csv||";
	CFileDialog dlg(FALSE, NULL, NULL, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT | OFN_EXPLORER, cFilter);
  dlg.m_ofn.nFilterIndex = 1;

	HRESULT hResult = (int)dlg.DoModal();
	if(FAILED(hResult))
		return;

	CString cFileName = dlg.m_ofn.lpstrFile, cExt;

// Add the file extension if the user didn't supply one
	if(dlg.m_ofn.nFileExtension == 0) 
		cFileName = cFileName + ".csv";

  std::string cName((const char*)cFileName);
  bool bFilled = cName.size() != 0;
  if(bFilled)
  {
    EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
    pObj->save_coulomb_field((const char*)cFileName);
  }
}

void CSaveFieldButton::OnDrawButton(CDC* pDC, CRect rectButton)
{
// TO DO: draw a non-trivial icon here.
  CProprtyListButton::OnDrawButton(pDC, rectButton);
}

//---------------------------------------------------------------------------------------
// CCheckBoxButton.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CCheckBoxButton, CProprtyListButton)

void CCheckBoxButton::OnClickButton(CPoint point)
{
  bool& bChecked = *(bool*)m_dwData;
  bChecked = !bChecked;
  SetValue((_variant_t)bChecked);
  Redraw();
}

void CCheckBoxButton::OnDrawButton(CDC* pDC, CRect rect)
{
  COLORREF crBkColor = pDC->GetBkColor();
  COLORREF crTextColor = pDC->GetTextColor();
  COLORREF crPenColor = pDC->GetDCPenColor();
  
  CRect box = rect;
	box.top++;
	box.left++;
  box.right++;

  CBrush OrigBrush;
	CBrush br1(clNavy);

  box.DeflateRect(1, 1);
  pDC->FrameRect(&box, &br1);

  bool bChecked = *(bool*)m_dwData;
  if(bChecked)
  {
    COLORREF clr = IsEnabled() ? RGB(0, 128, 0) : RGB(128, 128, 128);
    CBrush br2(clr);
    box.DeflateRect(3, 3);
    int nCx = (box.left + box.right) / 2;
    int nCy = (box.top + box.bottom) / 2;

    CRect rc = box;
    rc.top = nCy - 3;
    int nRight = nCx - 1;
    for(int i = box.left; i <= nRight; i++)
    {
      rc.left = i;
      rc.right = rc.left + 1;
      rc.bottom = i < nRight ? rc.top + 4 : rc.top + 3;

      pDC->FillRect(rc, &br2);

      rc.top += 1;
    }

    rc.top -= 2;
    nRight = box.right;
    for(int i = nCx - 1; i < nRight; i++)
    {
      rc.left = i;
      rc.right = rc.left + 1;
      rc.bottom = rc.top + 4;

      pDC->FillRect(rc, &br2);

      rc.top -= 1;
    }
  }

  pDC->SetTextColor(crTextColor);
  pDC->SetBkColor(crBkColor);
  pDC->SetDCPenColor(crPenColor);
}

//---------------------------------------------------------------------------------------
// CRedrawCheckBox - use this class if only an immediate redrwing is needed after the switch.
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CRedrawCheckBox, CCheckBoxButton)

void CRedrawCheckBox::OnClickButton(CPoint point)
{
  CCheckBoxButton::OnClickButton(point);

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->draw();
}

//---------------------------------------------------------------------------------------
// CContourRangeCheckBox
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CContourRangeCheckBox, CCheckBoxButton)

void CContourRangeCheckBox::OnClickButton(CPoint point)
{
  CCheckBoxButton::OnClickButton(point);

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  size_t nContoursCount = pDrawObj->get_contours_count();
  for(size_t i = 0; i < nContoursCount; i++)
  {
    EvaporatingParticle::CColorContour* pObj = pDrawObj->get_contour(i);
    if(pObj->get_enable_user_range_ptr() == m_dwData)
    {
      if(pObj->get_enable_user_range())
        pObj->restore_user_range();
      else
        pObj->get_min_max();
    }
  }

  m_pWndProp->set_update_all();
  pDrawObj->invalidate_contours();
  pDrawObj->draw();
}

//---------------------------------------------------------------------------------------
// CTrackRangeCheckBox
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CTrackRangeCheckBox, CCheckBoxButton)

void CTrackRangeCheckBox::OnClickButton(CPoint point)
{
  CCheckBoxButton::OnClickButton(point);

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pDrawObj->get_colored_tracks().get_enable_user_range())
    pDrawObj->get_colored_tracks().restore_user_range();
  else
    pDrawObj->get_colored_tracks().get_min_max();

  m_pWndProp->set_update_all();
  pDrawObj->draw();
}

//---------------------------------------------------------------------------------------
// CCrossSectCheckBox
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CCrossSectCheckBox, CCheckBoxButton)

void CCrossSectCheckBox::OnClickButton(CPoint point)
{
  CCheckBoxButton::OnClickButton(point);

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  m_pWndProp->set_update_all();

  pDrawObj->invalidate_contours();
  pDrawObj->invalidate_hidden();

  pDrawObj->draw();
}

//---------------------------------------------------------------------------------------
// CSourceCheckBox
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CSourceCheckBox, CCheckBoxButton)

void CSourceCheckBox::OnClickButton(CPoint point)
{
  CCheckBoxButton::OnClickButton(point);

  EvaporatingParticle::CSource* pSrc = CParticleTrackingApp::Get()->GetTracker()->get_src();
  pSrc->invalidate();
}

//---------------------------------------------------------------------------------------
// CInvalidateFieldCheckBox
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CInvalidateFieldCheckBox, CCheckBoxButton)

void CInvalidateFieldCheckBox::OnClickButton(CPoint point)
{
  CCheckBoxButton::OnClickButton(point);

  EvaporatingParticle::CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
  int nCurrFieldId = pFields->get_curr_field_index();
  EvaporatingParticle::CElectricFieldData* pData = nCurrFieldId >= 0 ? pFields->at(nCurrFieldId) : NULL;
  if(pData != NULL)
    pData->invalidate();
}

//---------------------------------------------------------------------------------------
// CPlaneYZCalcCheckBox
//---------------------------------------------------------------------------------------
IMPLEMENT_DYNAMIC(CPlaneYZCalcCheckBox, CCheckBoxButton)

void CPlaneYZCalcCheckBox::OnClickButton(CPoint point)
{
  CCheckBoxButton::OnClickButton(point);

  EvaporatingParticle::CCalcCollection* pCalcColl = CParticleTrackingApp::Get()->GetCalcs();
  EvaporatingParticle::CPlaneYZCalculator* pPlaneYZCalc = NULL;
  size_t nCalcCount = pCalcColl->size();
  for(size_t i = 0; i < nCalcCount; i++)
  {
    if(pCalcColl->at(i)->get_enable_ptr() == m_dwData)
    {
      pPlaneYZCalc = (EvaporatingParticle::CPlaneYZCalculator*)(pCalcColl->at(i));
      break;
    }
  }

  if(pPlaneYZCalc != NULL)
  {
    pPlaneYZCalc->update();
    EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
    pDrawObj->invalidate_hidden();
    pDrawObj->draw();
  }
}
