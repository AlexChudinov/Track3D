
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "ResponseProperty.h"
#include "Button.h"

using namespace EvaporatingParticle;

void CPropertiesWnd::add_areas_ctrls()
{
  CMFCPropertyGridProperty* pMainAreaGroup = new CMFCPropertyGridProperty(_T("Named Areas (Sets of 2D Regions)"));
  CAddNamedAreaButton* pAddBtn = new CAddNamedAreaButton(this, _T("Create New Named Area"), _T(""), _T("Click this button to add a new instance of named areas. You will have to select regions to fill this object with real contents."), NULL);
  pMainAreaGroup->AddSubItem(pAddBtn);
  m_wndPropList.AddProperty(pMainAreaGroup);

  CSelAreasColl* pSelAreasColl = CParticleTrackingApp::Get()->GetSelAreas();
  int nSelAreasCount = pSelAreasColl->size();
  char buff[8];
  for(int i = -1; i < nSelAreasCount; i++)
  {
    bool bSelectedArea = (i >= 0);
// Note: the default area will be populated after pSelAreasColl->get_default_area() call.
    CSelectedAreas* pSelAreas = bSelectedArea ? pSelAreasColl->at(i) : pSelAreasColl->get_default_area();
    CString sGroupName = bSelectedArea ? CString(_T("Named Regions Set #")) + CString(itoa(i + 1, buff, 10)) : CString(_T("Default Area"));

// Note: the GetData() of this group control will return, in fact, the pointer to pSelAreas. This will be used in CHideShowRegsCheckBox::OnClickButton().
    CMFCPropertyGridProperty* pNamedAreasGroup = new CMFCPropertyGridProperty(sGroupName, (DWORD_PTR)pSelAreas);

    CMFCPropertyGridProperty* pNameProp = new CMFCPropertyGridProperty(_T("Name"), pSelAreas->get_name(), _T("The name of the whole selection. This name will appear in the drop-down lists."), pSelAreas->get_name_ptr());
    pNameProp->AllowEdit(bSelectedArea);
    pNamedAreasGroup->AddSubItem(pNameProp);

    CString cRegNames = EvaporatingParticle::CObject::compile_string(*pSelAreas);
    if(pSelAreas->size() != 0)
      cRegNames += CString("  ") + pSelAreas->get_square_str() + CString(" mm^2");

    CSelectRegionButton* pSelRegButton = new CSelectRegionButton(this, _T("Select 2D Regions"), cRegNames, _T("Click to select 2D regions in the main view window for boundary conditions."), (DWORD_PTR)pSelAreas);
    pSelRegButton->Enable(bSelectedArea);
    pNamedAreasGroup->AddSubItem(pSelRegButton);
    
    COLORREF clr = pSelAreas->get_faces_color();
    CAreaFacesColorResponder* pColorBtn = new CAreaFacesColorResponder(_T("Area Faces Color"), clr, NULL, _T(""), pSelAreas->get_faces_color_ptr());
    pNamedAreasGroup->AddSubItem(pColorBtn);

    CHideShowRegsCheckBox* pHideShowBtn = new CHideShowRegsCheckBox(this, _T("Visibility"), (_variant_t)pSelAreas->get_visibility_flag(), _T("Change the visibility status of the selected regions"), pSelAreas->get_visibility_flag_ptr());
    pNamedAreasGroup->AddSubItem(pHideShowBtn);

    if(bSelectedArea)
    {
      CRemoveNamedAreaButton* pRemBtn = new CRemoveNamedAreaButton(this, _T("Remove Set of 2D Regions"), _T(""), _T("Click here to remove this set of selected regions. Note that removing a named set of regions will NOT affect any settings anywhere else."), (DWORD_PTR)pSelAreas);
      pNamedAreasGroup->AddSubItem(pRemBtn);
    }

    m_wndPropList.AddProperty(pNamedAreasGroup);
  }
}

void CPropertiesWnd::set_areas_data()
{
  CSelAreasColl* pSelAreasColl = CParticleTrackingApp::Get()->GetSelAreas();
  size_t nSelAreasCount = pSelAreasColl->size();
  for(size_t i = 0; i < nSelAreasCount; i++)
  {
    CSelectedAreas* pSelAreas = pSelAreasColl->at(i);
    CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pSelAreas->get_name_ptr());
    if(pProp != NULL)
    {
      CString sName = (CString)pProp->GetValue();
      pSelAreas->set_name(sName);
    }
  }
}

void CPropertiesWnd::update_areas_ctrls()
{
}
