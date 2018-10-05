
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "ResponseProperty.h"
#include "ColorContour.h"
#include "Button.h"

//---------------------------------------------------------------------------------------
// Drawing
//---------------------------------------------------------------------------------------
void CPropertiesWnd::add_draw_ctrls()
{
  EvaporatingParticle::CTrackDraw* pObj = CParticleTrackingApp::Get()->GetDrawObj();

// Drawing:
//  CMFCPropertyGridProperty* pDrawGroup = new CMFCPropertyGridProperty(_T("Drawing"));
  CMFCPropertyGridProperty* pFacesGroup = new CMFCPropertyGridProperty(_T("Geometry"));
// Geometry:
  COleVariant var(pObj->get_draw_mode_name(pObj->get_draw_mode()));
  CResponseProperty* pRespProp = new CResponseProperty(this, _T("Drawing Mode"), var, _T("Select the drowing mode."), pObj->get_draw_mode_ptr());
  for(int k = 0; k < EvaporatingParticle::CTrackDraw::dmCount; k++)
    pRespProp->AddOption(pObj->get_draw_mode_name(k));

  pRespProp->AllowEdit(FALSE);
  pFacesGroup->AddSubItem(pRespProp);

// Color and opacity
  COLORREF clr = pObj->get_faces_color();
  CColorResponseProperty* pColorBtn = new CColorResponseProperty(_T("Faces Color"), clr, NULL, _T(""), pObj->get_faces_color_ptr());
  pFacesGroup->AddSubItem(pColorBtn);

  pRespProp = new CResponseProperty(this, _T("Faces Opacity"), COleVariant(pObj->get_opacity()), _T("Opacity of the faces material. Unity gives opaque faces, zero means transparent faces"), pObj->get_opacity_ptr());
  pFacesGroup->AddSubItem(pRespProp);

  CString cHiddenRegNames = pObj->get_hidden_names_str();
  CSelectRegionButton* pHiddenRegButton = new CSelectRegionButton(this, _T("Hidden 2D Regions"), cHiddenRegNames, _T("Click to select 2D regions to hide."), pObj->get_hidden_reg_names_ptr());
  pFacesGroup->AddSubItem(pHiddenRegButton);

  CSelectRegionButton* pShowRegButton = new CSelectRegionButton(this, _T("Show All 2D Regions"), _T(""), _T("Click to show all 2D regions."), NULL);
  pFacesGroup->AddSubItem(pShowRegButton);

  m_wndPropList.AddProperty(pFacesGroup);

// Tracks:
  CMFCPropertyGridProperty* pTracksGroup = new CMFCPropertyGridProperty(_T("Tracks"));

  CRedrawCheckBox* pCheckBox = new CRedrawCheckBox(this, _T("Enable Tracks"), (_variant_t)pObj->get_enable_tracks(), _T("Turns ON/OFF drawing trajectories of particles"), pObj->get_enable_tracks_ptr());
  pTracksGroup->AddSubItem(pCheckBox);

  EvaporatingParticle::CColoredTracks& color_tracks = pObj->get_colored_tracks();
// Track color map:
  COleVariant var1(color_tracks.get_var_name(color_tracks.get_var_index()));
  CResponseProperty* pVarSelectProp = new CResponseProperty(this, _T("Variable"), var1, _T("Specify a variable the tracks will be colored by."), color_tracks.get_var_index_ptr());
  for(size_t j = 0; j < EvaporatingParticle::CColoredTracks::tcmCount; j++)
    pVarSelectProp->AddOption(color_tracks.get_var_name(j));

  pVarSelectProp->AllowEdit(FALSE);
  pTracksGroup->AddSubItem(pVarSelectProp);

// Scaling type (linear or logarithmic):
  COleVariant var2(color_tracks.get_scale_name(color_tracks.get_scale_type()));
  CResponseProperty* pScaleSelectProp = new CResponseProperty(this, _T("Scale Type"), var2, _T("Specify the color scaling to be either linear or logarithmic."), color_tracks.get_scale_type_ptr());
  for(size_t k = 0; k < EvaporatingParticle::CColorImage::scCount; k++)
    pScaleSelectProp->AddOption(color_tracks.get_scale_name(k));

  pScaleSelectProp->AllowEdit(FALSE);
  pTracksGroup->AddSubItem(pScaleSelectProp);

// Color map type (Rainbow or Extended Rainbow):
  COleVariant var3(color_tracks.get_clr_map_name(color_tracks.get_color_map_type()));
  CResponseProperty* pColorMapSelectProp = new CResponseProperty(this, _T("Color Map Type"), var3, _T("Specify the color mapping type."), color_tracks.get_color_map_type_ptr());
  for(size_t l = 0; l < EvaporatingParticle::CColorImage::cmCount; l++)
    pColorMapSelectProp->AddOption(color_tracks.get_clr_map_name(l));

  pColorMapSelectProp->AllowEdit(FALSE);
  pTracksGroup->AddSubItem(pColorMapSelectProp);

// Count of contours:
  CResponseProperty* pCountProp = new CResponseProperty(this, _T("Count of Levels"), COleVariant((long)color_tracks.get_levels_count()), _T("Total count of color levels."), color_tracks.get_levels_count_ptr());
  pTracksGroup->AddSubItem(pCountProp);

// Colored range:
  CTrackRangeCheckBox* pUserRangeProp = new CTrackRangeCheckBox(this, _T("User-Defined Range"), (_variant_t)color_tracks.get_enable_user_range(), _T("If this is ON the colored range is defined by the user."), color_tracks.get_enable_user_range_ptr());
  pTracksGroup->AddSubItem(pUserRangeProp);

  CResponseProperty* pMinValProp = new CResponseProperty(this, _T("Min (SI Units)"), COleVariant(color_tracks.get_min_val()), _T("Specify this to define the low margin of the colored range."), color_tracks.get_min_val_ptr());
  pTracksGroup->AddSubItem(pMinValProp);

  CResponseProperty* pMaxValProp = new CResponseProperty(this, _T("Max (SI Units)"), COleVariant(color_tracks.get_max_val()), _T("Specify this to define the upper margin of the colored range."), color_tracks.get_max_val_ptr());
  pTracksGroup->AddSubItem(pMaxValProp);

  m_wndPropList.AddProperty(pTracksGroup);

// General:
  CMFCPropertyGridProperty* pGeneralGroup = new CMFCPropertyGridProperty(_T("General"));

  clr = pObj->get_bkgr_color();
  CColorResponseProperty* pBkgColorBtn = new CColorResponseProperty(_T("Background Color"), clr, NULL, _T(""), pObj->get_bkgr_color_ptr());
  pGeneralGroup->AddSubItem(pBkgColorBtn);

  CCheckBoxButton* pGenericCheckBox = new CCheckBoxButton(this, _T("Rotate around Center"), (_variant_t)pObj->get_rot_center(), _T("If true the rotation takes place around the bounding box center; otherwise - around coordinate origin."), pObj->get_rot_center_ptr());
  pGeneralGroup->AddSubItem(pGenericCheckBox);

  pCheckBox = new CRedrawCheckBox(this, _T("Visualize Dirichlet Cells"), (_variant_t)pObj->get_enable_draw_norm(), _T("Enables / disables wireframe drawing a specified Dirichlet cell."), pObj->get_enable_draw_norm_ptr());
  pGeneralGroup->AddSubItem(pCheckBox);

  CResponseProperty* pIndexProp = new CResponseProperty(this, _T("Index of Dirichlet Cell"), COleVariant((long)pObj->get_cell_index()), _T("Specify node index for drawing the Dirichlet cell."), pObj->get_cell_index_ptr());
  pGeneralGroup->AddSubItem(pIndexProp);

  m_wndPropList.AddProperty(pGeneralGroup);

// Cross-section planes:
  add_cs_plane_ctrls();

// Contours:
  add_contour_ctrls();
}

void CPropertiesWnd::set_draw_data()
{
}

void CPropertiesWnd::add_cs_plane_ctrls()
{
  CMFCPropertyGridProperty* pCrossSectionGroup = new CMFCPropertyGridProperty(_T("Cross-Section Planes"));

  char buff[4];
  EvaporatingParticle::CCrossSectColl* pPlanes = CParticleTrackingApp::Get()->GetPlanes();
  size_t nPlanesCount = pPlanes->size();
  for(size_t i = 0; i < nPlanesCount; i++)
  {
    EvaporatingParticle::CDomainCrossSection* pObj = pPlanes->at(i);
    CString cName = CString(pObj->get_name().c_str());
    CMFCPropertyGridProperty* pCrossSectProp = new CMFCPropertyGridProperty(cName);

// Enable cross-section:
    CCrossSectCheckBox* pCheckBox = new CCrossSectCheckBox(this, _T("Enable Plane"), (_variant_t)pObj->get_enable(), _T("Turns ON/OFF drawing the cross-section plane"), pObj->get_enable_ptr());
    pCrossSectProp->AddSubItem(pCheckBox);

// Type of the plane:
    int nPlaneType = pObj->get_plane_type();
    COleVariant var(pObj->get_type_name(nPlaneType));
    CCrossSectionResponer* pTypeProp = new CCrossSectionResponer(this, _T("Type of the Plane"), var, _T("Specify the type of the cross-section plane."), pObj->get_plane_type_ptr());
    for(size_t j = 0; j < EvaporatingParticle::CDomainCrossSection::ptCount; j++)
      pTypeProp->AddOption(pObj->get_type_name(j));

    pTypeProp->AllowEdit(FALSE);
    pCrossSectProp->AddSubItem(pTypeProp);

    CCrossSectionResponer* pRespProp = NULL;
// Plane origin and normal:
    CMFCPropertyGridProperty* pPosGroup = new CMFCPropertyGridProperty(_T("Position"), pObj->get_plane_origin_ptr());
    switch(nPlaneType)
    {
      case EvaporatingParticle::CDomainCrossSection::ptPlaneXY:
      {
        pRespProp = new CCrossSectionResponer(this, _T("Z, mm"), COleVariant(10 * pObj->get_plane_origin().z), _T("Specify Z-coordinate of the XY plane."), pObj->get_plane_origin_ptr());
        pPosGroup->AddSubItem(pRespProp);
        break;
      }
      case EvaporatingParticle::CDomainCrossSection::ptPlaneXZ:
      {
        pRespProp = new CCrossSectionResponer(this, _T("Y, mm"), COleVariant(10 * pObj->get_plane_origin().y), _T("Specify Y-coordinate of the XZ plane."), pObj->get_plane_origin_ptr());
        pPosGroup->AddSubItem(pRespProp);
        break;
      }
      case EvaporatingParticle::CDomainCrossSection::ptPlaneYZ:
      {
        pRespProp = new CCrossSectionResponer(this, _T("X, mm"), COleVariant(10 * pObj->get_plane_origin().x), _T("Specify X-coordinate of the YZ plane."), pObj->get_plane_origin_ptr());
        pPosGroup->AddSubItem(pRespProp);
        break;
      }
    }

    pCrossSectProp->AddSubItem(pPosGroup);

// Remove cross-section button:
    COleVariant var1(_T(" "));
    CRemoveCrossSectionButton* pRemButton = new CRemoveCrossSectionButton(this, _T("Remove Cross-Section"), var1, _T("Click this to remove this cross-section plane from the scene."), (DWORD_PTR)pObj);
    pCrossSectProp->AddSubItem(pRemButton);

    pCrossSectionGroup->AddSubItem(pCrossSectProp);
  }

  m_wndPropList.AddProperty(pCrossSectionGroup);
}

void CPropertiesWnd::add_contour_ctrls()
{
  CMFCPropertyGridProperty* pContourGroup = new CMFCPropertyGridProperty(_T("Contours"));

  char buff[4];
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  size_t nCountourCount = pDrawObj->get_contours_count();
  for(size_t i = 0; i < nCountourCount; i++)
  {
    CString cName = CString("Contour # ") + CString(itoa(i + 1, buff, 10));
    CMFCPropertyGridProperty* pContourProp = new CMFCPropertyGridProperty(cName);

    EvaporatingParticle::CColorContour* pObj = pDrawObj->get_contour(i);

// Enable flag:
    CRedrawCheckBox* pRespProp = new CRedrawCheckBox(this, _T("Enable"), (_variant_t)pObj->get_enable_image(), _T("Turns ON/OFF drawing the contour"), pObj->get_enable_image_ptr());
    pContourProp->AddSubItem(pRespProp);

// Locations:
    CString cRegNames = pObj->get_drawn_reg_names();
    CSelectRegionButton* pSelectRegButton = new CSelectRegionButton(this, _T("Locations"), cRegNames, _T("Click to select 2D region(s) at which this contour will be plotted."), pObj->get_drawn_reg_names_ptr());
    pContourProp->AddSubItem(pSelectRegButton);

// Clear all locations:
    CClearLocationsButton* pClearLocBtn = new CClearLocationsButton(this, _T("Clear Locations"), _T(""), _T("Click to clear all locations from this contour."), (DWORD_PTR)pObj);
    pContourProp->AddSubItem(pClearLocBtn);

// Variable to be plotted:
    COleVariant var1(pObj->get_var_name(pObj->get_var_index()));
    CResponseProperty* pVarSelectProp = new CResponseProperty(this, _T("Variable"), var1, _T("Specify a variable to be contoured."), pObj->get_var_index_ptr());
    for(size_t j = 0; j < EvaporatingParticle::CColorContour::varCount; j++)
      pVarSelectProp->AddOption(pObj->get_var_name(j));

    pVarSelectProp->AllowEdit(FALSE);
    pContourProp->AddSubItem(pVarSelectProp);

// Scaling type (linear or logarithmic):
    COleVariant var2(pObj->get_scale_name(pObj->get_scale_type()));
    CResponseProperty* pScaleSelectProp = new CResponseProperty(this, _T("Scale Type"), var2, _T("Specify the color scaling to be either linear or logarithmic."), pObj->get_scale_type_ptr());
    for(size_t k = 0; k < EvaporatingParticle::CColorImage::scCount; k++)
      pScaleSelectProp->AddOption(pObj->get_scale_name(k));

    pScaleSelectProp->AllowEdit(FALSE);
    pContourProp->AddSubItem(pScaleSelectProp);

// Color map type (Rainbow or Extended Rainbow):
    COleVariant var3(pObj->get_clr_map_name(pObj->get_color_map_type()));
    CResponseProperty* pColorMapSelectProp = new CResponseProperty(this, _T("Color Map Type"), var3, _T("Specify the color mapping."), pObj->get_color_map_type_ptr());
    for(size_t l = 0; l < EvaporatingParticle::CColorImage::cmCount; l++)
      pColorMapSelectProp->AddOption(pObj->get_clr_map_name(l));

    pColorMapSelectProp->AllowEdit(FALSE);
    pContourProp->AddSubItem(pColorMapSelectProp);

// Count of contours:
    CResponseProperty* pCountProp = new CResponseProperty(this, _T("Count of Levels"), COleVariant((long)pObj->get_levels_count()), _T("Total count of levels."), pObj->get_levels_count_ptr());
    pContourProp->AddSubItem(pCountProp);

// Colored range:
    CContourRangeCheckBox* pUserRangeProp = new CContourRangeCheckBox(this, _T("User-Defined Range"), (_variant_t)pObj->get_enable_user_range(), _T("If this is ON the colored range is defined by the user."), pObj->get_enable_user_range_ptr());
    pContourProp->AddSubItem(pUserRangeProp);

    CResponseProperty* pMinValProp = new CResponseProperty(this, _T("Min (SI Units)"), COleVariant(pObj->get_min_val()), _T("Specify this to define the low margin of the colored range."), pObj->get_min_val_ptr());
    pContourProp->AddSubItem(pMinValProp);

    CResponseProperty* pMaxValProp = new CResponseProperty(this, _T("Max (SI Units)"), COleVariant(pObj->get_max_val()), _T("Specify this to define the upper margin of the colored range."), pObj->get_max_val_ptr());
    pContourProp->AddSubItem(pMaxValProp);

// Enable drawing contour lines:
    CRedrawCheckBox* pDrawLinesProp = new CRedrawCheckBox(this, _T("Draw Contour Lines"), (_variant_t)pObj->get_enable_lines(), _T("Turns ON/OFF drawing the contour lines"), pObj->get_enable_lines_ptr());
    pContourProp->AddSubItem(pDrawLinesProp);

// Remove contour button.
    COleVariant var4(_T(""));
    CRemoveContourButton* pRemoveButton = new CRemoveContourButton(this, _T("Remove Contour"), var4, _T("Click to remove this contour."), (DWORD_PTR)pObj);
    pContourProp->AddSubItem(pRemoveButton);

    pContourGroup->AddSubItem(pContourProp);
  }

  m_wndPropList.AddProperty(pContourGroup);
}

void CPropertiesWnd::update_draw_ctrls()
{
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  bool bEnableClr = pDrawObj->get_draw_mode() == EvaporatingParticle::CTrackDraw::dmFlatAndWire || pDrawObj->get_draw_mode() == EvaporatingParticle::CTrackDraw::dmFlatOnly;
  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pDrawObj->get_faces_color_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableClr);

  pProp = m_wndPropList.FindItemByData(pDrawObj->get_opacity_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableClr);

  EvaporatingParticle::CColoredTracks* pObj = &(CParticleTrackingApp::Get()->GetDrawObj()->get_colored_tracks());

  bool bEnable = false;
  pProp = m_wndPropList.FindItemByData(pObj->get_enable_user_range_ptr());
  if(pProp != NULL)
    bEnable = pProp->GetValue().boolVal;

  pProp = m_wndPropList.FindItemByData(pObj->get_min_val_ptr());
  if(pProp != NULL)
  {
    pProp->Enable(bEnable);
    if(!bEnable)
      pProp->SetValue(COleVariant(pObj->get_min_val()));
  }

  pProp = m_wndPropList.FindItemByData(pObj->get_max_val_ptr());
  if(pProp != NULL)
  {
    pProp->Enable(bEnable);
    if(!bEnable)
      pProp->SetValue(COleVariant(pObj->get_max_val()));
  }

  pProp = m_wndPropList.FindItemByData(pDrawObj->get_cell_index_ptr());
  if(pProp != NULL)
    pProp->Enable(pDrawObj->get_enable_draw_norm());

  update_contour_ctrls();
}

void CPropertiesWnd::update_contour_ctrls()
{
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  size_t nCountourCount = pDrawObj->get_contours_count();
  for(size_t i = 0; i < nCountourCount; i++)
  {
    EvaporatingParticle::CColorContour* pObj = pDrawObj->get_contour(i);

    bool bEnable = false;
    CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pObj->get_enable_user_range_ptr());
    if(pProp != NULL)
      bEnable = pProp->GetValue().boolVal;

    pProp = m_wndPropList.FindItemByData(pObj->get_min_val_ptr());
    if(pProp != NULL)
    {
      pProp->Enable(bEnable);
      if(!bEnable)
        pProp->SetValue(COleVariant(pObj->get_min_val()));
    }

    pProp = m_wndPropList.FindItemByData(pObj->get_max_val_ptr());
    if(pProp != NULL)
    {
      pProp->Enable(bEnable);
      if(!bEnable)
        pProp->SetValue(COleVariant(pObj->get_max_val()));
    }
  }
}
