
#include "stdafx.h"

#include "ResponseProperty.h"
#include "ParticleTracking.h"
#include "PropertiesWnd.h"

#include "ColorContour.h"

using namespace EvaporatingParticle;
//---------------------------------------------------------------------------------------
//  CResponseProperty
//---------------------------------------------------------------------------------------
BOOL CResponseProperty::OnUpdateValue()
{
  BOOL bRes = CMFCPropertyGridProperty::OnUpdateValue();

  DWORD_PTR pData = GetData();
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pDrawObj->get_draw_mode_ptr() == pData)
  {
    CString cMode = (CString)GetValue();
    for(int i = 0; i < EvaporatingParticle::CTrackDraw::dmCount; i++)
    {
      if(cMode == pDrawObj->get_draw_mode_name(i))
      {
        pDrawObj->set_draw_mode(i);
        pDrawObj->draw();
        break;
      }
    }
  }
  else if(pDrawObj->get_opacity_ptr() == pData)
  {
    pDrawObj->set_opacity(GetValue().dblVal);
    pDrawObj->draw();
  }
  else if(pDrawObj->get_cell_index_ptr() == pData)
  {
    pDrawObj->set_cell_index(GetValue().lVal);
    pDrawObj->draw();
  }
  else if(pDrawObj->get_colored_tracks().get_var_index_ptr() == pData)
  {
    CString cType = (CString)GetValue();
    for(UINT j = 0; j < EvaporatingParticle::CColorContour::varCount; j++)
    {
      if(cType == pDrawObj->get_colored_tracks().get_var_name(j))
      {
        pDrawObj->get_colored_tracks().set_var_index(j);
        m_pWndProp->set_update_all();
        pDrawObj->draw();
        break;
      }
    }
  }
  else if(pDrawObj->get_colored_tracks().get_scale_type_ptr() == pData)
  {
    CString cType = (CString)GetValue();
    for(UINT j = 0; j < EvaporatingParticle::CColorContour::varCount; j++)
    {
      if(cType == pDrawObj->get_colored_tracks().get_scale_name(j))
      {
        pDrawObj->get_colored_tracks().set_scale_type(j);
        m_pWndProp->set_update_all();
        pDrawObj->draw();
        break;
      }
    }
  }
  else if(pDrawObj->get_colored_tracks().get_color_map_type_ptr() == pData)
  {
    CString cType = (CString)GetValue();
    for(UINT j = 0; j < EvaporatingParticle::CColorContour::varCount; j++)
    {
      if(cType == pDrawObj->get_colored_tracks().get_clr_map_name(j))
      {
        pDrawObj->get_colored_tracks().set_color_map_type(j);
        m_pWndProp->set_update_all();
        pDrawObj->draw();
        break;
      }
    }
  }
  else if(pDrawObj->get_colored_tracks().get_levels_count_ptr() == pData)
  {
    pDrawObj->get_colored_tracks().set_levels_count(GetValue().lVal);
    pDrawObj->draw();
  }
// Range controls:
  else if(pDrawObj->get_colored_tracks().get_min_val_ptr() == pData)
  {
    pDrawObj->get_colored_tracks().set_min_val(GetValue().dblVal);
    pDrawObj->draw();
  }
  else if(pDrawObj->get_colored_tracks().get_max_val_ptr() == pData)
  {
    pDrawObj->get_colored_tracks().set_max_val(GetValue().dblVal);
    pDrawObj->draw();
  }
  else  // Contours kitchen. The changes in these controls must be processed immediately.
  {
    on_update_contour_ctrls();
  }

  return bRes;
}

void CResponseProperty::on_update_contour_ctrls()
{
  DWORD_PTR pData = GetData();
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  size_t nContoursCount = pDrawObj->get_contours_count();
  for(size_t i = 0; i < nContoursCount; i++)
  {
    EvaporatingParticle::CColorContour* pObj = pDrawObj->get_contour(i);
// Variable to be plotted:
    if(pObj->get_var_index_ptr() == pData)
    {
      CString cType = (CString)GetValue();
      for(UINT j = 0; j < EvaporatingParticle::CColorContour::varCount; j++)
      {
        if(cType == pObj->get_var_name(j))
        {
          pObj->set_var_index(j);
          m_pWndProp->set_update_all();
          pDrawObj->draw();
          return;
        }
      }
    }
// Scaling type (linear or logarithmic):
    if(pObj->get_scale_type_ptr() == pData)
    {
      CString cType = (CString)GetValue();
      for(UINT j = 0; j < EvaporatingParticle::CColorImage::scCount; j++)
      {
        if(cType == pObj->get_scale_name(j))
        {
          pObj->set_scale_type(j);
          pDrawObj->draw();
          return;
        }
      }
    }
// Color map type (Rainbow or Extended Rainbow):
    if(pObj->get_color_map_type_ptr() == pData)
    {
      CString cType = (CString)GetValue();
      for(UINT j = 0; j < EvaporatingParticle::CColorImage::cmCount; j++)
      {
        if(cType == pObj->get_clr_map_name(j))
        {
          pObj->set_color_map_type(j);
          pDrawObj->draw();
          return;
        }
      }
    }
// Count of contours:
    if(pObj->get_levels_count_ptr() == pData)
    {
      pObj->set_levels_count(GetValue().lVal);
      pDrawObj->draw();
      return;
    }
// Range controls:
    if(pObj->get_min_val_ptr() == pData)
    {
      pObj->set_min_val(GetValue().dblVal);
      pDrawObj->draw();
      return;
    }
    if(pObj->get_max_val_ptr() == pData)
    {
      pObj->set_max_val(GetValue().dblVal);
      pDrawObj->draw();
      return;
    }
  }
}

//---------------------------------------------------------------------------------------
//  CCrossSectionResponer is a response type of property for cross-sections.
//---------------------------------------------------------------------------------------
BOOL CCrossSectionResponer::OnUpdateValue()
{
  BOOL bRes = CMFCPropertyGridProperty::OnUpdateValue();

  m_pWndProp->set_data_to_model();
  DWORD_PTR pData = GetData();

  bool bNeedRedraw = false;
  EvaporatingParticle::CRegion* pReg = NULL;

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  EvaporatingParticle::CCrossSectColl* pPlanes = CParticleTrackingApp::Get()->GetPlanes();
  size_t nPlanesCount = pPlanes->size();
  for(size_t i = 0; i < nPlanesCount; i++)
  {
    EvaporatingParticle::CDomainCrossSection* pObj = pPlanes->at(i);
    pReg = pObj->get_region();
    if(pReg == NULL)
      continue;

    if(pData == pObj->get_enable_ptr())
    {
      bNeedRedraw = true;
      break;
    }

    if(pData == pObj->get_plane_type_ptr())
    {
      CString cType = (CString)GetValue();
      for(int j = 0; j < EvaporatingParticle::CDomainCrossSection::ptCount; j++)
      {
        if(cType == pObj->get_type_name(j))
        {
          bNeedRedraw = pObj->set_plane_type(j);
          pReg = pObj->get_region();  // the cross-section faces must be re-built.
          break;
        }
      }
      break;
    }

    if(pData == pObj->get_plane_origin_ptr())
    {
      int nPlaneType = pObj->get_plane_type();
      EvaporatingParticle::Vector3D vPlaneOrigin = pObj->get_plane_origin();
      switch(nPlaneType)
      {
        case EvaporatingParticle::CDomainCrossSection::ptPlaneXY: vPlaneOrigin.z = 0.1 * GetValue().dblVal; break;
        case EvaporatingParticle::CDomainCrossSection::ptPlaneXZ: vPlaneOrigin.y = 0.1 * GetValue().dblVal; break;
        case EvaporatingParticle::CDomainCrossSection::ptPlaneYZ: vPlaneOrigin.x = 0.1 * GetValue().dblVal; break;
      }

      bNeedRedraw = pObj->set_plane_origin(vPlaneOrigin);
      pReg = pObj->get_region();  // the cross-section faces must be re-built.
      break;
    }
  }

  if(bNeedRedraw)
  {
    m_pWndProp->set_update_all();
    pDrawObj->invalidate_contours();
    pDrawObj->invalidate_faces();
    pDrawObj->invalidate_aux();

    pDrawObj->draw();
  }

  return bRes;
}

//---------------------------------------------------------------------------------------
//  CCalcResponseProperty  - a response type of property specially used for calculators.
//---------------------------------------------------------------------------------------
BOOL CCalcResponseProperty::OnUpdateValue()
{
  BOOL bRes = CMFCPropertyGridProperty::OnUpdateValue();

  m_pWndProp->set_data_to_model();

  DWORD_PTR pData = GetData();
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  EvaporatingParticle::CCalcCollection* pCalcColl = CParticleTrackingApp::Get()->GetCalcs();
  size_t nCalcCount = pCalcColl->size();
  for(size_t i = 0; i < nCalcCount; i++)
  {
// Common properties of all types of calculators:
    EvaporatingParticle::CCalculator* pCalc = pCalcColl->at(i);

    if(pCalc->get_enable_ptr() == pData)
    {
      pCalc->set_enable(GetValue().boolVal);
      if(pCalc->get_update_flag())
      {
        pCalc->update();
        pCalc->do_calculate();
        pDrawObj->draw();
      }
    }
    else if(pCalc->get_clc_var_type_ptr() == pData)
    {
      CString cType = (CString)GetValue();
      for(UINT j = 0; j < pCalc->calc_vars_count(); j++)
      {
        if(cType == pCalc->get_var_name(j))
        {
          pCalc->set_clc_var_type(j);
          if(pCalc->get_update_flag())
            pCalc->do_calculate();

          if(pCalc->type() == EvaporatingParticle::CCalculator::ctPlaneYZ)
          {
            EvaporatingParticle::CPlaneYZCalculator* pPlaneYZCalc = (EvaporatingParticle::CPlaneYZCalculator*)pCalc;
            pPlaneYZCalc->update();
            pDrawObj->draw();
          }
        }
      }
    }

    int nType = pCalc->type();
    switch(nType)
    {
      case EvaporatingParticle::CCalculator::ctPlaneYZ:
      {
        EvaporatingParticle::CPlaneYZCalculator* pPlaneYZCalc = (EvaporatingParticle::CPlaneYZCalculator*)pCalc;
        if(pPlaneYZCalc->get_plane_x_ptr() == pData || pPlaneYZCalc->get_seq_clc_dir_ptr() == pData)
        {
          pPlaneYZCalc->update();
          pPlaneYZCalc->do_calculate();
          pDrawObj->draw();
        }

        break;
      }
      case EvaporatingParticle::CCalculator::ctTrackCalc:
      {
        EvaporatingParticle::CTrackCalculator* pTrackCalc = (EvaporatingParticle::CTrackCalculator*)pCalc;
        if(pTrackCalc->get_cs_pos_ptr() == pData)
        {
          pTrackCalc->set_cs_pos(0.1 * GetValue().dblVal);
          pTrackCalc->do_calculate();
        }

        break;
      }
      case EvaporatingParticle::CCalculator::ctAlongLine:
      {
        EvaporatingParticle::CLineCalculator* pLineCalc = (EvaporatingParticle::CLineCalculator*)pCalc;
// CLineCalculator::get_update_flag() returns false (and that is correct), so that the code above will not work for this calculator. 
// Here I want only the read-only "Average" edit-line to become zero if the variable type has been changed.
        if(pLineCalc->get_clc_var_type_ptr() == pData)
          pLineCalc->set_line_average(0);

        break;
      }
    }
  }

  m_pWndProp->set_update_all();
  return bRes;
}

//---------------------------------------------------------------------------------------
//  CColorResponseProperty
//---------------------------------------------------------------------------------------
BOOL CColorResponseProperty::OnUpdateValue()
{
  BOOL bRes = CMFCPropertyGridColorProperty::OnUpdateValue();

  DWORD_PTR pData = GetData();
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  if(pDrawObj->get_faces_color_ptr() == pData)
  {
    pDrawObj->set_faces_color(GetColor());
    pDrawObj->draw();
  }
  else if(pDrawObj->get_bkgr_color_ptr() == pData)
  {
    pDrawObj->set_bkgr_color(GetColor());
    pDrawObj->draw();
  }

  return bRes;
}

//---------------------------------------------------------------------------------------
//  CAreaFacesColorResponder
//---------------------------------------------------------------------------------------
BOOL CAreaFacesColorResponder::OnUpdateValue()
{
  BOOL bRes = CMFCPropertyGridColorProperty::OnUpdateValue();

  DWORD_PTR pData = GetData();
  EvaporatingParticle::CSelAreasColl* pSelAreasColl = CParticleTrackingApp::Get()->GetSelAreas();
  int nAreasCount = pSelAreasColl->size();
  for(int i = -1; i < nAreasCount; i++)
  {
    EvaporatingParticle::CSelectedAreas* pSelArea = (i >= 0) ? pSelAreasColl->at(i) : pSelAreasColl->get_default_area();
    if(pData == pSelArea->get_faces_color_ptr())
    {
      pSelArea->set_faces_color(GetColor());
      EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
      pDrawObj->draw();
      break;
    }
  }

  return bRes;
}

//---------------------------------------------------------------------------------------
//  CGeneralResponseProperty
//---------------------------------------------------------------------------------------
BOOL CGeneralResponseProperty::OnUpdateValue()
{
  BOOL bRes = CMFCPropertyGridProperty::OnUpdateValue();
  m_pWndProp->set_data_to_model();
  m_pWndProp->set_update_all();
  return bRes;
}

//---------------------------------------------------------------------------------------
//  CElectricFieldResponder
//---------------------------------------------------------------------------------------
BOOL CElectricFieldResponder::OnUpdateValue()
{
  BOOL bRes = CMFCPropertyGridProperty::OnUpdateValue();

  m_pWndProp->set_data_to_model();
  DWORD_PTR pData = GetData();

  CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
  pFields->set_curr_field_index(-1);
  if((DWORD_PTR)pFields == pData)
  {
    CString cFieldName = (CString)GetValue();
    size_t nFieldsCount = pFields->size();
    for(size_t i = 0; i < nFieldsCount; i++)
    {
      if(pFields->at(i)->get_field_name() == cFieldName)
      {
        pFields->set_curr_field_index(i);
        break;
      }
    }

    m_pWndProp->set_update_all();
  }

  return bRes;
}

//---------------------------------------------------------------------------------------
//  CNamedAreasSelResponder
//---------------------------------------------------------------------------------------
static CPotentialBoundCond* get_bc_object(DWORD_PTR pRegNames)
{
  CFieldDataColl* pFieldColl = CParticleTrackingApp::Get()->GetFields();
  size_t nFieldsCount = pFieldColl->size();
  for(size_t i = 0; i < nFieldsCount; i++)
  {
    CElectricFieldData* pField = pFieldColl->at(i);
    size_t nCount = pField->get_bc_count();
    for(size_t j = 0; j < nCount; j++)
    {
      CPotentialBoundCond* pBC = pField->get_bc(j);
      if(pRegNames == pBC->get_region_names_ptr())
        return pBC;
    }
  }

  return NULL;
}

static CAnalytRFField* get_analyt_rf_object(DWORD_PTR pRegNames)
{
  CFieldPtbCollection& coll = CParticleTrackingApp::Get()->GetTracker()->get_field_ptb();
  size_t nPtbCount = coll.size();
  for(size_t i = 0; i < nPtbCount; i++)
  {
    CFieldPerturbation* pPtb = coll.at(i);
    int nType = pPtb->type();
    if(nType == CFieldPerturbation::ptbFlatChannelRF || nType == CFieldPerturbation::ptbCylSubstrateRF || nType == CFieldPerturbation::ptbElliptSubstrRF)
    {
      CAnalytRFField* pAnalytPtb = (CAnalytRFField*)pPtb;
      if(pRegNames == pAnalytPtb->get_region_names_ptr())
        return pAnalytPtb;
    }
  }

  return NULL;
}

BOOL CNamedAreasSelResponder::OnUpdateValue()
{
  BOOL bRes = CMFCPropertyGridProperty::OnUpdateValue();

  m_pWndProp->set_data_to_model();

  CStringVector* pRegNames = (CStringVector*)GetData();
  CString sNamedArea = (CString)GetValue();   // user selected this string, it can be "None".

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();

  CSelAreasColl* pSelAreasColl = CParticleTrackingApp::Get()->GetSelAreas();
  size_t nSelAreasCount = pSelAreasColl->size();
  for(size_t i = 0; i < nSelAreasCount; i++)
  {
    CSelectedAreas* pSelArea = pSelAreasColl->at(i);
    CString sName = pSelArea->get_name();
    if(sName == sNamedArea)
    {
// Two possible types of objects who can posess the region names collection:
      CPotentialBoundCond* pBC = get_bc_object(GetData());
      CAnalytRFField* pAnalytPtb = get_analyt_rf_object(GetData());
      if(pBC != NULL)
      {
        pSelArea->merge_items(*pRegNames, pBC->get_merge_option());
        pDrawObj->enter_sel_context(pRegNames);
        pDrawObj->exit_sel_context(pRegNames);
        pBC->set_last_merged(sNamedArea);
        break;
      }
      else if(pAnalytPtb != NULL)
      {
        pSelArea->merge_items(*pRegNames, pAnalytPtb->get_merge_option());
        pDrawObj->enter_sel_context(pRegNames);
        pDrawObj->exit_sel_context(pRegNames);
        pAnalytPtb->set_last_merged(sNamedArea);
        break;
      }
      else
      {
        return FALSE;
      }
    }
  }

  m_pWndProp->set_update_all();
  pDrawObj->draw();
  return bRes;
}

//---------------------------------------------------------------------------------------
//  CSetAndRedrawResponder. 
//---------------------------------------------------------------------------------------
BOOL CSetAndRedrawResponder::OnUpdateValue()
{
  BOOL bRes = CGeneralResponseProperty::OnUpdateValue();
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->invalidate_contours();
  pDrawObj->draw();
  return bRes;
}




