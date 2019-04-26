
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
  m_pWndProp->set_data_to_model();

  BOOL bRes = CMFCPropertyGridProperty::OnUpdateValue();

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
        if(pPlaneYZCalc->get_plane_x_ptr() == pData)
        {
          pPlaneYZCalc->set_plane_x(0.1 * GetValue().dblVal);
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
  m_pWndProp->set_data_to_model();

  BOOL bRes = CMFCPropertyGridProperty::OnUpdateValue();
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

