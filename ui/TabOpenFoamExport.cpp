
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "ExportOpenFOAM.h"
#include "ResponseProperty.h"
#include "Button.h"

//---------------------------------------------------------------------------------------
// OpenFOAM export
//---------------------------------------------------------------------------------------
void CPropertiesWnd::add_export_ctrls()
{
// OpenFOAM export:
  EvaporatingParticle::CExportOpenFOAM* pExpObj = CParticleTrackingApp::Get()->GetExporter();

  CMFCPropertyGridProperty* pOpenFoamGroup = new CMFCPropertyGridProperty(_T("Export Mesh in OpenFOAM Format"));

  CCheckBoxButton* pEnableBoundCond = new CCheckBoxButton(this, _T("Export Boundary Conditions"), (_variant_t)pExpObj->get_enable_bound_cond(), _T("Turns ON/OFF boundary conditions export."), pExpObj->get_enable_bound_cond_ptr());
  pOpenFoamGroup->AddSubItem(pEnableBoundCond);

  CCheckBoxButton* pExportInternal = new CCheckBoxButton(this, _T("Export Internal Data"), (_variant_t)pExpObj->get_export_internal(), _T("When this flag is true the internal data from ANSYS CFX solution are exported."), pExpObj->get_export_internal_ptr());
  pOpenFoamGroup->AddSubItem(pExportInternal);

  CMFCPropertyGridProperty* pShiftProp = new CMFCPropertyGridProperty(_T("X-Coordinate Shift, mm"), COleVariant(pExpObj->get_coord_shift()), _T("Defines the shift of the exported part of the mesh with respect to the whole mesh."), pExpObj->get_coord_shift_ptr());
  pOpenFoamGroup->AddSubItem(pShiftProp);

  CMFCPropertyGridProperty* pDefParamGroup = new CMFCPropertyGridProperty(_T("Default Values"));

  CMFCPropertyGridProperty* pDefPressProp = new CMFCPropertyGridProperty(_T("Pressure, Pa"), COleVariant(pExpObj->get_def_press()), _T("Default value of pressure if no data in the data file are available."), pExpObj->get_def_press_ptr());
  pDefParamGroup->AddSubItem(pDefPressProp);

  CMFCPropertyGridProperty* pDefTempProp = new CMFCPropertyGridProperty(_T("Temperature, K"), COleVariant(pExpObj->get_def_temp()), _T("Default value of temperature if no data in the data file are available."), pExpObj->get_def_temp_ptr());
  pDefParamGroup->AddSubItem(pDefTempProp);

  CMFCPropertyGridProperty* pDefVxProp = new CMFCPropertyGridProperty(_T("Vx, m/s"), COleVariant(pExpObj->get_def_vx()), _T("Default value of axial velocity if no data in the data file are available."), pExpObj->get_def_vx_ptr());
  pDefParamGroup->AddSubItem(pDefVxProp);

  pOpenFoamGroup->AddSubItem(pDefParamGroup);

// Dynamic creation of boundary conditions:
  add_bc_ctrls(pOpenFoamGroup);

// Boundary data filename:
  CMFCPropertyGridProperty* pDataGroup = new CMFCPropertyGridProperty(_T("Data Files"));
  static TCHAR BASED_CODE szFilter[] = _T("Geometry Data Files(*.geom)|*.geom|All Files(*.*)|*.*||");
  CMFCPropertyGridFileProperty* pBoundCondFileProp = new CMFCPropertyGridFileProperty(_T("Boundary Data File"), TRUE, pExpObj->get_boundary_cond_file(), _T("geom"),
    OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, _T("Specify location of the boundary conditions data file for OpenFOAM export."), pExpObj->get_boundary_cond_file_ptr());

  pDataGroup->AddSubItem(pBoundCondFileProp);
  pOpenFoamGroup->AddSubItem(pDataGroup);

  m_wndPropList.AddProperty(pOpenFoamGroup);
}

static const size_t nTypeCount = 4;
static const char* cBoundCondType[nTypeCount] = { "Wall", "Inlet", "Patch", "Symmetry" };

void CPropertiesWnd::add_bc_ctrls(CMFCPropertyGridProperty* pOpenFOAMGroup)
{
  CMFCPropertyGridProperty* pBoundCondGroup = new CMFCPropertyGridProperty(_T("Boundary Conditions"));

  char buff[4];
  EvaporatingParticle::CExportOpenFOAM* pExpObj = CParticleTrackingApp::Get()->GetExporter();
  size_t nBoundCondCount = pExpObj->get_bc_count();
  for(size_t i = 0; i < nBoundCondCount; i++)
  {
    CString cName = CString("Boundary Conditions # ") + CString(itoa(i + 1, buff, 10));
    CMFCPropertyGridProperty* pBoundCondProp = new CMFCPropertyGridProperty(cName);

    EvaporatingParticle::CBoundaryConditions* pBC = pExpObj->get_bound_cond(i);
// Boundary conditions type:
    COleVariant var(_T(cBoundCondType[pBC->nType]));
    CGeneralResponseProperty* pType = new CGeneralResponseProperty(this, _T("Type"), var, _T("Specify type of boundary conditions: Wall, Inlet, Patch or Symmetry."), (DWORD_PTR)&(pBC->nType));
    for(size_t j = 0; j < nTypeCount; j++)
      pType->AddOption(_T(cBoundCondType[j]));

    pType->AllowEdit(FALSE);
    pBoundCondProp->AddSubItem(pType);

// Pressure:
    CMFCPropertyGridProperty* pPressProp = new CMFCPropertyGridProperty(_T("Pressure, Pa"), COleVariant(pBC->fPress), _T("Use this edit control to set pressure at the boundary."), (DWORD_PTR)&(pBC->fPress));
    pBoundCondProp->AddSubItem(pPressProp);

// Temperature:
    CMFCPropertyGridProperty* pTempProp = new CMFCPropertyGridProperty(_T("Temperature, K"), COleVariant(pBC->fTemp), _T("Use this edit control to set temperature at the boundary."), (DWORD_PTR)&(pBC->fTemp));
    pBoundCondProp->AddSubItem(pTempProp);

// 2D regions selector:
    CString cRegNames = EvaporatingParticle::CObject::compile_string(pBC->vRegNames);
    CSelectRegionButton* pSelectRegButton = new CSelectRegionButton(this, _T("2D Regions"), cRegNames, _T("Click to select 2D regions at which these boundary conditions are set."), (DWORD_PTR)&(pBC->vRegNames));
    pBoundCondProp->AddSubItem(pSelectRegButton);

// Remove boundary conditions button:
    COleVariant var1(_T(""));
    CRemoveBoundCondButton* pRemoveBCBtn = new CRemoveBoundCondButton(this, _T("Remove Boundary Condition"), var1, _T("Click this button to remove the boundary condition."), (DWORD_PTR)pBC);
    pBoundCondProp->AddSubItem(pRemoveBCBtn);

    pBoundCondGroup->AddSubItem(pBoundCondProp);
  }

  pOpenFOAMGroup->AddSubItem(pBoundCondGroup);
}

void CPropertiesWnd::set_export_data()
{
  EvaporatingParticle::CExportOpenFOAM* pExportObj = CParticleTrackingApp::Get()->GetExporter();

  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pExportObj->get_enable_bound_cond_ptr());
  if(pProp != NULL)
    pExportObj->set_enable_bound_cond(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(pExportObj->get_export_internal_ptr());
  if(pProp != NULL)
    pExportObj->set_export_internal(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(pExportObj->get_coord_shift_ptr());
  if(pProp != NULL)
    pExportObj->set_coord_shift(pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pExportObj->get_def_press_ptr());
  if(pProp != NULL)
    pExportObj->set_def_press(pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pExportObj->get_def_temp_ptr());
  if(pProp != NULL)
    pExportObj->set_def_temp(pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pExportObj->get_def_vx_ptr());
  if(pProp != NULL)
    pExportObj->set_def_vx(pProp->GetValue().dblVal);

  set_bc_data();

  pProp = m_wndPropList.FindItemByData(pExportObj->get_boundary_cond_file_ptr());
  if(pProp != NULL)
  {
    CString cFile = (CString)pProp->GetValue();
    std::string str = std::string(CT2CA(cFile));
    pExportObj->set_boundary_cond_file(str.c_str());
  }
}

void CPropertiesWnd::set_bc_data()
{
  EvaporatingParticle::CExportOpenFOAM* pExpObj = CParticleTrackingApp::Get()->GetExporter();
  size_t nBoundCondCount = pExpObj->get_bc_count();
  for(size_t i = 0; i < nBoundCondCount; i++)
  {
    EvaporatingParticle::CBoundaryConditions* pBC = pExpObj->get_bound_cond(i);

    DWORD_PTR pData = (DWORD_PTR)&(pBC->nType);
    CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pData);
    if(pProp != NULL)
    {
      CString cType = (CString)pProp->GetValue();
      for(size_t j = 0; j < nTypeCount; j++)
      {
        if(cType == cBoundCondType[j])
        {
          pBC->nType = j;
          break;
        }
      }
    }

    pData = (DWORD_PTR)&(pBC->fPress);
    pProp = m_wndPropList.FindItemByData(pData);
    if(pProp != NULL)
      pBC->fPress = pProp->GetValue().dblVal;

    pData = (DWORD_PTR)&(pBC->fTemp);
    pProp = m_wndPropList.FindItemByData(pData);
    if(pProp != NULL)
      pBC->fTemp = pProp->GetValue().dblVal;
  }
}

void CPropertiesWnd::update_export_ctrls()
{
  bool bMeshReady = !m_bBusy && CParticleTrackingApp::Get()->GetTracker()->is_ready();

  EvaporatingParticle::CExportOpenFOAM* pExportObj = CParticleTrackingApp::Get()->GetExporter();

  bool bEnable = false;
  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pExportObj->get_enable_bound_cond_ptr());
  if(pProp != NULL)
    bEnable = bMeshReady && pProp->GetValue().boolVal;

  pProp = m_wndPropList.FindItemByData(pExportObj->get_export_internal_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(pExportObj->get_coord_shift_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(pExportObj->get_def_press_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(pExportObj->get_def_temp_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(pExportObj->get_def_vx_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  update_bc_ctrls(bEnable);
  
  pProp = m_wndPropList.FindItemByData(pExportObj->get_boundary_cond_file_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);
}

void CPropertiesWnd::update_bc_ctrls(bool bEnable)
{
  bool bPressEnabled = false, bTempEnabled = false;
  EvaporatingParticle::CExportOpenFOAM* pExpObj = CParticleTrackingApp::Get()->GetExporter();
  size_t nBoundCondCount = pExpObj->get_bc_count();
  for(size_t i = 0; i < nBoundCondCount; i++)
  {
    EvaporatingParticle::CBoundaryConditions* pBC = pExpObj->get_bound_cond(i);

    DWORD_PTR pData = (DWORD_PTR)&(pBC->nType);
    CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pData);
    if(pProp != NULL)
    {
      pProp->Enable(bEnable);
      if(bEnable)
      {
        CString cType = (CString)pProp->GetValue();
        bPressEnabled = cType == cBoundCondType[2];
        bTempEnabled = cType == cBoundCondType[0] || cType == cBoundCondType[2];
      }
    }

    pData = (DWORD_PTR)&(pBC->fPress);
    pProp = m_wndPropList.FindItemByData(pData);
    if(pProp != NULL)
      pProp->Enable(bEnable && bPressEnabled);

    pData = (DWORD_PTR)&(pBC->fTemp);
    pProp = m_wndPropList.FindItemByData(pData);
    if(pProp != NULL)
      pProp->Enable(bEnable && bTempEnabled);

    pData = (DWORD_PTR)&(pBC->vRegNames);
    pProp = m_wndPropList.FindItemByData(pData);
    if(pProp != NULL)
      pProp->Enable(bEnable);
  }
}
