
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "ResponseProperty.h"
#include "Button.h"

static const size_t scnStepCount = 3;
static const char* sccImportSteps[scnStepCount] = { "First (Spade-Work)", "Second (Run-Time)", "Restore ANSYS Data" };
//---------------------------------------------------------------------------------------
//  Import OpenFOAM and Calculators
//---------------------------------------------------------------------------------------
void CPropertiesWnd::add_import_ctrls()
{
// Import OpenFOAM:
  EvaporatingParticle::CImportOpenFOAM& ImpObj = CParticleTrackingApp::Get()->GetTracker()->get_importer();

  CMFCPropertyGridProperty* pImportGroup = new CMFCPropertyGridProperty(_T("Import OpenFOAM Data"));

// Import consists of two steps:
  COleVariant var(_T(sccImportSteps[ImpObj.get_step()]));
  CMFCPropertyGridProperty* pStepProp = new CMFCPropertyGridProperty(_T("Step of Import"), var, _T("Specify the step of the import: 'None', 'First (Spade-Work)' or 'Second (Run-Time)'."), ImpObj.get_step_ptr());
  for(size_t i = 0; i < scnStepCount; i++)
    pStepProp->AddOption(_T(sccImportSteps[i]));

  pStepProp->AllowEdit(FALSE);
  pImportGroup->AddSubItem(pStepProp);

  CMFCPropertyGridProperty* pShiftProp = new CMFCPropertyGridProperty(_T("X-Coordinate Shift, mm"), COleVariant(10 * ImpObj.get_x_shift()), _T("Defines the shift of the exported part of the mesh with respect to the whole mesh."), ImpObj.get_x_shift_ptr());
  pImportGroup->AddSubItem(pShiftProp);

// Data filename:
  CMFCPropertyGridProperty* pDataGroup = new CMFCPropertyGridProperty(_T("Import Data File"));
  static TCHAR BASED_CODE szFilter[] = _T("Data Files(*.csv)|*.csv|Geometry Files(*.geom)|*.geom|All Files(*.*)|*.*||");
  CMFCPropertyGridFileProperty* pImportFileProp = new CMFCPropertyGridFileProperty(_T("Data File"), TRUE, ImpObj.get_filename(), _T("geom"),
    OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, _T("Specify location of the data file for import. Use '*.csv' files for the first step of the import and '*.geom' files for the second step."), ImpObj.get_filename_ptr());

  pDataGroup->AddSubItem(pImportFileProp);
  pImportGroup->AddSubItem(pDataGroup);

  m_wndPropList.AddProperty(pImportGroup);
}

void CPropertiesWnd::set_import_data()
{
  EvaporatingParticle::CImportOpenFOAM& ImpObj = CParticleTrackingApp::Get()->GetTracker()->get_importer();

  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(ImpObj.get_step_ptr());
  if(pProp != NULL)
  {
    CString cStep = (CString)pProp->GetValue();
    for(size_t i = 0; i < scnStepCount; i++)
    {
      if(cStep == sccImportSteps[i])
      {
        ImpObj.set_step(i);
        break;
      }
    }
  }

  pProp = m_wndPropList.FindItemByData(ImpObj.get_x_shift_ptr());
  if(pProp != NULL)
    ImpObj.set_x_shift(0.1 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(ImpObj.get_filename_ptr());
  if(pProp != NULL)
  {
    CString cFile = (CString)pProp->GetValue();
    std::string str = std::string(CT2CA(cFile));
    ImpObj.set_filename(str.c_str());
  }
}

void CPropertiesWnd::update_import_ctrls()
{
  bool bEnableFileName = false, bEnableShift = false;;
  bool bReady = !m_bBusy && CParticleTrackingApp::Get()->GetTracker()->is_ready();

  EvaporatingParticle::CImportOpenFOAM& ImpObj = CParticleTrackingApp::Get()->GetTracker()->get_importer();

  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(ImpObj.get_step_ptr());
  if(pProp != NULL)
  {
    CString cStep = (CString)pProp->GetValue();
    bEnableShift = cStep == sccImportSteps[0];
    bEnableFileName = cStep != sccImportSteps[2];
    pProp->Enable(bReady);
  }

  pProp = m_wndPropList.FindItemByData(ImpObj.get_x_shift_ptr());
  if(pProp != NULL)
    pProp->Enable(bReady && bEnableShift);

  pProp = m_wndPropList.FindItemByData(ImpObj.get_filename_ptr());
  if(pProp != NULL)
    pProp->Enable(bReady && bEnableFileName);
}
