
#include "stdafx.h"

#include "../field_solver/MeshData.h"
#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "ResponseProperty.h"
#include "Button.h"

using namespace EvaporatingParticle;

//---------------------------------------------------------------------------------------
// Built-in electric fields computation.
//---------------------------------------------------------------------------------------
void CPropertiesWnd::add_field_ctrls()
{
  CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
  size_t nFieldsCount = pFields->size();
  int nCurrFieldId = pFields->get_curr_field_index();
  CString sCurrFieldName = nCurrFieldId >= 0 ? pFields->at(nCurrFieldId)->get_field_name() : CString("None");
  
  CElectricFieldResponder* pFieldSelector = new CElectricFieldResponder(this, _T("Selected Electrtic Field"), COleVariant(sCurrFieldName), _T("Select electric field to view and modify its parameters."), (DWORD_PTR)pFields);
  for(size_t i = 0; i < nFieldsCount; i++)
  {
    CElectricFieldData* pData = pFields->at(i);
    pFieldSelector->AddOption(pData->get_field_name());
  }

  pFieldSelector->AllowEdit(TRUE);
  m_wndPropList.AddProperty(pFieldSelector);

  CElectricFieldData* pData = nCurrFieldId >= 0 ? pFields->at(nCurrFieldId) : NULL;
  if(pData != NULL)
  {
    CCheckBoxButton* pCheckBox = new CCheckBoxButton(this, _T("Enable Field"), (_variant_t)pData->get_enable_field(), _T("Click this check-box to enable / disable this field."), pData->get_enable_field_ptr());
    m_wndPropList.AddProperty(pCheckBox);

    COleVariant var(CElectricFieldData::get_field_type_name(pData->get_type()));
    CMFCPropertyGridProperty* pFieldType = new CMFCPropertyGridProperty(_T("Field Type (DC or RF)"), var, _T("Specify the electric field type. It can be either Direct Current or Radio-Frequency field."), pData->get_type_ptr());
    for(int k = CElectricFieldData::typeFieldDC; k < CElectricFieldData::typeCount; k++)
      pFieldType->AddOption(pData->get_field_type_name(k));

    m_wndPropList.AddProperty(pFieldType);

    CMFCPropertyGridProperty* pProp = new CMFCPropertyGridProperty(_T("Voltage Scale, V"), COleVariant(pData->get_scale()), _T("Set here the desirable value of voltage at the selected electrodes."), pData->get_scale_ptr());
    m_wndPropList.AddProperty(pProp);

    pProp = new CMFCPropertyGridProperty(_T("Frequency, kHz"), COleVariant(0.001 * pData->get_freq()), _T("Set the frequency of the radio-frequency field."), pData->get_freq_ptr());
    m_wndPropList.AddProperty(pProp);

    CString cDummy(_T(" "));
// Calculate field button:
    CCalcFieldButton* pCalcField = new CCalcFieldButton(this, _T("Calculate Field"), cDummy, _T("Click this button to start field calculation."), (DWORD_PTR)pData);
    m_wndPropList.AddProperty(pCalcField);

    pProp = new CMFCPropertyGridProperty(_T("Iterations Count"), COleVariant((long)pData->get_iter_count()), _T("Set the count of iterations for the Laplacian Soilver."), pData->get_iter_count_ptr());
    m_wndPropList.AddProperty(pProp);

// An attempt to simulate analytic RF field in the flatapole. Alpha version.
    CMFCPropertyGridProperty* pAnalytGroup = new CMFCPropertyGridProperty(_T("Analytic RF Field"));
    CCheckBoxButton* pEnableAnalytCheck = new CCheckBoxButton(this, _T("Enable Analyt RF Field"), (_variant_t)pData->get_analyt_field(), _T("Enable / disable the analytic expression of RF field (in flatapole). This is the alpha version."), pData->get_analyt_field_ptr());
    pAnalytGroup->AddSubItem(pEnableAnalytCheck);

    pProp = new CMFCPropertyGridProperty(_T("Low X Limit, mm"), COleVariant(10 * pData->get_low_analyt_lim()), _T("Currently the analytic model of the RF field is applied when X is greater than this value."), pData->get_low_analyt_lim_ptr());
    pAnalytGroup->AddSubItem(pProp);

    pProp = new CMFCPropertyGridProperty(_T("Inscribed Radius, mm"), COleVariant(10 * pData->get_inscr_radius()), _T("Currently the analytic model of the RF field simulates an ideal quadrupole field. This is the radius of an inscribed circle between the electrodes."), pData->get_inscr_radius_ptr());
    pAnalytGroup->AddSubItem(pProp);

    m_wndPropList.AddProperty(pAnalytGroup);

// Remove field button:
    CRemoveFieldButton* pRemField = new CRemoveFieldButton(this, _T("Remove this Field"), cDummy, _T("Click this button to remove the field from the scene."), (DWORD_PTR)pData);
    m_wndPropList.AddProperty(pRemField);

// Add boundary conditions button:
    CAddFieldBoundCondButton* pAddBCButon = new CAddFieldBoundCondButton(this, _T("Add Boundary Conditions"), cDummy, _T("Click to add boundary conditions for this field."), (DWORD_PTR)pData);
    m_wndPropList.AddProperty(pAddBCButon);

// Boundary conditions themselves:
    size_t nBoundCondCount = pData->get_bc_count();
    for(size_t j = 0; j < nBoundCondCount; j++)
    {
      CPotentialBoundCond* pBC = pData->get_bc(j);
      CMFCPropertyGridProperty* pBoundCondGroup = new CMFCPropertyGridProperty(pBC->sName.c_str());

      COleVariant var1(CPotentialBoundCond::get_bc_type_name(pBC->nType));
      CMFCPropertyGridProperty* pBCType = new CMFCPropertyGridProperty(_T("Boundary Conditions Type"), var1, _T("Set the boundary conditions type."), (DWORD_PTR)&(pBC->nType));
      pBCType->AddOption(CPotentialBoundCond::get_bc_type_name(BoundaryMesh::FIXED_VAL));
      pBCType->AddOption(CPotentialBoundCond::get_bc_type_name(BoundaryMesh::ZERO_GRAD));
      pBoundCondGroup->AddSubItem(pBCType);

      COleVariant var2(CPotentialBoundCond::get_fixed_value_name(pBC->nFixedValType));
//      CMFCPropertyGridProperty* pBCValue = new CMFCPropertyGridProperty(_T("Boundary Value"), var2, _T("Set the boundary conditions value."), (DWORD_PTR)&(pBC->nFixedValType));
      CGeneralResponseProperty* pBCValue = new CGeneralResponseProperty(this, _T("Boundary Value"), var2, _T("Set the boundary conditions value."), (DWORD_PTR)&(pBC->nFixedValType));
      pBCValue->AddOption(CPotentialBoundCond::get_fixed_value_name(CPotentialBoundCond::fvPlusUnity));
      pBCValue->AddOption(CPotentialBoundCond::get_fixed_value_name(CPotentialBoundCond::fvMinusUnity));
      pBCValue->AddOption(CPotentialBoundCond::get_fixed_value_name(CPotentialBoundCond::fvStepLike));
      pBoundCondGroup->AddSubItem(pBCValue);

      if(pBC->nFixedValType == CPotentialBoundCond::fvStepLike)
      {
        CMFCPropertyGridProperty* pStepWiseGroup = new CMFCPropertyGridProperty(_T("Step-Wise Boundary Conditions"));

        pProp = new CMFCPropertyGridProperty(_T("Start X, mm"), COleVariant(10 * pBC->fStartX), _T("Specify the X-coordinate of beginning of the first step in the step-wise potentials set."), (DWORD_PTR)&(pBC->fStartX));
        pStepWiseGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Step X, mm"), COleVariant(10 * pBC->fStepX), _T("Specify the width of the X-step in the step-wise potentials set."), (DWORD_PTR)&(pBC->fStepX));
        pStepWiseGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("End X, mm"), COleVariant(10 * pBC->fEndX), _T("Specify the X-coordinate of end of the last step in the step-wise potentials set."), (DWORD_PTR)&(pBC->fEndX));
        pStepWiseGroup->AddSubItem(pProp);

        pBoundCondGroup->AddSubItem(pStepWiseGroup);
      }

      CString cRegNames = EvaporatingParticle::CObject::compile_string(pBC->vRegNames);
      CSelectRegionButton* pSelRegButton = new CSelectRegionButton(this, _T("Boundary Regions"), cRegNames, _T("Click to select 2D regions for boundary conditions."), (DWORD_PTR)&(pBC->vRegNames));
      pBoundCondGroup->AddSubItem(pSelRegButton);

// Remove boundary condition button:
      CRemoveFieldBoundCondButton* pRemBCButton = new CRemoveFieldBoundCondButton(this, _T("Remove Boundary Condition"), cDummy, _T("Click to remove this boundary condition."), (DWORD_PTR)pBC);
      pBoundCondGroup->AddSubItem(pRemBCButton);

      m_wndPropList.AddProperty(pBoundCondGroup);
    }
  }
}

void CPropertiesWnd::set_field_data()
{
  CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
  int nCurrFieldId = pFields->get_curr_field_index();
  CElectricFieldData* pData = nCurrFieldId >= 0 ? pFields->at(nCurrFieldId) : NULL;
  if(pData != NULL)
  {
    CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pData->get_type_ptr());
    if(pProp != NULL)
    {
      CString cFieldType = (CString)pProp->GetValue();
      for(int j = CElectricFieldData::typeFieldDC; j < CElectricFieldData::typeCount; j++)
      {
        if(cFieldType == pData->get_field_type_name(j))
        {
          pData->set_type(j);
          break;
        }
      }
    }

    pProp = m_wndPropList.FindItemByData(pData->get_scale_ptr());
    if(pProp != NULL)
      pData->set_scale(pProp->GetValue().dblVal);

    pProp = m_wndPropList.FindItemByData(pData->get_freq_ptr());
    if(pProp != NULL)
      pData->set_freq(1000 * pProp->GetValue().dblVal);

    pProp = m_wndPropList.FindItemByData(pData->get_iter_count_ptr());
    if(pProp != NULL)
      pData->set_iter_count(pProp->GetValue().lVal);

    pProp = m_wndPropList.FindItemByData(pData->get_low_analyt_lim_ptr());
    if(pProp != NULL)
      pData->set_low_analyt_lim(0.1 * pProp->GetValue().dblVal);

    pProp = m_wndPropList.FindItemByData(pData->get_inscr_radius_ptr());
    if(pProp != NULL)
      pData->set_inscr_radius(0.1 * pProp->GetValue().dblVal);

// Boundary conditions:
    size_t nBoundCondCount = pData->get_bc_count();
    for(size_t k = 0; k < nBoundCondCount; k++)
    {
      CPotentialBoundCond* pBC = pData->get_bc(k);

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->nType));
      if(pProp != NULL)
      {
        CString cBCType = (CString)pProp->GetValue();
        if(cBCType == CPotentialBoundCond::get_bc_type_name(BoundaryMesh::FIXED_VAL))
          pBC->nType = BoundaryMesh::FIXED_VAL;
        if(cBCType == CPotentialBoundCond::get_bc_type_name(BoundaryMesh::ZERO_GRAD))
          pBC->nType = BoundaryMesh::ZERO_GRAD;
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->nFixedValType));
      if(pProp != NULL)
      {
        CString cBCValue = (CString)pProp->GetValue();
        if(cBCValue == CPotentialBoundCond::get_fixed_value_name(CPotentialBoundCond::fvPlusUnity))
          pBC->nFixedValType = CPotentialBoundCond::fvPlusUnity;
        if(cBCValue == CPotentialBoundCond::get_fixed_value_name(CPotentialBoundCond::fvMinusUnity))
          pBC->nFixedValType = CPotentialBoundCond::fvMinusUnity;
        if(cBCValue == CPotentialBoundCond::get_fixed_value_name(CPotentialBoundCond::fvStepLike))
          pBC->nFixedValType = CPotentialBoundCond::fvStepLike;
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->fStartX));
      if(pProp != NULL)
        pBC->fStartX = 0.1 * pProp->GetValue().dblVal;

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->fStepX));
      if(pProp != NULL)
        pBC->fStepX = 0.1 * pProp->GetValue().dblVal;

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->fEndX));
      if(pProp != NULL)
        pBC->fEndX = 0.1 * pProp->GetValue().dblVal;
    }
  }
}

void CPropertiesWnd::update_field_ctrls()
{
  CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
  int nCurrFieldId = pFields->get_curr_field_index();
  CElectricFieldData* pData = nCurrFieldId >= 0 ? pFields->at(nCurrFieldId) : NULL;
  if(pData != NULL)
  {
    bool bFieldRF = false;
    CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pData->get_type_ptr());
    if(pProp != NULL)
    {
      CString cFieldType = (CString)pProp->GetValue();
      bFieldRF = cFieldType == pData->get_field_type_name(CElectricFieldData::typeFieldRF);
    }

    pProp = m_wndPropList.FindItemByData(pData->get_freq_ptr());
    if(pProp != NULL)
      pProp->Enable(bFieldRF);

    pProp = m_wndPropList.FindItemByData(pData->get_analyt_field_ptr());
    if(pProp != NULL)
      pProp->Enable(bFieldRF);

    pProp = m_wndPropList.FindItemByData(pData->get_low_analyt_lim_ptr());
    if(pProp != NULL)
      pProp->Enable(bFieldRF);

    pProp = m_wndPropList.FindItemByData(pData->get_inscr_radius_ptr());
    if(pProp != NULL)
      pProp->Enable(bFieldRF);

// Boundary conditions:
    size_t nBoundCondCount = pData->get_bc_count();
    for(size_t k = 0; k < nBoundCondCount; k++)
    {
      CPotentialBoundCond* pBC = pData->get_bc(k);

      bool bZeroGrad = false;
      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->nType));
      if(pProp != NULL)
      {
        CString cBCType = (CString)pProp->GetValue();
        if(cBCType == CPotentialBoundCond::get_bc_type_name(BoundaryMesh::ZERO_GRAD))
          bZeroGrad = true;
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->nFixedValType));
      if(pProp != NULL)
        pProp->Enable(!bZeroGrad);
    }
  }
}
