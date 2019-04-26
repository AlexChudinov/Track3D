
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
    CMFCPropertyGridProperty* pProp = NULL;
    CMFCPropertyGridProperty* pGenGroup = new CMFCPropertyGridProperty(_T("General Field Properties"));

    CCheckBoxButton* pCheckBox = new CCheckBoxButton(this, _T("Enable Field"), (_variant_t)pData->get_enable_field(), _T("Click this check-box to enable / disable this field."), pData->get_enable_field_ptr());
    pGenGroup->AddSubItem(pCheckBox);

    COleVariant var(CElectricFieldData::get_field_type_name(pData->get_type()));
    CGeneralResponseProperty* pFieldType = new CGeneralResponseProperty(this, _T("Field Type (DC, RF or Mirror)"), var, _T("Specify the electric field type. It can be either Direct Current or Radio-Frequency or Coulomb Mirror field."), pData->get_type_ptr());
    for(int k = CElectricFieldData::typeFieldDC; k < CElectricFieldData::typeCount; k++)
      pFieldType->AddOption(CElectricFieldData::get_field_type_name(k));

    pGenGroup->AddSubItem(pFieldType);

    pProp = new CMFCPropertyGridProperty(_T("Voltage Scale, V"), COleVariant(pData->get_scale()), _T("Set here the desirable value of voltage at the selected electrodes."), pData->get_scale_ptr());
    pGenGroup->AddSubItem(pProp);

    if(pData->get_type() == CElectricFieldData::typeFieldRF)
    {
      pProp = new CMFCPropertyGridProperty(_T("Frequency, kHz"), COleVariant(0.001 * pData->get_freq()), _T("Set the frequency of the radio-frequency field."), pData->get_freq_ptr());
      pGenGroup->AddSubItem(pProp);
    }

// Enable / disable field visualization:
    CContourRangeCheckBox* pVisCheckBox = new CContourRangeCheckBox(this, _T("Enable Visualization"), (_variant_t)pData->get_enable_vis(), _T("Enable / Disable this field as a component in the Electric Potential contour."), pData->get_enable_vis_ptr());
    pGenGroup->AddSubItem(pVisCheckBox);

// Field calculation method:
    COleVariant varCalc(CElectricFieldData::get_calc_method_name(pData->get_calc_method()));
    CMFCPropertyGridProperty* pCalcMethod = new CMFCPropertyGridProperty(_T("Field Calc. Method"), varCalc, _T("Specify the field calculation method."), pData->get_calc_method_ptr());
    for(int j = CElectricFieldData::cmLaplacian3; j < CElectricFieldData::cmCount; j++)
      pCalcMethod->AddOption(CElectricFieldData::get_calc_method_name(j));

    pGenGroup->AddSubItem(pCalcMethod);

    CString cDummy(_T(" "));
// Calculate field button:
    CCalcFieldButton* pCalcField = new CCalcFieldButton(this, _T("Calculate Field"), cDummy, _T("Click this button to start field calculation."), (DWORD_PTR)pData);
    bool bEnableCalc = pData->get_type() != CElectricFieldData::typeMirror;
    pCalcField->Enable(bEnableCalc);
    pGenGroup->AddSubItem(pCalcField);

    m_wndPropList.AddProperty(pGenGroup);

// Solver parameters (iterations count, tolerance, multithreading):
    CMFCPropertyGridProperty* pSolverParamGroup = new CMFCPropertyGridProperty(_T("Solver Parameters"));

    pProp = new CMFCPropertyGridProperty(_T("Max. Iterations Count"), COleVariant((long)pData->get_iter_count()), _T("Set the maximal count of iterations for the Laplacian Solver."), pData->get_iter_count_ptr());
    pSolverParamGroup->AddSubItem(pProp);

    pProp = new CMFCPropertyGridProperty(_T("Tolerance"), COleVariant(pData->get_tol()), _T("Specify tolerance for the Laplacian Solver. Note: if the maxinal relative error becomes less than the tolerance, the solution terminates."), pData->get_tol_ptr());
    pSolverParamGroup->AddSubItem(pProp);

// Enable / disable multithreading during field calculation.
    CCheckBoxButton* pEnableMultiThread = new CCheckBoxButton(this, _T("Enable Multithreading"), (_variant_t)pData->get_enable_multithread(), _T("Click this check-box to enable / disable multithreading during calculation of this field."), pData->get_enable_multithread_ptr());
    pSolverParamGroup->AddSubItem(pEnableMultiThread);

    m_wndPropList.AddProperty(pSolverParamGroup);

// An attempt to simulate analytic RF field in the flatapole. Alpha version.
    if(pData->get_type() == CElectricFieldData::typeFieldRF)
    {
      CMFCPropertyGridProperty* pAnalytGroup = new CMFCPropertyGridProperty(_T("Analytic RF Field"));
      CCheckBoxButton* pEnableAnalytCheck = new CCheckBoxButton(this, _T("Enable Analyt RF Field"), (_variant_t)pData->get_analyt_field(), _T("Enable / disable the analytic expression of RF field (in flatapole). This is the alpha version."), pData->get_analyt_field_ptr());
      pAnalytGroup->AddSubItem(pEnableAnalytCheck);

      pProp = new CMFCPropertyGridProperty(_T("Low X Limit, mm"), COleVariant(10 * pData->get_low_analyt_lim()), _T("Currently the analytic model of the RF field is applied when X is greater than this value."), pData->get_low_analyt_lim_ptr());
      pAnalytGroup->AddSubItem(pProp);

      pProp = new CMFCPropertyGridProperty(_T("Inscribed Radius, mm"), COleVariant(10 * pData->get_inscr_radius()), _T("Currently the analytic model of the RF field simulates an ideal quadrupole field. This is the radius of an inscribed circle between the electrodes."), pData->get_inscr_radius_ptr());
      pAnalytGroup->AddSubItem(pProp);

      m_wndPropList.AddProperty(pAnalytGroup);
    }

    CMFCPropertyGridProperty* pButtonsGroup = new CMFCPropertyGridProperty(_T("Actions"));

// Remove field button:
    CRemoveFieldButton* pRemField = new CRemoveFieldButton(this, _T("Remove this Field"), cDummy, _T("Click this button to remove the field from the scene."), (DWORD_PTR)pData);
    pButtonsGroup->AddSubItem(pRemField);

// Add boundary conditions button:
    CAddFieldBoundCondButton* pAddBCButon = new CAddFieldBoundCondButton(this, _T("Add Boundary Conditions"), cDummy, _T("Click to add boundary conditions for this field."), (DWORD_PTR)pData);
    pButtonsGroup->AddSubItem(pAddBCButon);

    m_wndPropList.AddProperty(pButtonsGroup);

// Boundary conditions themselves:
    size_t nBoundCondCount = pData->get_bc_count();
    for(size_t j = 0; j < nBoundCondCount; j++)
    {
      CPotentialBoundCond* pBC = pData->get_bc(j);
      CMFCPropertyGridProperty* pBoundCondGroup = new CMFCPropertyGridProperty(pBC->sName.c_str());

      COleVariant var1(CPotentialBoundCond::get_bc_type_name(pBC->nType));
      CGeneralResponseProperty* pBCType = new CGeneralResponseProperty(this, _T("Boundary Conditions Type"), var1, _T("Set the boundary conditions type (fixed value or zero normal derivative)."), (DWORD_PTR)&(pBC->nType));
      pBCType->AddOption(CPotentialBoundCond::get_bc_type_name(BoundaryMesh::FIXED_VAL));
      pBCType->AddOption(CPotentialBoundCond::get_bc_type_name(BoundaryMesh::ZERO_GRAD));
      pBoundCondGroup->AddSubItem(pBCType);

      if(pBC->nType == BoundaryMesh::FIXED_VAL)
      {
        COleVariant var2(CPotentialBoundCond::get_fixed_value_name(pBC->nFixedValType));
        CGeneralResponseProperty* pBCValue = new CGeneralResponseProperty(this, _T("Boundary Value"), var2, _T("Set the boundary conditions value."), (DWORD_PTR)&(pBC->nFixedValType));
        for(int k = CPotentialBoundCond::fvPlusUnity; k < CPotentialBoundCond::fvCount; k++)
          pBCValue->AddOption(CPotentialBoundCond::get_fixed_value_name(k));

        pBoundCondGroup->AddSubItem(pBCValue);
      }

      if(pBC->nFixedValType == CPotentialBoundCond::fvLinearStepsX || pBC->nFixedValType == CPotentialBoundCond::fvLinearStepsY)
      {
        CMFCPropertyGridProperty* pStepWiseGroup = new CMFCPropertyGridProperty(_T("Linear Step-Wise Potential"));
        CMFCPropertyGridProperty* pLinGradGroup = new CMFCPropertyGridProperty(_T("Linear Gradient"));

        CString cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitStartX);
        CString cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitStartX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->fStartX), cHint, (DWORD_PTR)&(pBC->fStartX));
        pLinGradGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitEndX);
        cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitEndX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->fEndX), cHint, (DWORD_PTR)&(pBC->fEndX));
        pLinGradGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitEndPhi);
        cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitEndPhi);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(pBC->fEndPhi), cHint, (DWORD_PTR)&(pBC->fEndPhi));
        pLinGradGroup->AddSubItem(pProp);

        pStepWiseGroup->AddSubItem(pLinGradGroup);
        CMFCPropertyGridProperty* pElectrGroup = new CMFCPropertyGridProperty(_T("Electrodes"));

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitCenterFirst);
        cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitCenterFirst);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->fCenterFirstElectr), cHint, (DWORD_PTR)&(pBC->fCenterFirstElectr));
        pElectrGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitStepX);
        cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitStepX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->fStepX), cHint, (DWORD_PTR)&(pBC->fStepX));
        pElectrGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitStepsCount);
        cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitStepsCount);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant((long)(pBC->nStepsCount)), cHint, (DWORD_PTR)&(pBC->nStepsCount));
        pElectrGroup->AddSubItem(pProp);

        pStepWiseGroup->AddSubItem(pElectrGroup);
        pBoundCondGroup->AddSubItem(pStepWiseGroup);
      }

      if(pBC->nFixedValType == CPotentialBoundCond::fvQuadricStepsX || pBC->nFixedValType == CPotentialBoundCond::fvQuadricStepsY)
      {
        CMFCPropertyGridProperty* pStepWiseGroup = new CMFCPropertyGridProperty(_T("Quadric Step-Wise Potential"));

        CString cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitStartX);
        CString cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitStartX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->fStartX), cHint, (DWORD_PTR)&(pBC->fStartX));
        pStepWiseGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitStepX);
        cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitStepX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->fStepX), cHint, (DWORD_PTR)&(pBC->fStepX));
        pStepWiseGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitEndX);
        cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitEndX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->fEndX), cHint, (DWORD_PTR)&(pBC->fEndX));
        pStepWiseGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitStartPhi);
        cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitStartPhi);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(pBC->fStartPhi), cHint, (DWORD_PTR)&(pBC->fStartPhi));
        pStepWiseGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitFirstStepPhi);
        cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitFirstStepPhi);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(pBC->fFirstStepPhi), cHint, (DWORD_PTR)&(pBC->fFirstStepPhi));
        pStepWiseGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->nFixedValType, CPotentialBoundCond::uitLastStepPhi);
        cHint = CPotentialBoundCond::get_hint(pBC->nFixedValType, CPotentialBoundCond::uitLastStepPhi);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(pBC->fLastStepPhi), cHint, (DWORD_PTR)&(pBC->fLastStepPhi));
        pStepWiseGroup->AddSubItem(pProp);

        double fEndPhi = pData->quadric_step_potential(pBC, Vector3D(pBC->fEndX, pBC->fEndX, 0));
        CGeneralResponseProperty* pReadOnly = new CGeneralResponseProperty(this, _T("End Potential, V"), COleVariant(fEndPhi), _T(""), NULL);
        pReadOnly->AllowEdit(FALSE);
        pReadOnly->Enable(FALSE);
        pStepWiseGroup->AddSubItem(pReadOnly);

        pBoundCondGroup->AddSubItem(pStepWiseGroup);
      }

// Select boundary regions group:
      CMFCPropertyGridProperty* pBoundRegsGroup = new CMFCPropertyGridProperty(_T("Boundary Regions"));

      CString cRegNames = EvaporatingParticle::CObject::compile_string(pBC->vRegNames);
      CSelectRegionButton* pSelRegButton = new CSelectRegionButton(this, _T("Select Regions"), cRegNames, _T("Click to select 2D regions in the main view window for boundary conditions."), (DWORD_PTR)&(pBC->vRegNames));
      pBoundRegsGroup->AddSubItem(pSelRegButton);
      
      CHideShowRegsCheckBox* pHideShowBtn = new CHideShowRegsCheckBox(this, _T("Visibility"), (_variant_t)pBC->bVisible, _T("Change the visibility status of the selected regions"), (DWORD_PTR)&(pBC->bVisible));
      pBoundRegsGroup->AddSubItem(pHideShowBtn);
      
      pBoundCondGroup->AddSubItem(pBoundRegsGroup);

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
// Field type:
    CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pData->get_type_ptr());
    if(pProp != NULL)
    {
      CString cFieldType = (CString)pProp->GetValue();
      for(int j = CElectricFieldData::typeFieldDC; j < CElectricFieldData::typeCount; j++)
      {
        if(cFieldType == CElectricFieldData::get_field_type_name(j))
        {
          pData->set_type(j);
          break;
        }
      }
    }

// Field calculation method:
    pProp = m_wndPropList.FindItemByData(pData->get_calc_method_ptr());
    if(pProp != NULL)
    {
      CString cCalcMethod = (CString)pProp->GetValue();
      for(int j = CElectricFieldData::cmLaplacian3; j < CElectricFieldData::cmCount; j++)
      {
        if(cCalcMethod == CElectricFieldData::get_calc_method_name(j))
        {
          pData->set_calc_method(j);
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

    pProp = m_wndPropList.FindItemByData(pData->get_tol_ptr());
    if(pProp != NULL)
      pData->set_tol(pProp->GetValue().dblVal);

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
        bool bRecalc = false;
        CString cBCType = (CString)pProp->GetValue();
        if(cBCType == CPotentialBoundCond::get_bc_type_name(BoundaryMesh::FIXED_VAL))
          bRecalc = pBC->set_bc_type(BoundaryMesh::FIXED_VAL);
        if(cBCType == CPotentialBoundCond::get_bc_type_name(BoundaryMesh::ZERO_GRAD))
          bRecalc = pBC->set_bc_type(BoundaryMesh::ZERO_GRAD);

        if(bRecalc)
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->nFixedValType));
      if(pProp != NULL)
      {
        CString cBCValue = (CString)pProp->GetValue();
        for(int i = CPotentialBoundCond::fvPlusUnity; i < CPotentialBoundCond::fvCount; i++)
        {
          if(cBCValue == CPotentialBoundCond::get_fixed_value_name(i))
          {
            if(pBC->set_fixed_val_type(i))
              pData->invalidate();
            break;
          }
        }
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->fStartX));
      if(pProp != NULL)
      {
        if(pBC->set_start_coord(0.1 * pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->fStepX));
      if(pProp != NULL)
      {
        if(pBC->set_step_coord(0.1 * pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->fEndX));
      if(pProp != NULL)
      {
        if(pBC->set_end_coord(0.1 * pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->fEndPhi));
      if(pProp != NULL)
      {
        if(pBC->set_end_phi(pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->fStartPhi));
      if(pProp != NULL)
      {
        if(pBC->set_start_phi(pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->fFirstStepPhi));
      if(pProp != NULL)
      {
        if(pBC->set_first_dphi(pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->fLastStepPhi));
      if(pProp != NULL)
      {
        if(pBC->set_last_dphi(pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->fCenterFirstElectr));
      if(pProp != NULL)
      {
        if(pBC->set_center_first_electr(0.1 * pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData((DWORD_PTR)&(pBC->nStepsCount));
      if(pProp != NULL)
      {
        if(pBC->set_steps_count(pProp->GetValue().lVal))
          pData->invalidate();
      }
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
    bool bMirrorField = false;
    CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pData->get_type_ptr());
    if(pProp != NULL)
    {
      CString cFieldType = (CString)pProp->GetValue();
      bFieldRF = cFieldType == pData->get_field_type_name(CElectricFieldData::typeFieldRF);
      bMirrorField = cFieldType == pData->get_field_type_name(CElectricFieldData::typeMirror);
    }

    pProp = m_wndPropList.FindItemByData(pData->get_scale_ptr());
    if(pProp != NULL)
      pProp->Enable(!bMirrorField);

    pProp = m_wndPropList.FindItemByData(pData->get_freq_ptr());
    if(pProp != NULL)
      pProp->Enable(bFieldRF);

    pProp = m_wndPropList.FindItemByData(pData->get_analyt_field_ptr());
    if(pProp != NULL)
      pProp->Enable(bFieldRF);

    bool bAnalyt = pData->get_analyt_field();
    pProp = m_wndPropList.FindItemByData(pData->get_low_analyt_lim_ptr());
    if(pProp != NULL)
      pProp->Enable(bFieldRF && bAnalyt);

    pProp = m_wndPropList.FindItemByData(pData->get_inscr_radius_ptr());
    if(pProp != NULL)
      pProp->Enable(bFieldRF && bAnalyt);

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
