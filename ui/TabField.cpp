
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
  CMFCPropertyGridProperty* pMainFieldProp = new CMFCPropertyGridProperty(_T("Electric Fields"));
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
  pMainFieldProp->AddSubItem(pFieldSelector);
  m_wndPropList.AddProperty(pMainFieldProp);

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
      CMFCPropertyGridProperty* pBoundCondGroup = new CMFCPropertyGridProperty(pBC->get_name());

      CString sType(CPotentialBoundCond::get_bc_type_name(pBC->get_bc_type()));
      CGeneralResponseProperty* pBCType = new CGeneralResponseProperty(this, _T("Boundary Conditions Type"), sType, _T("Set the boundary conditions type (fixed value or zero normal derivative)."), pBC->get_bc_type_ptr());
      pBCType->AddOption(CPotentialBoundCond::get_bc_type_name(BoundaryMesh::FIXED_VAL));
      pBCType->AddOption(CPotentialBoundCond::get_bc_type_name(BoundaryMesh::ZERO_GRAD));
      pBoundCondGroup->AddSubItem(pBCType);

      if(pBC->get_bc_type() == BoundaryMesh::FIXED_VAL)
      {
        CString sValType(CPotentialBoundCond::get_fixed_value_name(pBC->get_fixed_val_type()));
        CGeneralResponseProperty* pBCValue = new CGeneralResponseProperty(this, _T("Boundary Value"), sValType, _T("Set the boundary conditions value."), pBC->get_fixed_val_type_ptr());
        for(int k = CPotentialBoundCond::fvPlusUnity; k < CPotentialBoundCond::fvCount; k++)
          pBCValue->AddOption(CPotentialBoundCond::get_fixed_value_name(k));

        pBoundCondGroup->AddSubItem(pBCValue);
      }

      if(pBC->get_fixed_val_type() == CPotentialBoundCond::fvLinearStepsX || pBC->get_fixed_val_type() == CPotentialBoundCond::fvLinearStepsY)
      {
        CMFCPropertyGridProperty* pStepWiseGroup = new CMFCPropertyGridProperty(_T("Linear Step-Wise Potential"));
        CMFCPropertyGridProperty* pLinGradGroup = new CMFCPropertyGridProperty(_T("Linear Gradient"));

        CString cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStartX);
        CString cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStartX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->get_start_coord()), cHint, pBC->get_start_coord_ptr());
        pLinGradGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitEndX);
        cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitEndX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->get_end_coord()), cHint, pBC->get_end_coord_ptr());
        pLinGradGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitEndPhi);
        cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitEndPhi);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(pBC->get_end_phi()), cHint, pBC->get_end_phi_ptr());
        pLinGradGroup->AddSubItem(pProp);

        pStepWiseGroup->AddSubItem(pLinGradGroup);
        CMFCPropertyGridProperty* pElectrGroup = new CMFCPropertyGridProperty(_T("Electrodes"));

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitCenterFirst);
        cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitCenterFirst);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->get_center_first_electr()), cHint, pBC->get_center_first_electr_ptr());
        pElectrGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStepX);
        cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStepX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->get_step_coord()), cHint, pBC->get_step_coord_ptr());
        pElectrGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStepsCount);
        cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStepsCount);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant((long)(pBC->get_steps_count())), cHint, pBC->get_steps_count_ptr());
        pElectrGroup->AddSubItem(pProp);

        pStepWiseGroup->AddSubItem(pElectrGroup);
        pBoundCondGroup->AddSubItem(pStepWiseGroup);
      }

      if(pBC->get_fixed_val_type() == CPotentialBoundCond::fvQuadricStepsX || pBC->get_fixed_val_type() == CPotentialBoundCond::fvQuadricStepsY)
      {
        CMFCPropertyGridProperty* pStepWiseGroup = new CMFCPropertyGridProperty(_T("Quadric Step-Wise Potential"));

        CString cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStartX);
        CString cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStartX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->get_start_coord()), cHint, pBC->get_start_coord_ptr());
        pStepWiseGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStepX);
        cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStepX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->get_step_coord()), cHint, pBC->get_step_coord_ptr());
        pStepWiseGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitEndX);
        cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitEndX);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(10 * pBC->get_end_coord()), cHint, pBC->get_end_coord_ptr());
        pStepWiseGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStartPhi);
        cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitStartPhi);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(pBC->get_start_phi()), cHint, pBC->get_start_phi_ptr());
        pStepWiseGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitFirstStepPhi);
        cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitFirstStepPhi);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(pBC->get_first_dphi()), cHint, pBC->get_first_dphi_ptr());
        pStepWiseGroup->AddSubItem(pProp);

        cCtrlTitle = CPotentialBoundCond::get_control_title(pBC->get_fixed_val_type(), CPotentialBoundCond::uitLastStepPhi);
        cHint = CPotentialBoundCond::get_hint(pBC->get_fixed_val_type(), CPotentialBoundCond::uitLastStepPhi);
        pProp = new CMFCPropertyGridProperty(cCtrlTitle, COleVariant(pBC->get_last_dphi()), cHint, pBC->get_last_dphi_ptr());
        pStepWiseGroup->AddSubItem(pProp);

        double fEndPhi = pData->quadric_step_potential(pBC, Vector3D(pBC->get_end_coord(), pBC->get_end_coord(), 0));
        CGeneralResponseProperty* pReadOnly = new CGeneralResponseProperty(this, _T("End Potential, V"), COleVariant(fEndPhi), _T(""), NULL);
        pReadOnly->AllowEdit(FALSE);
        pReadOnly->Enable(FALSE);
        pStepWiseGroup->AddSubItem(pReadOnly);

        pBoundCondGroup->AddSubItem(pStepWiseGroup);
      }

// Select boundary regions group. NOTE: the GetData() of this group control will return, in fact, the pointer to pBC->vRegNames. This will be used in CHideShowRegsCheckBox.
      CMFCPropertyGridProperty* pBoundRegsGroup = new CMFCPropertyGridProperty(_T("Boundary Regions"), pBC->get_region_names_ptr());

// Manual selection:
      CString cRegNames = EvaporatingParticle::CObject::compile_string(pBC->get_region_names());
      CSelectRegionButton* pSelRegButton = new CSelectRegionButton(this, _T("Select Regions Manually"), cRegNames, _T("Click to select 2D regions in the main view window."), pBC->get_region_names_ptr());
      pBoundRegsGroup->AddSubItem(pSelRegButton);

// Merging with Named Areas:
      CNamedAreasSelResponder* pNamedAreasSelector = new CNamedAreasSelResponder(this, _T("Merge with Named Areas"), pBC->get_last_merged(), _T("Select the existing Named Areas to use these surfaces for boundary condtions setting."), pBC->get_region_names_ptr());
      pNamedAreasSelector->AllowEdit(FALSE);
      pNamedAreasSelector->AddOption(_T("None"));
      EvaporatingParticle::CSelAreasColl* pSelAreasColl = CParticleTrackingApp::Get()->GetSelAreas();
      size_t nSelAreasCount = pSelAreasColl->size();
      for(size_t k = 0; k < nSelAreasCount; k++)
        pNamedAreasSelector->AddOption(pSelAreasColl->at(k)->get_name());

      pBoundRegsGroup->AddSubItem(pNamedAreasSelector);

// Merging options:
      CString sMergeOpt = CSelectedAreas::merge_opt_name(pBC->get_merge_option());
      CMFCPropertyGridProperty* pMergeOptSelector = new CMFCPropertyGridProperty(_T("Merge Options"), sMergeOpt, _T("Select one of three allowed merge options: add, substitute and subtract."), pBC->get_merge_option_ptr());
      for(int l = EvaporatingParticle::CSelectedAreas::optAdd; l < EvaporatingParticle::CSelectedAreas::optCount; l++)
        pMergeOptSelector->AddOption(EvaporatingParticle::CSelectedAreas::merge_opt_name(l));

      pBoundRegsGroup->AddSubItem(pMergeOptSelector);
      

// Hide/Show selected regions:
      CHideShowRegsCheckBox* pHideShowBtn = new CHideShowRegsCheckBox(this, _T("Visibility"), (_variant_t)pBC->get_visibility_flag(), _T("Change the visibility status of the selected regions"), pBC->get_visibility_flag_ptr());
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

      pProp = m_wndPropList.FindItemByData(pBC->get_bc_type_ptr());
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

      pProp = m_wndPropList.FindItemByData(pBC->get_fixed_val_type_ptr());
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

      pProp = m_wndPropList.FindItemByData(pBC->get_start_coord_ptr());
      if(pProp != NULL)
      {
        if(pBC->set_start_coord(0.1 * pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData(pBC->get_step_coord_ptr());
      if(pProp != NULL)
      {
        if(pBC->set_step_coord(0.1 * pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData(pBC->get_end_coord_ptr());
      if(pProp != NULL)
      {
        if(pBC->set_end_coord(0.1 * pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData(pBC->get_end_phi_ptr());
      if(pProp != NULL)
      {
        if(pBC->set_end_phi(pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData(pBC->get_start_phi_ptr());
      if(pProp != NULL)
      {
        if(pBC->set_start_phi(pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData(pBC->get_first_dphi_ptr());
      if(pProp != NULL)
      {
        if(pBC->set_first_dphi(pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData(pBC->get_last_dphi_ptr());
      if(pProp != NULL)
      {
        if(pBC->set_last_dphi(pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData(pBC->get_center_first_electr_ptr());
      if(pProp != NULL)
      {
        if(pBC->set_center_first_electr(0.1 * pProp->GetValue().dblVal))
          pData->invalidate();
      }

      pProp = m_wndPropList.FindItemByData(pBC->get_steps_count_ptr());
      if(pProp != NULL)
      {
        if(pBC->set_steps_count(pProp->GetValue().lVal))
          pData->invalidate();
      }

// Merge options:
      pProp = m_wndPropList.FindItemByData(pBC->get_merge_option_ptr());
      if(pProp != NULL)
      {
        CString sOptSel = pProp->GetValue();
        for(int i = CSelectedAreas::optAdd; i < CSelectedAreas::optCount; i++)
        {
          if(sOptSel == CSelectedAreas::merge_opt_name(i))
          {
            pBC->set_merge_option(i);
            break;
          }
        }
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
      pProp = m_wndPropList.FindItemByData(pBC->get_bc_type_ptr());
      if(pProp != NULL)
      {
        CString cBCType = (CString)pProp->GetValue();
        if(cBCType == CPotentialBoundCond::get_bc_type_name(BoundaryMesh::ZERO_GRAD))
          bZeroGrad = true;
      }

      pProp = m_wndPropList.FindItemByData(pBC->get_fixed_val_type_ptr());
      if(pProp != NULL)
        pProp->Enable(!bZeroGrad);
    }
  }
}
