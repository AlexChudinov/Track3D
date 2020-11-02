
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "ResponseProperty.h"
#include "Button.h"

using namespace EvaporatingParticle;

void CPropertiesWnd::add_ptb_ctrls()
{
  CFieldPtbCollection& coll = CParticleTrackingApp::Get()->GetTracker()->get_field_ptb();
  size_t nPtbCount = coll.size();
  for(size_t i = 0; i < nPtbCount; i++)
  {
    CFieldPerturbation* pPtb = coll.at(i);
    int nType = pPtb->type();
    switch(nType)
    {
      case CFieldPerturbation::ptbRing:
      {
        CChargedRingPerturbation* pChargedRing = (CChargedRingPerturbation*)pPtb;
        CMFCPropertyGridProperty* pChargedRingGroup = new CMFCPropertyGridProperty(pChargedRing->name());
// Enable:
        CCheckBoxButton* pCheckBox = new CCheckBoxButton(this, _T("Enable"), (_variant_t)pChargedRing->get_enable(), _T("Turns ON/OFF the field perturbation."), pChargedRing->get_enable_ptr());
        pChargedRingGroup->AddSubItem(pCheckBox);
// Ring center position:
        CMFCPropertyGridProperty* pPosGroup = new CMFCPropertyGridProperty(_T("Position, mm"), pChargedRing->get_ring_pos_ptr());
        CMFCPropertyGridProperty* pProp = new CMFCPropertyGridProperty(_T("X"), COleVariant(10 * pChargedRing->get_ring_pos().x), _T("X coordinate of the charged ring center."));
        pPosGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Y"), COleVariant(10 * pChargedRing->get_ring_pos().y), _T("Y coordinate of the charged ring center."));
        pPosGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Z"), COleVariant(10 * pChargedRing->get_ring_pos().z), _T("Z coordinate of the charged ring center."));
        pPosGroup->AddSubItem(pProp);
        pChargedRingGroup->AddSubItem(pPosGroup);
// Ring radius and charge:
        CMFCPropertyGridProperty* pChargeRadiusGroup = new CMFCPropertyGridProperty(_T("Charge and Radius"));
        pProp = new CMFCPropertyGridProperty(_T("Charge, 10^6 elem. charges"), COleVariant(1e-6 * pChargedRing->get_ring_charge() / Const_Charge_CGS), _T("Set the charge at the ring in elementary charges."), pChargedRing->get_ring_charge_ptr());
        pChargeRadiusGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Ring Radius, mm"), COleVariant(10 * pChargedRing->get_ring_radius()), _T("Set the radius of the ring."), pChargedRing->get_ring_radius_ptr());
        pChargeRadiusGroup->AddSubItem(pProp);
        pChargedRingGroup->AddSubItem(pChargeRadiusGroup);
// Remove the perturbation:
        COleVariant var(_T(""));
        CRemovePerturbationButton* pRemBtn = new CRemovePerturbationButton(this, _T("Remove Perturbation"), var, _T("Click to delete this perturbation."), (DWORD_PTR)pPtb);
        pChargedRingGroup->AddSubItem(pRemBtn);

        m_wndPropList.AddProperty(pChargedRingGroup);
        break;
      }
      case CFieldPerturbation::ptbStackOfRings:
      {
        CStackRingPerturbation* pStackPtb = (CStackRingPerturbation*)pPtb;
        CMFCPropertyGridProperty* pStackRingGroup = new CMFCPropertyGridProperty(pStackPtb->name());
// Enable:
        CCheckBoxButton* pCheckBox = new CCheckBoxButton(this, _T("Enable"), (_variant_t)pStackPtb->get_enable(), _T("Turns ON/OFF the field perturbation."), pStackPtb->get_enable_ptr());
        pStackRingGroup->AddSubItem(pCheckBox);
// Count of rings, ring radius and the sum charge:
        CMFCPropertyGridProperty* pProp = new CMFCPropertyGridProperty(_T("Rings Count"), COleVariant((long)pStackPtb->get_rings_count()), _T("Set the count of rings in the stack."), pStackPtb->get_rings_count_ptr());
        pStackRingGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Ring Radius, mm"), COleVariant(10 * pStackPtb->get_ring_radius()), _T("Set the radius of each ring in the stack."), pStackPtb->get_ring_radius_ptr());
        pStackRingGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Full Charge, 10^6 elem. charges"), COleVariant(1e-6 * pStackPtb->get_sum_charge() / Const_Charge_CGS), _T("Set the summary charge in the stack."), pStackPtb->get_sum_charge_ptr());
        pStackRingGroup->AddSubItem(pProp);
// Charge distribution type in the stack:
        COleVariant var(pStackPtb->get_distr_type_name(pStackPtb->get_charge_distr_type()));
        CMFCPropertyGridProperty* pType = new CMFCPropertyGridProperty(_T("Stack Charge Distribution"), var, _T("Specify type of charges distribution in the stack: uniform or regressive."), pStackPtb->get_charge_distr_type_ptr());
        for(int i = 0; i < EvaporatingParticle::CStackRingPerturbation::distrCount; i++)
          pType->AddOption(pStackPtb->get_distr_type_name(i));

        pType->AllowEdit(FALSE);
        pStackRingGroup->AddSubItem(pType);

// First ring position:
        CMFCPropertyGridProperty* pFirstPosGroup = new CMFCPropertyGridProperty(_T("First Ring Position, mm"), pStackPtb->get_stack_beg_pos_ptr());
        pProp = new CMFCPropertyGridProperty(_T("X"), COleVariant(10 * pStackPtb->get_stack_beg_pos().x), _T("X coordinate of the first ring center."));
        pFirstPosGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Y"), COleVariant(10 * pStackPtb->get_stack_beg_pos().y), _T("Y coordinate of the first ring center."));
        pFirstPosGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Z"), COleVariant(10 * pStackPtb->get_stack_beg_pos().z), _T("Z coordinate of the first ring center."));
        pFirstPosGroup->AddSubItem(pProp);
        pStackRingGroup->AddSubItem(pFirstPosGroup);
// Last ring position:
        CMFCPropertyGridProperty* pLastPosGroup = new CMFCPropertyGridProperty(_T("Last Ring Position, mm"), pStackPtb->get_stack_end_pos_ptr());
        pProp = new CMFCPropertyGridProperty(_T("X"), COleVariant(10 * pStackPtb->get_stack_end_pos().x), _T("X coordinate of the last ring center."));
        pLastPosGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Y"), COleVariant(10 * pStackPtb->get_stack_end_pos().y), _T("Y coordinate of the last ring center."));
        pLastPosGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Z"), COleVariant(10 * pStackPtb->get_stack_end_pos().z), _T("Z coordinate of the last ring center."));
        pLastPosGroup->AddSubItem(pProp);
        pStackRingGroup->AddSubItem(pLastPosGroup);
// Remove the perturbation:
        COleVariant del(_T(""));
        CRemovePerturbationButton* pRemBtn = new CRemovePerturbationButton(this, _T("Remove Perturbation"), del, _T("Click to delete this perturbation."), (DWORD_PTR)pPtb);
        pStackRingGroup->AddSubItem(pRemBtn);

        m_wndPropList.AddProperty(pStackRingGroup);
        break;
      }
      case CFieldPerturbation::ptbUniform:
      {
        CUniformAddField* pAddField = (CUniformAddField*)pPtb;
        CMFCPropertyGridProperty* pAddFieldGroup = new CMFCPropertyGridProperty(pAddField->name());
// Enable:
        CCheckBoxButton* pCheckBox = new CCheckBoxButton(this, _T("Enable"), (_variant_t)pAddField->get_enable(), _T("Turns ON/OFF the field perturbation."), pAddField->get_enable_ptr());
        pAddFieldGroup->AddSubItem(pCheckBox);
// Additional Ex parameters:
        CMFCPropertyGridProperty* pComponentsGroup = new CMFCPropertyGridProperty(_T("Additional Field, V/cm"), pAddField->get_add_Edc_ptr());
        CMFCPropertyGridProperty* pProp = new CMFCPropertyGridProperty(_T("Ex"), COleVariant(pAddField->get_add_Edc().x / SI_to_CGS_Voltage), _T("X-component of the additional electric field."));
        pComponentsGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Ey"), COleVariant(pAddField->get_add_Edc().y / SI_to_CGS_Voltage), _T("Y-component of the additional electric field."));
        pComponentsGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Ez"), COleVariant(pAddField->get_add_Edc().z / SI_to_CGS_Voltage), _T("Z-component of the additional electric field."));
        pComponentsGroup->AddSubItem(pProp);
        pAddFieldGroup->AddSubItem(pComponentsGroup);
// Restriction group:
        CMFCPropertyGridProperty* pLimitsGroup = new CMFCPropertyGridProperty(_T("Limits"));
        pProp = new CMFCPropertyGridProperty(_T("X min, mm"), COleVariant(10 * pAddField->get_add_Edc_beg_x()), _T("The additional electric field is applied in the range from Xmin to Xmax."), pAddField->get_add_Edc_beg_x_ptr());
        pLimitsGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("X max, mm"), COleVariant(10 * pAddField->get_add_Edc_end_x()), _T("The additional electric field is applied in the range from Xmin to Xmax."), pAddField->get_add_Edc_end_x_ptr());
        pLimitsGroup->AddSubItem(pProp);
        pAddFieldGroup->AddSubItem(pLimitsGroup);
// Remove the perturbation:
        COleVariant var(_T(""));
        CRemovePerturbationButton* pRemBtn = new CRemovePerturbationButton(this, _T("Remove Perturbation"), var, _T("Click to delete this perturbation."), (DWORD_PTR)pPtb);
        pAddFieldGroup->AddSubItem(pRemBtn);

        m_wndPropList.AddProperty(pAddFieldGroup);
        break;
      }
      case CFieldPerturbation::ptbDoubleLayer:
      {
        CDoubleLayerField* pDblLayer = (CDoubleLayerField*)pPtb;
        CMFCPropertyGridProperty* pDblLayerGroup = new CMFCPropertyGridProperty(pDblLayer->name());
// Enable:
        CCheckBoxButton* pCheckBox = new CCheckBoxButton(this, _T("Enable"), (_variant_t)pDblLayer->get_enable(), _T("Turns ON/OFF the field perturbation."), pDblLayer->get_enable_ptr());
        pDblLayerGroup->AddSubItem(pCheckBox);
// Film thickness ... 
        CMFCPropertyGridProperty* pProp = new CMFCPropertyGridProperty(_T("Film Thickness, mcm"), COleVariant(1e+4 * pDblLayer->get_film_depth()), _T("The thin dielectric film thickness."), pDblLayer->get_film_depth_ptr());
        pDblLayerGroup->AddSubItem(pProp);
// ... and surface charge density:
        pProp = new CMFCPropertyGridProperty(_T("Surface Charge, nA*hour/mm2"), COleVariant(pDblLayer->get_charge_srf_dens() * Const_Srf_Charge_Dens), _T("The surface charge density in elem. charges per square millimeter."), pDblLayer->get_charge_srf_dens_ptr());
        pDblLayerGroup->AddSubItem(pProp);
// Enable/disable multithreading:
        pCheckBox = new CCheckBoxButton(this, _T("Enable Multithreading"), (_variant_t)pDblLayer->get_enable_multi_thread(), _T("Click this check-box to enable/disable multithreading during the double layer field calculation."), pDblLayer->get_enable_multi_thread_ptr());
        pDblLayerGroup->AddSubItem(pCheckBox);

// Faces selection:
        CMFCPropertyGridProperty* pFacesSelGroup = new CMFCPropertyGridProperty(_T("Contaminated Area"));
        CSelectFacesButton* pSelFaceBtn = new CSelectFacesButton(this, _T("Select Faces"), _T(""), _T("Click to enter the faces selection context. Select proper faces in the draw window and click this button again to exit the context and to confirm the selection."), (DWORD_PTR)pPtb);
        pSelFaceBtn->SetValue(pSelFaceBtn->ButtonValue());
        pFacesSelGroup->AddSubItem(pSelFaceBtn);

        EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
        CRedrawCheckBox* pEnableDrawBtn = new CRedrawCheckBox(this, _T("Enable Faces Drawing"), (_variant_t)pDrawObj->get_enable_sel_faces(), _T("Click to enable/disable selected faces drawing."), pDrawObj->get_enable_sel_faces_ptr());
        pFacesSelGroup->AddSubItem(pEnableDrawBtn);

        CString cDummy(_T(" "));
// Clear selection of all faces:
        CClearSelectedFacesButton* pClearSelBtn = new CClearSelectedFacesButton(this, _T("Clear Faces Selection"), cDummy, _T("Click this button to deselect all previously tagged faces."), (DWORD_PTR)pPtb);
        pFacesSelGroup->AddSubItem(pClearSelBtn);
// Save and load selected faces:
        CSaveLoadSelFacesButton* pSaveFacesBtn = new CSaveLoadSelFacesButton(true, this, _T("Save Selected Faces"), cDummy, _T("Click this button to save all selected faces."), (DWORD_PTR)pPtb);
        pFacesSelGroup->AddSubItem(pSaveFacesBtn);
        CSaveLoadSelFacesButton* pLoadFacesBtn = new CSaveLoadSelFacesButton(false, this, _T("Import Selected Faces"), cDummy, _T("Click this button to import previously saved faces from a file. Note: all current selections will be cleaned."), (DWORD_PTR)pPtb);
        pFacesSelGroup->AddSubItem(pLoadFacesBtn);
        pDblLayerGroup->AddSubItem(pFacesSelGroup);
// Actions:
        CMFCPropertyGridProperty* pActionsGroup = new CMFCPropertyGridProperty(_T("Actions"));
// Calculate field button:
        CCalcFieldPtbButton* pCalcPtb = new CCalcFieldPtbButton(this, _T("Calculate Field"), cDummy, _T("Click this button to start perturbation field calculation."), (DWORD_PTR)pPtb);
        pCalcPtb->Enable(pPtb->get_enable());
        pActionsGroup->AddSubItem(pCalcPtb);
// Remove perturbation:
        CRemovePerturbationButton* pRemBtn = new CRemovePerturbationButton(this, _T("Remove Perturbation"), cDummy, _T("Click to delete this perturbation."), (DWORD_PTR)pPtb);
        pActionsGroup->AddSubItem(pRemBtn);

        pDblLayerGroup->AddSubItem(pActionsGroup);

        m_wndPropList.AddProperty(pDblLayerGroup);
        break;
      }
      case CFieldPerturbation::ptbFlatChannelRF:
      case CFieldPerturbation::ptbCylSubstrateRF:
      case CFieldPerturbation::ptbElliptSubstrRF:
      {
        CAnalytRFField* pAnalytRF = (CAnalytRFField*)pPtb;
        CMFCPropertyGridProperty* pAnalytRFGroup = new CMFCPropertyGridProperty(pAnalytRF->name());
// Enable:
        CCheckBoxButton* pCheckBox = new CCheckBoxButton(this, _T("Enable"), (_variant_t)pAnalytRF->get_enable(), _T("Turns ON/OFF the field perturbation."), pAnalytRF->get_enable_ptr());
        pAnalytRFGroup->AddSubItem(pCheckBox);
// Amplitude and frequency:
        CMFCPropertyGridProperty* pProp = new CMFCPropertyGridProperty(_T("Amplitude, V"), COleVariant(pAnalytRF->get_rf_ampl() / SI_to_CGS_Voltage), _T("Zero-to-Peak radio-frequency amplitude."), pAnalytRF->get_rf_ampl_ptr());
        pAnalytRFGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Frequency, kHz"), COleVariant(0.001 * pAnalytRF->get_rf_freq()), _T("Frequency of simulated RF field in kHz."), pAnalytRF->get_rf_freq_ptr());
        pAnalytRFGroup->AddSubItem(pProp);
// Size of PCB stripes and Fourier decomposition parameters:
        CMFCPropertyGridProperty* pPCBGroup = new CMFCPropertyGridProperty(_T("Parameters of PCB stripes"));
        pProp = new CMFCPropertyGridProperty(_T("Stripe Width, mm"), COleVariant(10 * pAnalytRF->get_stripe_width()), _T("Width of a virtual stripe in mm. Note: the character width of the MEMS stripes is tens of microns."), pAnalytRF->get_stripe_width_ptr());
        pPCBGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("Gap Width, mm"), COleVariant(10 * pAnalytRF->get_gap_width()), _T("Width of gaps between stripes in mm. Note: the character width of the MEMS stripes is tens of microns."), pAnalytRF->get_gap_width_ptr());
        pPCBGroup->AddSubItem(pProp);

        CString sDir(pAnalytRF->get_trans_dir_name(pAnalytRF->get_trans_dir()));
        pProp = new CMFCPropertyGridProperty(_T("Stripes Stretching Dir."), COleVariant(sDir), _T("Translational direction, along which the virtual stripes are stretched."), pAnalytRF->get_trans_dir_ptr());
        for(int i = 0; i < CAnalytRFField::strCount; i++)
          pProp->AddOption(pAnalytRF->get_trans_dir_name(i));

        pPCBGroup->AddSubItem(pProp);

        pProp = new CMFCPropertyGridProperty(_T("Count of Terms in Fourier Decomposition"), COleVariant(long(pAnalytRF->get_decomp_count())), _T("Count of terms in the Fourier decomposition of the bottom boundary function."), pAnalytRF->get_decomp_count_ptr());
        pPCBGroup->AddSubItem(pProp);
        pAnalytRFGroup->AddSubItem(pPCBGroup);

// Round cylinder parameters:
        if(nType == CFieldPerturbation::ptbCylSubstrateRF)
        {
          CCurvedSubstrateRF* pCylSubPtb = (CCurvedSubstrateRF*)pAnalytRF;
          CMFCPropertyGridProperty* pCylGroup = new CMFCPropertyGridProperty(_T("Round Cylinder Parameters"));

          CMFCPropertyGridProperty* pCylx0y0Group = new CMFCPropertyGridProperty(_T("Cylinder Axis Coordinates (x0, y0) and Cylinder Radius"));
          pProp = new CMFCPropertyGridProperty(_T(" x0, mm"), COleVariant(10 * pCylSubPtb->get_x0()), _T("X-coordinate of the ellipse center."), pCylSubPtb->get_x0_ptr());
          pCylx0y0Group->AddSubItem(pProp);
          pProp = new CMFCPropertyGridProperty(_T(" y0, mm"), COleVariant(10 * pCylSubPtb->get_y0()), _T("Y-coordinate of the ellipse center."), pCylSubPtb->get_y0_ptr());
          pCylx0y0Group->AddSubItem(pProp);
          pProp = new CMFCPropertyGridProperty(_T(" Radius, mm"), COleVariant(10 * pCylSubPtb->get_radius()), _T("Y-coordinate of the ellipse center."), pCylSubPtb->get_radius_ptr());
          pCylx0y0Group->AddSubItem(pProp);
          pCylGroup->AddSubItem(pCylx0y0Group);

          pAnalytRFGroup->AddSubItem(pCylGroup);
        }

// Ellipse parameters:
        if(nType == CFieldPerturbation::ptbElliptSubstrRF)
        {
          CEllipticalSubstrateRF* pElliptSubPtb = (CEllipticalSubstrateRF*)pAnalytRF;
          CMFCPropertyGridProperty* pEllipseGroup = new CMFCPropertyGridProperty(_T("Elliptical Shape Parameters"));

          CMFCPropertyGridProperty* pEllipsex0y0Group = new CMFCPropertyGridProperty(_T("Origin of Coordinates x0 and y0"));
          pProp = new CMFCPropertyGridProperty(_T(" x0, mm"), COleVariant(10 * pElliptSubPtb->get_x0()), _T("X-coordinate of the ellipse center."), pElliptSubPtb->get_x0_ptr());
          pEllipsex0y0Group->AddSubItem(pProp);
          pProp = new CMFCPropertyGridProperty(_T(" y0, mm"), COleVariant(10 * pElliptSubPtb->get_y0()), _T("Y-coordinate of the ellipse center."), pElliptSubPtb->get_y0_ptr());
          pEllipsex0y0Group->AddSubItem(pProp);
          pEllipseGroup->AddSubItem(pEllipsex0y0Group);

          CMFCPropertyGridProperty* pEllipseABGroup = new CMFCPropertyGridProperty(_T("Ellipse Semi-Axses a and b"));
          pProp = new CMFCPropertyGridProperty(_T(" a (along x), mm"), COleVariant(10 * pElliptSubPtb->get_a()), _T("The X-semi-axis of the ellipse."), pElliptSubPtb->get_a_ptr());
          pEllipseABGroup->AddSubItem(pProp);
          pProp = new CMFCPropertyGridProperty(_T(" b (along y), mm"), COleVariant(10 * pElliptSubPtb->get_b()), _T("The Y-semi-axis of the ellipse."), pElliptSubPtb->get_b_ptr());
          pEllipseABGroup->AddSubItem(pProp);
          pEllipseGroup->AddSubItem(pEllipseABGroup);
          
          pAnalytRFGroup->AddSubItem(pEllipseGroup);
        }

// Substrate selection:
        CMFCPropertyGridProperty* pSubstrRegsGroup = new CMFCPropertyGridProperty(_T("Substrate Regions"), pAnalytRF->get_region_names_ptr());

// Manual selection:
        CString cRegNames = EvaporatingParticle::CObject::compile_string(pAnalytRF->get_region_names());
        CSelectRegionButton* pSelRegButton = new CSelectRegionButton(this, _T("Select Regions Manually"), cRegNames, _T("Click to select 2D regions in the main view window."), pAnalytRF->get_region_names_ptr());
        pSubstrRegsGroup->AddSubItem(pSelRegButton);

// Merging with Named Areas:
        CNamedAreasSelResponder* pNamedAreasSelector = new CNamedAreasSelResponder(this, _T("Merge with Named Areas"), pAnalytRF->get_last_merged(), _T("Select the existing Named Areas to use these surfaces for boundary condtions setting."), pAnalytRF->get_region_names_ptr());
        pNamedAreasSelector->AllowEdit(FALSE);
        pNamedAreasSelector->AddOption(_T("None"));
        EvaporatingParticle::CSelAreasColl* pSelAreasColl = CParticleTrackingApp::Get()->GetSelAreas();
        size_t nSelAreasCount = pSelAreasColl->size();
        for(size_t k = 0; k < nSelAreasCount; k++)
          pNamedAreasSelector->AddOption(pSelAreasColl->at(k)->get_name());

        pSubstrRegsGroup->AddSubItem(pNamedAreasSelector);

// Merging options:
        CString sMergeOpt = CSelectedAreas::merge_opt_name(pAnalytRF->get_merge_option());
        CMFCPropertyGridProperty* pMergeOptSelector = new CMFCPropertyGridProperty(_T("Merge Options"), sMergeOpt, _T("Select one of three allowed merge options: add, substitute and subtract."), pAnalytRF->get_merge_option_ptr());
        for(int l = EvaporatingParticle::CSelectedAreas::optAdd; l < EvaporatingParticle::CSelectedAreas::optCount; l++)
          pMergeOptSelector->AddOption(EvaporatingParticle::CSelectedAreas::merge_opt_name(l));

        pSubstrRegsGroup->AddSubItem(pMergeOptSelector);
      
// Hide/Show selected regions:
        CHideShowRegsCheckBox* pHideShowBtn = new CHideShowRegsCheckBox(this, _T("Visibility"), (_variant_t)pAnalytRF->get_visibility_flag(), _T("Change the visibility status of the selected regions"), pAnalytRF->get_visibility_flag_ptr());
        pSubstrRegsGroup->AddSubItem(pHideShowBtn);
        pAnalytRFGroup->AddSubItem(pSubstrRegsGroup);

// Analytical domain size:
        CMFCPropertyGridProperty* pSizeGroup = new CMFCPropertyGridProperty(_T("Size of Analytical Domain"));
        pProp = new CMFCPropertyGridProperty(_T("Height, mm"), COleVariant(10 * pAnalytRF->get_channel_height()), _T("Domain height in Y direction."), pAnalytRF->get_channel_height_ptr());
        pSizeGroup->AddSubItem(pProp);
        pAnalytRFGroup->AddSubItem(pSizeGroup);
// Remove perturbation:
        CRemovePerturbationButton* pRemBtn = new CRemovePerturbationButton(this, _T("Remove Perturbation"), _T(" "), _T("Click to delete this perturbation."), (DWORD_PTR)pPtb);
        pAnalytRFGroup->AddSubItem(pRemBtn);

        m_wndPropList.AddProperty(pAnalytRFGroup);
        break;
      }
    }
  }
}

void CPropertiesWnd::set_ptb_data()
{
  CFieldPtbCollection& coll = CParticleTrackingApp::Get()->GetTracker()->get_field_ptb();
  size_t nPtbCount = coll.size();
  for(size_t i = 0; i < nPtbCount; i++)
  {
    CFieldPerturbation* pPtb = coll.at(i);
    int nType = pPtb->type();
    switch(nType)
    {
      case CFieldPerturbation::ptbRing:
      {
        CChargedRingPerturbation* pChargedRing = (CChargedRingPerturbation*)pPtb;

        Vector3D vPos;
        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pChargedRing->get_ring_pos_ptr());
        if(pProp != NULL)
        {
          vPos.x = 0.1 * pProp->GetSubItem(0)->GetValue().dblVal;
          vPos.y = 0.1 * pProp->GetSubItem(1)->GetValue().dblVal;
          vPos.z = 0.1 * pProp->GetSubItem(2)->GetValue().dblVal;
          pChargedRing->set_ring_pos(vPos);
        }

        pProp = m_wndPropList.FindItemByData(pChargedRing->get_ring_charge_ptr());
        if(pProp != NULL)
          pChargedRing->set_ring_charge(1e+6 * Const_Charge_CGS * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pChargedRing->get_ring_radius_ptr());
        if(pProp != NULL)
          pChargedRing->set_ring_radius(0.1 * pProp->GetValue().dblVal);

        break;
      }
      case CFieldPerturbation::ptbStackOfRings:
      {
        CStackRingPerturbation* pStackPtb = (CStackRingPerturbation*)pPtb;

        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pStackPtb->get_rings_count_ptr());
        if(pProp != NULL)
          pStackPtb->set_rings_count(pProp->GetValue().lVal);

        pProp = m_wndPropList.FindItemByData(pStackPtb->get_ring_radius_ptr());
        if(pProp != NULL)
          pStackPtb->set_ring_radius(0.1 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pStackPtb->get_sum_charge_ptr());
        if(pProp != NULL)
          pStackPtb->set_sum_charge(1e+6 * Const_Charge_CGS * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pStackPtb->get_charge_distr_type_ptr());
        if(pProp != NULL)
        {
          CString cType = (CString)pProp->GetValue();
          for(int i = 0; i < EvaporatingParticle::CStackRingPerturbation::distrCount; i++)
          {
            if(cType == pStackPtb->get_distr_type_name(i))
            {
              pStackPtb->set_charge_distr_type(i);
              break;
            }
          }
        }

        Vector3D vPos;
        pProp = m_wndPropList.FindItemByData(pStackPtb->get_stack_beg_pos_ptr());
        if(pProp != NULL)
        {
          vPos.x = 0.1 * pProp->GetSubItem(0)->GetValue().dblVal;
          vPos.y = 0.1 * pProp->GetSubItem(1)->GetValue().dblVal;
          vPos.z = 0.1 * pProp->GetSubItem(2)->GetValue().dblVal;
          pStackPtb->set_stack_beg_pos(vPos);
        }

        pProp = m_wndPropList.FindItemByData(pStackPtb->get_stack_end_pos_ptr());
        if(pProp != NULL)
        {
          vPos.x = 0.1 * pProp->GetSubItem(0)->GetValue().dblVal;
          vPos.y = 0.1 * pProp->GetSubItem(1)->GetValue().dblVal;
          vPos.z = 0.1 * pProp->GetSubItem(2)->GetValue().dblVal;
          pStackPtb->set_stack_end_pos(vPos);
        }
      }
      case CFieldPerturbation::ptbUniform:
      {
        CUniformAddField* pAddField = (CUniformAddField*)pPtb;

        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pAddField->get_add_Edc_ptr());
        if(pProp != NULL)
        {
          EvaporatingParticle::Vector3D vE;
          vE.x = SI_to_CGS_Voltage * pProp->GetSubItem(0)->GetValue().dblVal;
          vE.y = SI_to_CGS_Voltage * pProp->GetSubItem(1)->GetValue().dblVal;
          vE.z = SI_to_CGS_Voltage * pProp->GetSubItem(2)->GetValue().dblVal;
          pAddField->set_add_Edc(vE);
        }

        pProp = m_wndPropList.FindItemByData(pAddField->get_add_Edc_beg_x_ptr());
        if(pProp != NULL)
          pAddField->set_add_Edc_beg_x(0.1 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pAddField->get_add_Edc_end_x_ptr());
        if(pProp != NULL)
          pAddField->set_add_Edc_end_x(0.1 * pProp->GetValue().dblVal);

        break;
      }
      case CFieldPerturbation::ptbDoubleLayer:
      {
        CDoubleLayerField* pDblLayer = (CDoubleLayerField*)pPtb;

        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pDblLayer->get_film_depth_ptr());
        if(pProp != NULL)
          pDblLayer->set_film_depth(1e-4 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pDblLayer->get_charge_srf_dens_ptr());
        if(pProp != NULL)
          pDblLayer->set_charge_srf_dens(pProp->GetValue().dblVal / Const_Srf_Charge_Dens);

        break;
      }
      case CFieldPerturbation::ptbFlatChannelRF:
      case CFieldPerturbation::ptbCylSubstrateRF:
      case CFieldPerturbation::ptbElliptSubstrRF:
      {
        CAnalytRFField* pAnalytRF = (CAnalytRFField*)pPtb;
        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pAnalytRF->get_rf_ampl_ptr());
        if(pProp != NULL)
          pAnalytRF->set_rf_ampl(SI_to_CGS_Voltage * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_rf_freq_ptr());
        if(pProp != NULL)
          pAnalytRF->set_rf_freq(1000 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_stripe_width_ptr());
        if(pProp != NULL)
          pAnalytRF->set_stripe_width(0.1 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_gap_width_ptr());
        if(pProp != NULL)
          pAnalytRF->set_gap_width(0.1 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_decomp_count_ptr());
        if(pProp != NULL)
          pAnalytRF->set_decomp_count(pProp->GetValue().lVal);

        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_channel_height_ptr());
        if(pProp != NULL)
          pAnalytRF->set_channel_height(0.1 * pProp->GetValue().dblVal);
        
        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_trans_dir_ptr());
        if(pProp != NULL)
        {
          CString sDir = (CString)pProp->GetValue();
          for(int i = 0; i < CAnalytRFField::strCount; i++)
          {
            if((CString)pProp->GetOption(i) == sDir)
            {
              pAnalytRF->set_trans_dir(i);
              break;
            }
          }
        }

        if(nType == CFieldPerturbation::ptbCylSubstrateRF)
        {
          CCurvedSubstrateRF* pCylSubPtb = (CCurvedSubstrateRF*)pAnalytRF;
          pProp = m_wndPropList.FindItemByData(pCylSubPtb->get_x0_ptr());
          if(pProp != NULL)
            pCylSubPtb->set_x0(0.1 * pProp->GetValue().dblVal);

          pProp = m_wndPropList.FindItemByData(pCylSubPtb->get_y0_ptr());
          if(pProp != NULL)
            pCylSubPtb->set_y0(0.1 * pProp->GetValue().dblVal);

          pProp = m_wndPropList.FindItemByData(pCylSubPtb->get_radius_ptr());
          if(pProp != NULL)
            pCylSubPtb->set_radius(0.1 * pProp->GetValue().dblVal);
        }

        if(nType == CFieldPerturbation::ptbElliptSubstrRF)
        {
          CEllipticalSubstrateRF* pElliptSubPtb = (CEllipticalSubstrateRF*)pAnalytRF;
          pProp = m_wndPropList.FindItemByData(pElliptSubPtb->get_x0_ptr());
          if(pProp != NULL)
            pElliptSubPtb->set_x0(0.1 * pProp->GetValue().dblVal);

          pProp = m_wndPropList.FindItemByData(pElliptSubPtb->get_y0_ptr());
          if(pProp != NULL)
            pElliptSubPtb->set_y0(0.1 * pProp->GetValue().dblVal);

          pProp = m_wndPropList.FindItemByData(pElliptSubPtb->get_a_ptr());
          if(pProp != NULL)
            pElliptSubPtb->set_a(0.1 * pProp->GetValue().dblVal);

          pProp = m_wndPropList.FindItemByData(pElliptSubPtb->get_b_ptr());
          if(pProp != NULL)
            pElliptSubPtb->set_b(0.1 * pProp->GetValue().dblVal);
        }

// Merge options:
        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_merge_option_ptr());
        if(pProp != NULL)
        {
          CString sOptSel = pProp->GetValue();
          for(int i = CSelectedAreas::optAdd; i < CSelectedAreas::optCount; i++)
          {
            if(sOptSel == CSelectedAreas::merge_opt_name(i))
            {
              pAnalytRF->set_merge_option(i);
              break;
            }
          }
        }

        break;
      }
    }
  }
}

void CPropertiesWnd::update_ptb_ctrls()
{
  CFieldPtbCollection& coll = CParticleTrackingApp::Get()->GetTracker()->get_field_ptb();
  size_t nPtbCount = coll.size();
  for(size_t i = 0; i < nPtbCount; i++)
  {
    CFieldPerturbation* pPtb = coll.at(i);
    int nType = pPtb->type();
    bool bEnable = pPtb->get_enable();
    switch(nType)
    {
      case CFieldPerturbation::ptbRing:
      {
        CChargedRingPerturbation* pChargedRing = (CChargedRingPerturbation*)pPtb;
        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pChargedRing->get_ring_pos_ptr());
        if(pProp != NULL)
        {
          pProp->GetSubItem(0)->Enable(bEnable);
          pProp->GetSubItem(1)->Enable(bEnable);
          pProp->GetSubItem(2)->Enable(bEnable);
        }

        pProp = m_wndPropList.FindItemByData(pChargedRing->get_ring_charge_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        pProp = m_wndPropList.FindItemByData(pChargedRing->get_ring_radius_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        break;
      }
      case CFieldPerturbation::ptbStackOfRings:
      {
        CStackRingPerturbation* pStackPtb = (CStackRingPerturbation*)pPtb;
        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pStackPtb->get_rings_count_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        pProp = m_wndPropList.FindItemByData(pStackPtb->get_ring_radius_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        pProp = m_wndPropList.FindItemByData(pStackPtb->get_sum_charge_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        pProp = m_wndPropList.FindItemByData(pStackPtb->get_charge_distr_type_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        Vector3D vPos;
        pProp = m_wndPropList.FindItemByData(pStackPtb->get_stack_beg_pos_ptr());
        if(pProp != NULL)
        {
          pProp->GetSubItem(0)->Enable(bEnable);
          pProp->GetSubItem(1)->Enable(bEnable);
          pProp->GetSubItem(2)->Enable(bEnable);
        }

        pProp = m_wndPropList.FindItemByData(pStackPtb->get_stack_end_pos_ptr());
        if(pProp != NULL)
        {
          pProp->GetSubItem(0)->Enable(bEnable);
          pProp->GetSubItem(1)->Enable(bEnable);
          pProp->GetSubItem(2)->Enable(bEnable);
        }
      }
      case CFieldPerturbation::ptbUniform:
      {
        CUniformAddField* pAddField = (CUniformAddField*)pPtb;
        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pAddField->get_add_Edc_ptr());
        if(pProp != NULL)
        {
          pProp->GetSubItem(0)->Enable(bEnable);
          pProp->GetSubItem(1)->Enable(bEnable);
          pProp->GetSubItem(2)->Enable(bEnable);
        }

        pProp = m_wndPropList.FindItemByData(pAddField->get_add_Edc_beg_x_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        pProp = m_wndPropList.FindItemByData(pAddField->get_add_Edc_end_x_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        break;
      }
      case CFieldPerturbation::ptbFlatChannelRF:
      case CFieldPerturbation::ptbCylSubstrateRF:
      case CFieldPerturbation::ptbElliptSubstrRF:
      {
        CAnalytRFField* pAnalytRF = (CAnalytRFField*)pPtb;
        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pAnalytRF->get_rf_ampl_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_rf_freq_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_stripe_width_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_gap_width_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_decomp_count_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_channel_height_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);
        
        pProp = m_wndPropList.FindItemByData(pAnalytRF->get_trans_dir_ptr());
        if(pProp != NULL)
          pProp->Enable(bEnable);

        if(nType == CFieldPerturbation::ptbCylSubstrateRF)
        {
          CCurvedSubstrateRF* pCylSubPtb = (CCurvedSubstrateRF*)pAnalytRF;
          pProp = m_wndPropList.FindItemByData(pCylSubPtb->get_x0_ptr());
          if(pProp != NULL)
            pProp->Enable(bEnable);

          pProp = m_wndPropList.FindItemByData(pCylSubPtb->get_y0_ptr());
          if(pProp != NULL)
            pProp->Enable(bEnable);

          pProp = m_wndPropList.FindItemByData(pCylSubPtb->get_radius_ptr());
          if(pProp != NULL)
            pProp->Enable(bEnable);
        }

        if(nType == CFieldPerturbation::ptbElliptSubstrRF)
        {
          CEllipticalSubstrateRF* pElliptSubPtb = (CEllipticalSubstrateRF*)pAnalytRF;
          pProp = m_wndPropList.FindItemByData(pElliptSubPtb->get_x0_ptr());
          if(pProp != NULL)
            pProp->Enable(bEnable);

          pProp = m_wndPropList.FindItemByData(pElliptSubPtb->get_y0_ptr());
          if(pProp != NULL)
            pProp->Enable(bEnable);

          pProp = m_wndPropList.FindItemByData(pElliptSubPtb->get_a_ptr());
          if(pProp != NULL)
            pProp->Enable(bEnable);

          pProp = m_wndPropList.FindItemByData(pElliptSubPtb->get_b_ptr());
          if(pProp != NULL)
            pProp->Enable(bEnable);
        }

        break;
      }
    }
  }
}
