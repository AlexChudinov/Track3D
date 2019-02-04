
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
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
        CSelectFacesButton* pSelFaceBtn = new CSelectFacesButton(this, _T("Select Faces"), _T(""), _T("Click to enter the faces selection context. Select proper faces in the draw window and click this button again to exit the confirm the selection."), (DWORD_PTR)pPtb);
        pSelFaceBtn->SetValue(pSelFaceBtn->ButtonValue());
        pFacesSelGroup->AddSubItem(pSelFaceBtn);

        EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
        CRedrawCheckBox* pEnableDrawBtn = new CRedrawCheckBox(this, _T("Enable Faces Drawing"), (_variant_t)pDrawObj->get_enable_sel_faces(), _T("Click to enable/disable selected faces drawing."), pDrawObj->get_enable_sel_faces_ptr());
        pFacesSelGroup->AddSubItem(pEnableDrawBtn);

        CString cDummy(_T(" "));
// Clear selection of all faces:
        CClearSelectedFacesButton* pClearSelBtn = new CClearSelectedFacesButton(this, _T("Clear Faces Selection"), cDummy, _T("Click this button to deselect all previously tagged faces."), (DWORD_PTR)pPtb);
        pFacesSelGroup->AddSubItem(pClearSelBtn);
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
    }
  }
}
