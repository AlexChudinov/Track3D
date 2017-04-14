
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "Button.h"

using namespace EvaporatingParticle;

void CPropertiesWnd::add_ptb_ctrls()
{
  CMFCPropertyGridProperty* pFieldPtbGroup = new CMFCPropertyGridProperty(_T("DC Field Perturbations"));

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

        pFieldPtbGroup->AddSubItem(pChargedRingGroup);
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

        pFieldPtbGroup->AddSubItem(pStackRingGroup);
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
        CMFCPropertyGridProperty* pProp = new CMFCPropertyGridProperty(_T("Additional Ex, V/cm"), COleVariant(pAddField->get_add_Edc() / SI_to_CGS_Voltage), _T("X-component of the additional electric field."), pAddField->get_add_Edc_ptr());
        pAddFieldGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("X min, mm"), COleVariant(10 * pAddField->get_add_Edc_beg_x()), _T("The additional electric field is applied in the range from Xmin to Xmax."), pAddField->get_add_Edc_beg_x_ptr());
        pAddFieldGroup->AddSubItem(pProp);
        pProp = new CMFCPropertyGridProperty(_T("X max, mm"), COleVariant(10 * pAddField->get_add_Edc_end_x()), _T("The additional electric field is applied in the range from Xmin to Xmax."), pAddField->get_add_Edc_end_x_ptr());
        pAddFieldGroup->AddSubItem(pProp);
// Remove the perturbation:
        COleVariant var(_T(""));
        CRemovePerturbationButton* pRemBtn = new CRemovePerturbationButton(this, _T("Remove Perturbation"), var, _T("Click to delete this perturbation."), (DWORD_PTR)pPtb);
        pAddFieldGroup->AddSubItem(pRemBtn);

        pFieldPtbGroup->AddSubItem(pAddFieldGroup);
        break;
      }
    }
  }

  m_wndPropList.AddProperty(pFieldPtbGroup);
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

        CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pAddField->get_enable_ptr());
        if(pProp != NULL)
          pAddField->set_enable(pProp->GetValue().boolVal);

        pProp = m_wndPropList.FindItemByData(pAddField->get_add_Edc_ptr());
        if(pProp != NULL)
          pAddField->set_add_Edc(SI_to_CGS_Voltage * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pAddField->get_add_Edc_beg_x_ptr());
        if(pProp != NULL)
          pAddField->set_add_Edc_beg_x(0.1 * pProp->GetValue().dblVal);

        pProp = m_wndPropList.FindItemByData(pAddField->get_add_Edc_end_x_ptr());
        if(pProp != NULL)
          pAddField->set_add_Edc_end_x(0.1 * pProp->GetValue().dblVal);

        break;
      }
    }
  }
}

void CPropertiesWnd::update_ptb_ctrls()
{
}
