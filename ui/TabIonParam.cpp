
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "BeamCrossSection.h"
#include "Button.h"

static const double scfA2 = Const_Angstrem_CGS * Const_Angstrem_CGS;
//---------------------------------------------------------------------------------------
// Ion type of particles
//---------------------------------------------------------------------------------------
void CPropertiesWnd::add_ion_ctrls()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CMFCPropertyGridProperty* pProp = NULL;

  CMFCPropertyGridProperty* pMassGroup = new CMFCPropertyGridProperty(_T("Mass, a.m.u."));
  pProp = new CMFCPropertyGridProperty(_T("Mass"), COleVariant(long(pObj->get_ion_mass() / Const_AMU_CGS)), _T("Ion mass, atomic mass units."), pObj->get_ion_mass_ptr());
  pMassGroup->AddSubItem(pProp);
  m_wndPropList.AddProperty(pMassGroup);

  CMFCPropertyGridProperty* pMobGroup = new CMFCPropertyGridProperty(_T("Mobility, cm2/s/V"));
  pProp = new CMFCPropertyGridProperty(_T("Mobility"), COleVariant(pObj->get_ion_mobility() * SI_to_CGS_Voltage), _T("Ion mobility at STP."), pObj->get_ion_mobility_ptr());
  pMobGroup->AddSubItem(pProp);
  m_wndPropList.AddProperty(pMobGroup);

  CMFCPropertyGridProperty* pCollisionGroup = new CMFCPropertyGridProperty(_T("Collision Parameters"));
  pProp = new CMFCPropertyGridProperty(_T("Cross-Section"), COleVariant(pObj->get_ion_cross_section() / scfA2), _T("Collision cross-section with the environmental gas, squared angstrem."), pObj->get_ion_cross_section_ptr());
  pCollisionGroup->AddSubItem(pProp);

  CCheckBoxButton* pCheckBox = new CCheckBoxButton(this, _T("Velocity-dependent"), (_variant_t)pObj->get_vel_depend_flag(), _T("If this is set to 'true' the cross-section is inversely proportional to the relative ion - air velocity."), pObj->get_vel_depend_flag_ptr());
  pCollisionGroup->AddSubItem(pCheckBox);

  pProp = new CMFCPropertyGridProperty(_T("Activation Energy, eV"), COleVariant(pObj->get_act_energy()), _T("Ion collisional activation energy, eV"), pObj->get_act_energy_ptr());
  pCollisionGroup->AddSubItem(pProp);
  m_wndPropList.AddProperty(pCollisionGroup);

// Coulomb:
  CMFCPropertyGridProperty* pCoulombGroup = new CMFCPropertyGridProperty(_T("Coulomb Repulsion"));
  pCheckBox = new CCheckBoxButton(this, _T("Enable Coulomb"), (_variant_t)pObj->get_enable_coulomb(), _T("Turns ON/OFF the Coulomb repulsion term in the ion momentum equation."), pObj->get_enable_coulomb_ptr());
  pCoulombGroup->AddSubItem(pCheckBox);
  pProp = new CMFCPropertyGridProperty(_T("Full Ion Current, nA"), COleVariant(pObj->get_full_current() / Const_nA_to_CGSE), _T("Full current carried by the ion bunch."), pObj->get_full_current_ptr());
  pCoulombGroup->AddSubItem(pProp);

  pCheckBox = new CCheckBoxButton(this, _T("Axial Symmetry"), (_variant_t)pObj->get_axial_symm(), _T("Must be 'true' for axially-symmetric systems. If set to 'false' the Coulomb repulsion is estimated for a flat model of letter-box capillary."), pObj->get_axial_symm_ptr());
  pCoulombGroup->AddSubItem(pCheckBox);
  pProp = new CMFCPropertyGridProperty(_T("Bunch Radius, mm"), COleVariant(10 * pObj->get_bunch_r0()), _T("Initial radius of the ion bunch."), pObj->get_bunch_r0_ptr());
  pCoulombGroup->AddSubItem(pProp);

// Iterational approach for the Coulomb repulsion computation. The general case (no axial symmetry):
  pProp = new CMFCPropertyGridProperty(_T("Iterations Count"), COleVariant((long)pObj->get_iter_count()), _T("Iterations count in the flux-tube iterational method of Coulomb repulsion calculation"), pObj->get_iter_count_ptr());
  pCoulombGroup->AddSubItem(pProp);

  EvaporatingParticle::CSpaceChargeDistrib& distrib = pObj->get_space_charge_dist();
  pProp = new CMFCPropertyGridProperty(_T("Subdivisions Count"), COleVariant((long)distrib.get_planes_count()), _T("Count of model volumes for the space charge distribution."), distrib.get_planes_count_ptr());
  pCoulombGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Space Charge Step, mm"), COleVariant(10 * distrib.get_space_charge_step()), _T("Distance between pseudo-charges placed into nodes of a constant cubical mesh."), distrib.get_space_charge_step_ptr());
  pCoulombGroup->AddSubItem(pProp);

  pCheckBox = new CCheckBoxButton(this, _T("Use Radial Coulomb"), (_variant_t)pObj->get_use_radial_coulomb(), _T("If this is true the radial Gabovich formula will be used for the space-charge field for x > Transition X."), pObj->get_use_radial_coulomb_ptr());
  pCoulombGroup->AddSubItem(pCheckBox);
  pProp = new CMFCPropertyGridProperty(_T("Radial Coulomb Transition X, mm"), COleVariant(10 * pObj->get_radial_coulomb_trans()), _T("The radial Gabovich formula will be used for the space-charge field for x > Transition X."), pObj->get_radial_coulomb_trans_ptr());
  pCoulombGroup->AddSubItem(pProp);

// Save Coulomb field to a file:
  COleVariant vrn(_T(""));
  CSaveFieldButton* pSaveBtn = new CSaveFieldButton(this, _T("Save Coulomb Field"), vrn, _T("Click this button to save the pre-calculated Coulomb field in a disk file"), (DWORD_PTR)pObj);
  pCoulombGroup->AddSubItem(pSaveBtn);

// Pre-calculated Coulomb field:
  pCheckBox = new CCheckBoxButton(this, _T("Use Pre-calculated Coulomb"), (_variant_t)pObj->get_use_pre_calc_coulomb(), _T("If this is true the Coulomb field will be read from the file specified below."), pObj->get_use_pre_calc_coulomb_ptr());
  pCoulombGroup->AddSubItem(pCheckBox);

// Coulomb data filename:
  static TCHAR BASED_CODE szFilter[] = _T("CSV Files(*.csv)|*.csv|All Files(*.*)|*.*||");
  CMFCPropertyGridFileProperty* pFileProp = new CMFCPropertyGridFileProperty(_T("Coulomb Field Data"), TRUE, pObj->get_pre_calc_clmb_file(), _T("csv"),
    OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, _T("Specify location of the pre-calculated Coulomb field data file."), pObj->get_pre_calc_clmb_file_ptr());

  pCoulombGroup->AddSubItem(pFileProp);

// Barnes-Hut group of controls:
  CMFCPropertyGridProperty* pBHGroup = new CMFCPropertyGridProperty(_T("Barnes-Hut Method Parameters"));
  pProp = new CMFCPropertyGridProperty(_T("Dimensionless Dist. Coeff."), COleVariant(pObj->get_BH_dist_par()), _T("The greater this parameter the farther from a cubical cell the exact Coulomb force is substituted by its approximate expression."), pObj->get_BH_dist_par_ptr());
  pBHGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Critical Radius, mm"), COleVariant(10 * pObj->get_crit_radius()), _T("At distances less than the critical radius the repulsion force linearly changes with the distance to the point charge."), pObj->get_crit_radius_ptr());
  pBHGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Max Recursion Depth"), COleVariant((long)pObj->get_max_rec_depth()), _T("The depth of sub-division of the initial cubical cell into sub-cells."), pObj->get_max_rec_depth_ptr());
  pBHGroup->AddSubItem(pProp);
  pCheckBox = new CCheckBoxButton(this, _T("Enable Quad Terms"), (_variant_t)pObj->get_enable_quad_terms(), _T("If this is set to 'true' the quadripole terms are taken into account when Coulomb force is calculated."), pObj->get_enable_quad_terms_ptr());
  pBHGroup->AddSubItem(pCheckBox);

// Pseudo-ions distribution type in the Barnes-Hut object:
  COleVariant var(EvaporatingParticle::CSpaceChargeDistrib::get_distrib_type_name(distrib.get_ion_distrib_type()));
  pProp = new CMFCPropertyGridProperty(_T("Pseudo-Ions Distrib. Type"), var, _T("The pseudo-ions for the space charge simulation can be distributed either along the trajectories or in the nodes of a constant cubical mesh."), distrib.get_ion_distrib_type_ptr());
  for(UINT i = 0; i < EvaporatingParticle::CSpaceChargeDistrib::distCount; i++)
    pProp->AddOption(EvaporatingParticle::CSpaceChargeDistrib::get_distrib_type_name(i));

  pProp->AllowEdit(FALSE);
  pBHGroup->AddSubItem(pProp);

  pCoulombGroup->AddSubItem(pBHGroup);

  m_wndPropList.AddProperty(pCoulombGroup);
}

void CPropertiesWnd::set_ion_data()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pObj->get_ion_mass_ptr());
  if(pProp != NULL)
    pObj->set_ion_mass(Const_AMU_CGS * pProp->GetValue().lVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_ion_mobility_ptr());
  if(pProp != NULL)
    pObj->set_ion_mobility(pProp->GetValue().dblVal / SI_to_CGS_Voltage);

  pProp = m_wndPropList.FindItemByData(pObj->get_ion_cross_section_ptr());
  if(pProp != NULL)
    pObj->set_ion_cross_section(pProp->GetValue().dblVal * scfA2);

  pProp = m_wndPropList.FindItemByData(pObj->get_vel_depend_flag_ptr());
  if(pProp != NULL)
    pObj->set_vel_depend_flag(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_act_energy_ptr());
  if(pProp != NULL)
    pObj->set_act_energy(pProp->GetValue().dblVal);

// Coulomb:
  pProp = m_wndPropList.FindItemByData(pObj->get_enable_coulomb_ptr());
  if(pProp != NULL)
    pObj->set_enable_coulomb(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_full_current_ptr());
  if(pProp != NULL)
    pObj->set_full_current(Const_nA_to_CGSE * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_bunch_r0_ptr());
  if(pProp != NULL)
    pObj->set_bunch_r0(0.1 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_axial_symm_ptr());
  if(pProp != NULL)
    pObj->set_axial_symm(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_iter_count_ptr());
  if(pProp != NULL)
    pObj->set_iter_count(pProp->GetValue().lVal);

  EvaporatingParticle::CSpaceChargeDistrib& distrib = pObj->get_space_charge_dist();
  pProp = m_wndPropList.FindItemByData(distrib.get_planes_count_ptr());
  if(pProp != NULL)
    distrib.set_planes_count(pProp->GetValue().lVal);

  pProp = m_wndPropList.FindItemByData(distrib.get_space_charge_step_ptr());
  if(pProp != NULL)
    distrib.set_space_charge_step(0.1 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_use_radial_coulomb_ptr());
  if(pProp != NULL)
    pObj->set_use_radial_coulomb(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_radial_coulomb_trans_ptr());
  if(pProp != NULL)
    pObj->set_radial_coulomb_trans(0.1 * pProp->GetValue().dblVal);

// Coulomb data filename:
  pProp = m_wndPropList.FindItemByData(pObj->get_pre_calc_clmb_file_ptr());
  if(pProp != NULL)
  {
    CString cFile = (CString)pProp->GetValue();
    std::string str = std::string(CT2CA(cFile));
    pObj->set_pre_calc_clmb_file(str.c_str());
  }

// Barnes-Hut group of controls: 
  pProp = m_wndPropList.FindItemByData(pObj->get_BH_dist_par_ptr());
  if(pProp != NULL)
    pObj->set_BH_dist_par(pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_crit_radius_ptr());
  if(pProp != NULL)
    pObj->set_crit_radius(0.1 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_max_rec_depth_ptr());
  if(pProp != NULL)
    pObj->set_max_rec_depth(pProp->GetValue().lVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_enable_quad_terms_ptr());
  if(pProp != NULL)
    pObj->set_enable_quad_terms(pProp->GetValue().boolVal);

   pProp = m_wndPropList.FindItemByData(distrib.get_ion_distrib_type_ptr());
  if(pProp != NULL)
  {
    CString cTypeName = (CString)pProp->GetValue();
    for(UINT i = 0; i < EvaporatingParticle::CSpaceChargeDistrib::distCount; i++)
    {
      if(cTypeName == EvaporatingParticle::CSpaceChargeDistrib::get_distrib_type_name(i))
      {
        distrib.set_ion_distrib_type(i);
        break;
      }
    }
  }
}

void CPropertiesWnd::update_ion_ctrls()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  bool bEnable = pObj->get_particle_type() == EvaporatingParticle::CTrack::ptIon;

  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pObj->get_ion_mass_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(pObj->get_ion_mobility_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(pObj->get_ion_cross_section_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(pObj->get_vel_depend_flag_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(pObj->get_act_energy_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

// Coulomb:
  bool bEnableCoulomb = FALSE;
  bool bAxialSymm = TRUE;
  pProp = m_wndPropList.FindItemByData(pObj->get_enable_coulomb_ptr());
  if(pProp != NULL)
  {
    pProp->Enable(bEnable);
    bEnableCoulomb = pProp->GetValue().boolVal;
  }

  pProp = m_wndPropList.FindItemByData(pObj->get_full_current_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable && bEnableCoulomb);

  pProp = m_wndPropList.FindItemByData(pObj->get_axial_symm_ptr());
  if(pProp != NULL)
  {
    pProp->Enable(bEnable && bEnableCoulomb);
    if(bEnable)
      bAxialSymm = pProp->GetValue().boolVal;
  }

  pProp = m_wndPropList.FindItemByData(pObj->get_bunch_r0_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable && bEnableCoulomb && bAxialSymm);

  pProp = m_wndPropList.FindItemByData(pObj->get_iter_count_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable && bEnableCoulomb && !bAxialSymm);

  bool bDistAlongTraject = false;
  EvaporatingParticle::CSpaceChargeDistrib& distrib = pObj->get_space_charge_dist();
  pProp = m_wndPropList.FindItemByData(distrib.get_ion_distrib_type_ptr());
  if(pProp != NULL)
  {
    bool bEnableSelector = bEnable && bEnableCoulomb && !bAxialSymm;
    pProp->Enable(bEnableSelector);
    if(bEnableSelector)
    {
      CString cTypeName = (CString)pProp->GetValue();
      int nAlongTraject = EvaporatingParticle::CSpaceChargeDistrib::distAlongTraject;
      bDistAlongTraject = cTypeName == EvaporatingParticle::CSpaceChargeDistrib::get_distrib_type_name(nAlongTraject);
    }
  }

  pProp = m_wndPropList.FindItemByData(distrib.get_planes_count_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable && bEnableCoulomb && !bAxialSymm && !bDistAlongTraject);

  pProp = m_wndPropList.FindItemByData(distrib.get_space_charge_step_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable && bEnableCoulomb && !bAxialSymm);

  pProp = m_wndPropList.FindItemByData(pObj->get_BH_dist_par_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable && bEnableCoulomb && !bAxialSymm);

  pProp = m_wndPropList.FindItemByData(pObj->get_crit_radius_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable && bEnableCoulomb && !bAxialSymm);

  pProp = m_wndPropList.FindItemByData(pObj->get_max_rec_depth_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable && bEnableCoulomb && !bAxialSymm);

  pProp = m_wndPropList.FindItemByData(pObj->get_enable_quad_terms_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable && bEnableCoulomb && !bAxialSymm);

  bool bEnableRadialCoulomb = false;
  pProp = m_wndPropList.FindItemByData(pObj->get_use_radial_coulomb_ptr());
  if(pProp != NULL)
  {
    bool bEnableSelector = bEnable && bEnableCoulomb && !bAxialSymm;
    pProp->Enable(bEnableSelector);
    bEnableRadialCoulomb = bEnableSelector && pProp->GetValue().boolVal;
  }

  pProp = m_wndPropList.FindItemByData(pObj->get_radial_coulomb_trans_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableRadialCoulomb);

  bool bEnablePreCalc = false;
  pProp = m_wndPropList.FindItemByData(pObj->get_use_pre_calc_coulomb_ptr());
  if(pProp != NULL)
  {
    bool bEnableCheckBox = bEnable && bEnableCoulomb && !bAxialSymm;
    pProp->Enable(bEnableCheckBox);
    if(bEnableCheckBox)
      bEnablePreCalc = pProp->GetValue().boolVal;
  }

  pProp = m_wndPropList.FindItemByData(pObj->get_pre_calc_clmb_file_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable && bEnableCoulomb && !bAxialSymm && bEnablePreCalc);

  pProp = m_wndPropList.FindItemByData((DWORD_PTR)pObj);
  if(pProp != NULL)
    pProp->Enable(bEnable && bEnableCoulomb && !bAxialSymm);
}
