
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "Button.h"

//---------------------------------------------------------------------------------------
// Tracking controls
//---------------------------------------------------------------------------------------
void CPropertiesWnd::add_tracking_ctrls()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

// Multi-threading support:
  CMFCPropertyGridProperty* pMultiGroup = new CMFCPropertyGridProperty(_T("Multi-Threading"));
  CCheckBoxButton* pCheckBox = new CCheckBoxButton(this, _T("Enable"), (_variant_t)pObj->get_use_multi_thread(), _T("Turns ON/OFF multi-threrading support."), pObj->get_use_multi_thread_ptr());
  pMultiGroup->AddSubItem(pCheckBox);

  m_wndPropList.AddProperty(pMultiGroup);

// Integrator type:
  CMFCPropertyGridProperty* pIntegrGroup = new CMFCPropertyGridProperty(_T("Integrators"));
  COleVariant var(pObj->get_integr_name(pObj->get_integr_type()));
  CMFCPropertyGridProperty* pProp = new CMFCPropertyGridProperty(_T("Integrator Type:"), var, _T("Select a method of particle's parameters integration along the trajectory."), pObj->get_integr_type_ptr());
  for(int i = 0; i < EvaporatingParticle::CTracker::intCount; i++)
    pProp->AddOption(pObj->get_integr_name(i));

  pIntegrGroup->AddSubItem(pProp);
  m_wndPropList.AddProperty(pIntegrGroup);

// Time parameters:
  CMFCPropertyGridProperty* pTimeGroup = new CMFCPropertyGridProperty(_T("Time Parameters"));
  pProp = new CMFCPropertyGridProperty(_T("Time Step, mcs"), COleVariant(1.e+6 * pObj->get_time_step()), _T("Integration time step."), pObj->get_time_step_ptr());
  pTimeGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Maximal Tracking Time, s"), COleVariant(pObj->get_max_track_time()), _T("Maximal allowed tracking time."), pObj->get_max_track_time_ptr());
  pTimeGroup->AddSubItem(pProp);

  EvaporatingParticle::COutputEngine& outEng = pObj->get_output_engine();
  pProp = new CMFCPropertyGridProperty(_T("Output Time Step, mcs"), COleVariant(1.e+6 * outEng.get_output_time_step()), _T("Time interval between consecutive points in output data."), outEng.get_output_time_step_ptr());
  pTimeGroup->AddSubItem(pProp);

  m_wndPropList.AddProperty(pTimeGroup);

// Electrostatics:
  CMFCPropertyGridProperty* pFieldGroup = new CMFCPropertyGridProperty(_T("Electrostatics"));
  pCheckBox = new CCheckBoxButton(this, _T("Enable DC Field"), (_variant_t)pObj->get_enable_field(), _T("Turns ON/OFF application of DC electric field in particles tracking."), pObj->get_enable_field_ptr());
  pFieldGroup->AddSubItem(pCheckBox);
  pProp = new CMFCPropertyGridProperty(_T("DC Voltage, V"), COleVariant(pObj->get_dc_amplitude()), _T("DC potential applied to the emitter"), pObj->get_dc_amplitude_ptr());
  pFieldGroup->AddSubItem(pProp);

  m_wndPropList.AddProperty(pFieldGroup);

// Radio-Frequency field:
  CMFCPropertyGridProperty* pRFGroup = new CMFCPropertyGridProperty(_T("Radio-Frequency Field"));
  pCheckBox = new CCheckBoxButton(this, _T("Enable RF"), (_variant_t)pObj->get_enable_rf(), _T("Turns ON/OFF radio-frequency field."), pObj->get_enable_rf_ptr());
  pRFGroup->AddSubItem(pCheckBox);
  pProp = new CMFCPropertyGridProperty(_T("Amplitude, V"), COleVariant(pObj->get_rf_amplitude()), _T("Amplitude of RF voltage applied to S-Lens electrodes."), pObj->get_rf_amplitude_ptr());
  pRFGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Frequency, kHz"), COleVariant(0.001 * pObj->get_rf_frequency()), _T("Frquency of RF voltage applied to S-Lens electrodes."), pObj->get_rf_frequency_ptr());
  pRFGroup->AddSubItem(pProp);

// RF in Q00:
  CMFCPropertyGridProperty* pRFQ00Group = new CMFCPropertyGridProperty(_T("RF in Q00 Region"));
  pProp = new CMFCPropertyGridProperty(_T("Amplitude, V"), COleVariant(pObj->get_rf_Q00_ampl()), _T("Amplitude of RF voltage applied to electrodes in Q00 region."), pObj->get_rf_Q00_ampl_ptr());
  pRFQ00Group->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Frequency, kHz"), COleVariant(0.001 * pObj->get_rf_Q00_freq()), _T("Frquency of RF voltage applied electrodes in Q00 region."), pObj->get_rf_Q00_freq_ptr());
  pRFQ00Group->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Transition to Q00, mm"), COleVariant(10 * pObj->get_Q00_trans()), _T("X-coordinate of transition to Q00 region."), pObj->get_Q00_trans_ptr());
  pRFQ00Group->AddSubItem(pProp);
  pRFGroup->AddSubItem(pRFQ00Group);

// RF in flatapole:
  CMFCPropertyGridProperty* pRFFlatapoleGroup = new CMFCPropertyGridProperty(_T("RF in Flatapole"));
  pProp = new CMFCPropertyGridProperty(_T("Amplitude, V"), COleVariant(pObj->get_rf_flatapole_ampl()), _T("Amplitude of RF voltage applied to electrodes in flatapole region."), pObj->get_rf_flatapole_ampl_ptr());
  pRFFlatapoleGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Frequency, kHz"), COleVariant(0.001 * pObj->get_rf_flatapole_freq()), _T("Frquency of RF voltage applied to electrodes in flatapole region."), pObj->get_rf_flatapole_freq_ptr());
  pRFFlatapoleGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Transition to Flatapole, mm"), COleVariant(10 * pObj->get_flatapole_trans()), _T("X-coordinate of transition to flatapole region."), pObj->get_flatapole_trans_ptr());
  pRFFlatapoleGroup->AddSubItem(pProp);
  pRFGroup->AddSubItem(pRFFlatapoleGroup);

  m_wndPropList.AddProperty(pRFGroup);

// Random diffusion group:
  CMFCPropertyGridProperty* pDiffusionGroup = new CMFCPropertyGridProperty(_T("Random Processes"));
  pCheckBox = new CCheckBoxButton(this, _T("Enable Diffusion"), (_variant_t)pObj->get_enable_diffusion(), _T("If this is set to 'true' the ion positions (or ion velocities) are disturbed by random variations at every time step. The random diffusion is applied if at X < Xc"), pObj->get_enable_diffusion_ptr());
  pDiffusionGroup->AddSubItem(pCheckBox);

  COleVariant var1(RandomProcess::rndProcName(pObj->get_rand_diff_type()));
  pProp = new CMFCPropertyGridProperty(_T("Random Diffusion Type"), var1, _T("Select the type of random diffusion model."), pObj->get_rand_diff_type_ptr());
  pProp->AddOption(RandomProcess::rndProcName(RandomProcess::DIFFUSION_VELOCITY_JUMP));
  pProp->AddOption(RandomProcess::rndProcName(RandomProcess::DIFFUSION_COORD_JUMP));
  pProp->AllowEdit(FALSE);
  pDiffusionGroup->AddSubItem(pProp);

  pCheckBox = new CCheckBoxButton(this, _T("Enable Collisions"), (_variant_t)pObj->get_enable_collisions(), _T("If this is set to 'true' the ion positions are disturbed by random variations once per several time steps. The random collisions are applied if at X > Xc"), pObj->get_enable_collisions_ptr());
  pDiffusionGroup->AddSubItem(pCheckBox);

  COleVariant var2(RandomProcess::rndProcName(pObj->get_rand_collision_type()));
  pProp = new CMFCPropertyGridProperty(_T("Random Collisions Model"), var2, _T("Select the type of random collisions model."), pObj->get_rand_collision_type_ptr());
  pProp->AddOption(RandomProcess::rndProcName(RandomProcess::COLLISION));
  pProp->AddOption(RandomProcess::rndProcName(RandomProcess::COLLISION_ANY_PRESS));
  pProp->AllowEdit(FALSE);
  pDiffusionGroup->AddSubItem(pProp);

  pProp = new CMFCPropertyGridProperty(_T("Limiting X-Coordinate, mm"), COleVariant(10 * pObj->get_diffusion_switch_cond()), _T("Xc the x-coordinate, limiting the application of both random diffusion and random collision models."), pObj->get_diffusion_switch_cond_ptr());
  pDiffusionGroup->AddSubItem(pProp);

  pProp = new CMFCPropertyGridProperty(_T("Random Seed"), COleVariant(pObj->get_random_seed()), _T("Seed for the random numbers generator."), pObj->get_random_seed_ptr());
  pDiffusionGroup->AddSubItem(pProp);

  m_wndPropList.AddProperty(pDiffusionGroup);

// Backward compatibility:
  CMFCPropertyGridProperty* pOldIntegratorGroup = new CMFCPropertyGridProperty(_T("Backward compatibility"));

  pCheckBox = new CCheckBoxButton(this, _T("Use ANSYS Electric Fields"), (_variant_t)pObj->get_enable_ansys_field(), _T("If this is true the electric fields computed in ANSYS are used."), pObj->get_enable_ansys_field_ptr());
  pOldIntegratorGroup->AddSubItem(pCheckBox);
  pCheckBox = new CCheckBoxButton(this, _T("Save Screen Image"), (_variant_t)pObj->get_save_image(), _T("If this is true the screen image will be captured and saved at every iteration."), pObj->get_save_image_ptr());
  pOldIntegratorGroup->AddSubItem(pCheckBox);

  m_wndPropList.AddProperty(pOldIntegratorGroup);
}

void CPropertiesWnd::set_tracking_data()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

// Integrator type:
  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pObj->get_integr_type_ptr());
  if(pProp != NULL)
  {
    CString cType = (CString)pProp->GetValue();
    for(int i = 0; i < EvaporatingParticle::CTracker::intCount; i++)
    {
      if(pObj->get_integr_name(i) == cType)
      {
        pObj->set_integr_type(i);
        break;
      }
    }
  }

// Time parameters:
  pProp = m_wndPropList.FindItemByData(pObj->get_time_step_ptr());
  if(pProp != NULL)
    pObj->set_time_step(1.e-6 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_max_track_time_ptr());
  if(pProp != NULL)
    pObj->set_max_track_time(pProp->GetValue().dblVal);

  EvaporatingParticle::COutputEngine& outEng = pObj->get_output_engine();
  pProp = m_wndPropList.FindItemByData(outEng.get_output_time_step_ptr());
  if(pProp != NULL)
   outEng.set_output_time_step(1.e-6 * pProp->GetValue().dblVal);

// Electrostatics:
  pProp = m_wndPropList.FindItemByData(pObj->get_dc_amplitude_ptr());
  if(pProp != NULL)
    pObj->set_dc_amplitude(pProp->GetValue().dblVal);

// Radio-Frequency field:
  pProp = m_wndPropList.FindItemByData(pObj->get_rf_amplitude_ptr());
  if(pProp != NULL)
    pObj->set_rf_amplitude(pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_rf_frequency_ptr());
  if(pProp != NULL)
    pObj->set_rf_frequency(1000 * pProp->GetValue().dblVal);

// RF in Q00:
  pProp = m_wndPropList.FindItemByData(pObj->get_rf_Q00_ampl_ptr());
  if(pProp != NULL)
    pObj->set_rf_Q00_ampl(pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_rf_Q00_freq_ptr());
  if(pProp != NULL)
    pObj->set_rf_Q00_freq(1000 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_Q00_trans_ptr());
  if(pProp != NULL)
    pObj->set_Q00_trans(0.1 * pProp->GetValue().dblVal);

// RF in flatapole:
  pProp = m_wndPropList.FindItemByData(pObj->get_rf_flatapole_ampl_ptr());
  if(pProp != NULL)
    pObj->set_rf_flatapole_ampl(pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_rf_flatapole_freq_ptr());
  if(pProp != NULL)
    pObj->set_rf_flatapole_freq(1000 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_flatapole_trans_ptr());
  if(pProp != NULL)
    pObj->set_flatapole_trans(0.1 * pProp->GetValue().dblVal);

// Random processes group:
  pProp = m_wndPropList.FindItemByData(pObj->get_diffusion_switch_cond_ptr());
  if(pProp != NULL)
    pObj->set_diffusion_switch_cond(0.1 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_random_seed_ptr());
  if(pProp != NULL)
    pObj->set_random_seed(pProp->GetValue().lVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_rand_diff_type_ptr());
  if(pProp != NULL)
  {
    CString cTypeName = (CString)pProp->GetValue();
    if(cTypeName == RandomProcess::rndProcName(RandomProcess::DIFFUSION_VELOCITY_JUMP))
      pObj->set_rand_diff_type(RandomProcess::DIFFUSION_VELOCITY_JUMP);
    else
      pObj->set_rand_diff_type(RandomProcess::DIFFUSION_COORD_JUMP);
  }

  pProp = m_wndPropList.FindItemByData(pObj->get_rand_collision_type_ptr());
  if(pProp != NULL)
  {
    CString cTypeName = (CString)pProp->GetValue();
    if(cTypeName == RandomProcess::rndProcName(RandomProcess::COLLISION))
      pObj->set_rand_collision_type(RandomProcess::COLLISION);
    else
      pObj->set_rand_collision_type(RandomProcess::COLLISION_ANY_PRESS);
  }
}

void CPropertiesWnd::update_tracking_ctrls()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  bool bUseAnsysFields = pObj->get_enable_ansys_field();

  BOOL bEnable = FALSE;
  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pObj->get_enable_field_ptr());
  if(pProp != NULL)
  {
    bEnable = bUseAnsysFields && (BOOL)(pProp->GetValue().boolVal);
    pProp->Enable(bUseAnsysFields);
  }

  pProp = m_wndPropList.FindItemByData(pObj->get_dc_amplitude_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

// Radio-Frequency field:
  bool bEnableRF = FALSE;
  pProp = m_wndPropList.FindItemByData(pObj->get_enable_rf_ptr());
  if(pProp != NULL)
  {
    bEnableRF = bUseAnsysFields && (BOOL)(pProp->GetValue().boolVal);
    pProp->Enable(bUseAnsysFields);
  }

  pProp = m_wndPropList.FindItemByData(pObj->get_rf_amplitude_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableRF);

  pProp = m_wndPropList.FindItemByData(pObj->get_rf_frequency_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableRF);

// RF in Q00:
  pProp = m_wndPropList.FindItemByData(pObj->get_rf_Q00_ampl_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableRF);

  pProp = m_wndPropList.FindItemByData(pObj->get_rf_Q00_freq_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableRF);

  pProp = m_wndPropList.FindItemByData(pObj->get_Q00_trans_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableRF);

// RF in flatapole:
  pProp = m_wndPropList.FindItemByData(pObj->get_rf_flatapole_ampl_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableRF);

  pProp = m_wndPropList.FindItemByData(pObj->get_rf_flatapole_freq_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableRF);

  pProp = m_wndPropList.FindItemByData(pObj->get_flatapole_trans_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableRF);

// Random processes:
  bool bDiffOn = pObj->get_enable_diffusion();
  bool bCollOn = pObj->get_enable_collisions();
  pProp = m_wndPropList.FindItemByData(pObj->get_enable_diffusion_ptr());
  if(pProp != NULL)
    bDiffOn = pProp->GetValue().boolVal;

  pProp = m_wndPropList.FindItemByData(pObj->get_enable_collisions_ptr());
  if(pProp != NULL)
    bCollOn = pProp->GetValue().boolVal;

  pProp = m_wndPropList.FindItemByData(pObj->get_diffusion_switch_cond_ptr());
  if(pProp != NULL)
    pProp->Enable(bDiffOn || bCollOn);

  pProp = m_wndPropList.FindItemByData(pObj->get_random_seed_ptr());
  if(pProp != NULL)
    pProp->Enable(bDiffOn || bCollOn);

  pProp = m_wndPropList.FindItemByData(pObj->get_rand_diff_type_ptr());
  if(pProp != NULL)
    pProp->Enable(bDiffOn);

  pProp = m_wndPropList.FindItemByData(pObj->get_rand_collision_type_ptr());
  if(pProp != NULL)
    pProp->Enable(bCollOn);
}
