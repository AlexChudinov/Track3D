
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"

//---------------------------------------------------------------------------------------
// Tracking controls
//---------------------------------------------------------------------------------------
void CPropertiesWnd::add_tracking_ctrls()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

// Multi-threading support:
  CMFCPropertyGridProperty* pMultiGroup = new CMFCPropertyGridProperty(_T("Multi-Threading"));
  CMFCPropertyGridProperty* pProp = new CMFCPropertyGridProperty(_T("Enable"), (_variant_t)pObj->get_use_multi_thread(), _T("Turns ON/OFF multi-threrading support."), pObj->get_use_multi_thread_ptr());
  pMultiGroup->AddSubItem(pProp);

  m_wndPropList.AddProperty(pMultiGroup);

// Integrator type:
  CMFCPropertyGridProperty* pIntegrGroup = new CMFCPropertyGridProperty(_T("Integrators"));
  COleVariant var(pObj->get_integr_name(pObj->get_integr_type()));
  pProp = new CMFCPropertyGridProperty(_T("Integrator Type:"), var, _T("Select a method of particle's parameters integration along the trajectory."), pObj->get_integr_type_ptr());
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
  pProp = new CMFCPropertyGridProperty(_T("Enable DC Field"), (_variant_t)pObj->get_enable_field(), _T("Turns ON/OFF application of DC electric field in particles tracking."), pObj->get_enable_field_ptr());
  pFieldGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("DC Voltage, V"), COleVariant(pObj->get_dc_amplitude()), _T("DC potential applied to the emitter"), pObj->get_dc_amplitude_ptr());
  pFieldGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Charge, elem. charges"), COleVariant(long(pObj->get_particle_charge() / Const_Charge_CGS)), _T("Electric charge carried by a particle."), pObj->get_particle_charge_ptr());
  pFieldGroup->AddSubItem(pProp);

  m_wndPropList.AddProperty(pFieldGroup);

// Radio-Frequency field:
  CMFCPropertyGridProperty* pRFGroup = new CMFCPropertyGridProperty(_T("Radio-Frequency Field"));
  pProp = new CMFCPropertyGridProperty(_T("Enable RF"), (_variant_t)pObj->get_enable_rf(), _T("Turns ON/OFF radio-frequency field."), pObj->get_enable_rf_ptr());
  pRFGroup->AddSubItem(pProp);
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
  pProp = new CMFCPropertyGridProperty(_T("Frequency, kHz"), COleVariant(0.001 * pObj->get_rf_flatapole_freq()), _T("Frquency of RF voltage applied electrodes in flatapole region."), pObj->get_rf_flatapole_freq_ptr());
  pRFFlatapoleGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Transition to Flatapole, mm"), COleVariant(10 * pObj->get_flatapole_trans()), _T("X-coordinate of transition to flatapole region."), pObj->get_flatapole_trans_ptr());
  pRFFlatapoleGroup->AddSubItem(pProp);
  pRFGroup->AddSubItem(pRFFlatapoleGroup);

  m_wndPropList.AddProperty(pRFGroup);

// Old Predictor-Corrector switcher:
  CMFCPropertyGridProperty* pOldIntegratorGroup = new CMFCPropertyGridProperty(_T("Use the Old Integrator"));
  pProp = new CMFCPropertyGridProperty(_T("Enable Old Integrator"), (_variant_t)pObj->get_use_old_integrator(), _T("If this is true the old integrator is used. Note that the droplet type of particles is supported by the old integrator only"), pObj->get_use_old_integrator_ptr());
  pOldIntegratorGroup->AddSubItem(pProp);

  m_wndPropList.AddProperty(pOldIntegratorGroup);
}

void CPropertiesWnd::set_tracking_data()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

// Multi-threading support:
  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pObj->get_use_multi_thread_ptr());
  if(pProp != NULL)
    pObj->set_use_multi_thread(pProp->GetValue().boolVal);

// Integrator type:
  pProp = m_wndPropList.FindItemByData(pObj->get_integr_type_ptr());
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
  pProp = m_wndPropList.FindItemByData(pObj->get_enable_field_ptr());
  if(pProp != NULL)
    pObj->set_enable_field(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_particle_charge_ptr());
  if(pProp != NULL)
    pObj->set_particle_charge(Const_Charge_CGS * pProp->GetValue().lVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_dc_amplitude_ptr());
  if(pProp != NULL)
    pObj->set_dc_amplitude(pProp->GetValue().dblVal);

// Radio-Frequency field:
    pProp = m_wndPropList.FindItemByData(pObj->get_enable_rf_ptr());
  if(pProp != NULL)
    pObj->set_enable_rf(pProp->GetValue().boolVal);

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

// Old Integrator switcher:
  pProp = m_wndPropList.FindItemByData(pObj->get_use_old_integrator_ptr());
  if(pProp != NULL)
    pObj->set_use_old_integrator(pProp->GetValue().boolVal);
}

void CPropertiesWnd::update_tracking_ctrls()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  BOOL bEnable = FALSE;
  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pObj->get_enable_field_ptr());
  if(pProp != NULL)
    bEnable = (BOOL)(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_dc_amplitude_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

// Radio-Frequency field:
    bool bEnableRF = FALSE;
  pProp = m_wndPropList.FindItemByData(pObj->get_enable_rf_ptr());
  if(pProp != NULL)
  {
    pProp->Enable(TRUE);
    bEnableRF = (BOOL)(pProp->GetValue().boolVal);
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

// Old Integrator switcher:
  BOOL bOldIntegr = FALSE;
  pProp = m_wndPropList.FindItemByData(pObj->get_use_old_integrator_ptr());
  if(pProp != NULL)
    bOldIntegr = (BOOL)(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(pObj->get_integr_type_ptr());
  if(pProp != NULL)
    pProp->Enable(!bOldIntegr);
}
