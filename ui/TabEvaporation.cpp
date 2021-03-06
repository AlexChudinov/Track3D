
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "EvaporationModel.h"
#include "ResponseProperty.h"
#include "DropletSizeGen.h"
#include "Button.h"

//---------------------------------------------------------------------------------------
// Evaporation
//---------------------------------------------------------------------------------------
void CPropertiesWnd::add_evapor_ctrls()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  CMFCPropertyGridProperty* pElectroGroup = new CMFCPropertyGridProperty(_T("Electrostatics"));
  CMFCPropertyGridProperty* pProp = new CMFCPropertyGridProperty(_T("Charge, elem. charges"), COleVariant(pObj->get_particle_charge() / Const_Charge_CGS), _T("Electric charge carried by a droplet."), pObj->get_particle_charge_ptr());
  pElectroGroup->AddSubItem(pProp);
  m_wndPropList.AddProperty(pElectroGroup);

// Initial diameter of droplets:
  CMFCPropertyGridProperty* pDiamGroup = new CMFCPropertyGridProperty(_T("Initial Droplet Diameter"));
  EvaporatingParticle::CDropletSizeGen& gen = pObj->get_size_gen_obj();
  int nSizeDistrType = gen.get_distr_type();
  CString sType(EvaporatingParticle::CDropletSizeGen::distr_name(nSizeDistrType));
  CGeneralResponseProperty* pDistrTypeProp = new CGeneralResponseProperty(this, _T("Distribution Type"), COleVariant(sType), _T("Select one of the following distribution types: Const, Homogeneous and Gaussian."), gen.get_distr_type_ptr());
  for(int i = EvaporatingParticle::CDropletSizeGen::dstConst; i < EvaporatingParticle::CDropletSizeGen::dstCount; i++)
    pDistrTypeProp->AddOption(EvaporatingParticle::CDropletSizeGen::distr_name(i));

  pDiamGroup->AddSubItem(pDistrTypeProp);

  CMFCPropertyGridProperty* pDistribGroup = new CMFCPropertyGridProperty(_T("Distribution Parameters"));

  switch(nSizeDistrType)
  {
    case EvaporatingParticle::CDropletSizeGen::dstConst:
    {
      pProp = new CMFCPropertyGridProperty(_T("Const. Droplet Diameter, mcm"), COleVariant(1.0e+4 * gen.get_mean_size()), _T("Mean droplets diameter. In the case of Constant Value distribution type this value will be initial diameter for all droplets."), gen.get_mean_size_ptr());
      pDistribGroup->AddSubItem(pProp);
      break;
    }
    case EvaporatingParticle::CDropletSizeGen::dstUniform:
    {
      pProp = new CMFCPropertyGridProperty(_T("Min Droplet Diameter, mcm"), COleVariant(1.0e+4 * gen.get_min_size()), _T("Minimal droplets diameter. Has effect for Homogeneous type of distribution only."), gen.get_min_size_ptr());
      pDistribGroup->AddSubItem(pProp);

      pProp = new CMFCPropertyGridProperty(_T("Max Droplet Diameter, mcm"), COleVariant(1.0e+4 * gen.get_max_size()), _T("Maximal droplets diameter. Has effect for Homogeneous type of distribution only."), gen.get_max_size_ptr());
      pDistribGroup->AddSubItem(pProp);
      break;
    }
    case EvaporatingParticle::CDropletSizeGen::dstGauss:
    {
      pProp = new CMFCPropertyGridProperty(_T("Mean Droplet Diameter, mcm"), COleVariant(1.0e+4 * gen.get_mean_size()), _T("Mean droplets diameter. In the case of Constant Value distribution type this value will be initial diameter for all droplets."), gen.get_mean_size_ptr());
      pDistribGroup->AddSubItem(pProp);

      pProp = new CMFCPropertyGridProperty(_T("Std. Dev. of Droplet Diameter, mcm"), COleVariant(1.0e+4 * gen.get_std_dev()), _T("Standard deviation of droplets diameter. Has effect for Gaussian type of distribution only."), gen.get_std_dev_ptr());
      pDistribGroup->AddSubItem(pProp);
    }
  }

  pDiamGroup->AddSubItem(pDistribGroup);
  m_wndPropList.AddProperty(pDiamGroup);

  CMFCPropertyGridProperty* pEvaporGroup = new CMFCPropertyGridProperty(_T("Evaporation"));

// Evaporation model:
  COleVariant var;
  switch(pObj->get_evapor_model_type())
  {
    case EvaporatingParticle::CTracker::emNone: var = CString("None"); break;
    case EvaporatingParticle::CTracker::emMaxwell: var = CString("Maxwell (isothermal)"); break;
    case EvaporatingParticle::CTracker::emSteadyDiffusive: var = CString("Diffusive (steady-state)"); break;
    case EvaporatingParticle::CTracker::emDiffusive: var = CString("Diffusive (transient)"); break;
  }

  CMFCPropertyGridProperty* pModGroup = new CMFCPropertyGridProperty(_T("Model"), var, _T("Model of particle's evaporation applied during tracking."), pObj->get_evapor_model_type_ptr());
  pModGroup->AddOption(_T("None"));
  pModGroup->AddOption(_T("Maxwell (isothermal)"));
  pModGroup->AddOption(_T("Diffusive (steady-state)"));
  pModGroup->AddOption(_T("Diffusive (transient)"));

  pEvaporGroup->AddSubItem(pModGroup);
/*
  CMFCPropertyGridProperty* pDiamGroup = new CMFCPropertyGridProperty(_T("Initial Diameter, mcm"));
  pProp = new CMFCPropertyGridProperty(_T("Initial Diameter"), COleVariant(1.e+4 * pObj->get_init_diameter()), _T("Diameter of a particle at the start point of the track."), pObj->get_init_diameter_ptr());
  pDiamGroup->AddSubItem(pProp);
  pEvaporGroup->AddSubItem(pDiamGroup);
*/
  EvaporatingParticle::CEvaporationModel* pMod = pObj->get_evapor_model();

  CMFCPropertyGridProperty* pHumidGroup = new CMFCPropertyGridProperty(_T("Environment Humidity, %"));
  pProp = new CMFCPropertyGridProperty(_T("Humidity"), COleVariant(100 * pMod->get_env_humidity()), _T("Environmental humidity"), pMod->get_env_humidity_ptr());
  pHumidGroup->AddSubItem(pProp);
  pEvaporGroup->AddSubItem(pHumidGroup);

  CMFCPropertyGridProperty* pTensGroup = new CMFCPropertyGridProperty(_T("Surface Tension"));
  CCheckBoxButton* pCheckBox = new CCheckBoxButton(this, _T("Enable"), (_variant_t)pMod->get_enable_surf_tens(), _T("Turns ON/OFF effect of additional pressure over a curved droplet's surface."), pMod->get_enable_surf_tens_ptr());
  pTensGroup->AddSubItem(pCheckBox);
  pEvaporGroup->AddSubItem(pTensGroup);

  switch(pMod->get_mass_trans_model())
  {
    case EvaporatingParticle::CEvaporationModel::mtmNone: var = CString("None"); break;
    case EvaporatingParticle::CEvaporationModel::mtmRanzMarshall: var = CString("Ranz-Marshall"); break;
  }

  pModGroup = new CMFCPropertyGridProperty(_T("Mass and Heat Transfer"));
  pProp = new CMFCPropertyGridProperty(_T("Model"), var, _T("Model of mass and heat transfer when a particle is moving through the environmental air. None implies that the relative movement is ignored."), pMod->get_mass_trans_model_ptr());
  pProp->AddOption(_T("None"));
  pProp->AddOption(_T("Ranz-Marshall"));
  pModGroup->AddSubItem(pProp);

  pEvaporGroup->AddSubItem(pModGroup);

  m_wndPropList.AddProperty(pEvaporGroup);
}

void CPropertiesWnd::set_evapor_data()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

// Electrostatics:
  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pObj->get_particle_charge_ptr());
  if(pProp != NULL)
    pObj->set_particle_charge(Const_Charge_CGS * pProp->GetValue().dblVal);

// Initial diameters distribution:
  EvaporatingParticle::CDropletSizeGen& gen = pObj->get_size_gen_obj();
  pProp = m_wndPropList.FindItemByData(gen.get_distr_type_ptr());
  if(pProp != NULL)
  {
    CString cDistrType = (CString)pProp->GetValue();
    for(int i = pProp->GetOptionCount() - 1; i >= 0; i--)
    {
      if((CString)pProp->GetOption(i) == cDistrType)
      {
        gen.set_distr_type(i);
        break;
      }
    }
  }

  pProp = m_wndPropList.FindItemByData(gen.get_mean_size_ptr());
  if(pProp != NULL)
    gen.set_mean_size(1.e-4 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(gen.get_std_dev_ptr());
  if(pProp != NULL)
    gen.set_std_dev(1.e-4 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(gen.get_min_size_ptr());
  if(pProp != NULL)
    gen.set_min_size(1.e-4 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(gen.get_max_size_ptr());
  if(pProp != NULL)
    gen.set_max_size(1.e-4 * pProp->GetValue().dblVal);

// Evaporation model:
  pProp = m_wndPropList.FindItemByData(pObj->get_evapor_model_type_ptr());
  if(pProp != NULL)
  {
    CString cMod = (CString)pProp->GetValue();
    for(int i = pProp->GetOptionCount() - 1; i >= 0; i--)
    {
      if((CString)pProp->GetOption(i) == cMod)
      {
        pObj->set_evapor_model_type(i);
        break;
      }
    }
  }

  EvaporatingParticle::CEvaporationModel* pMdl = pObj->get_evapor_model();

  pProp = m_wndPropList.FindItemByData(pMdl->get_env_humidity_ptr());
  if(pProp != NULL)
    pMdl->set_env_humidity(0.01 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pMdl->get_mass_trans_model_ptr());
  if(pProp != NULL)
  {
    CString cMassTransMod = (CString)pProp->GetValue();
    for(int i = pProp->GetOptionCount() - 1; i >= 0; i--)
    {
      if((CString)pProp->GetOption(i) == cMassTransMod)
      {
        pMdl->set_mass_trans_model(i);
        break;
      }
    }
  }
}

void CPropertiesWnd::update_evapor_ctrls()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  BOOL bEnable = pObj->get_particle_type() == EvaporatingParticle::CTrack::ptDroplet;

// Electrostatics:
  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pObj->get_particle_charge_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

// Initial diameters distribution:
  EvaporatingParticle::CDropletSizeGen& gen = pObj->get_size_gen_obj();
  pProp = m_wndPropList.FindItemByData(gen.get_distr_type_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(gen.get_mean_size_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(gen.get_std_dev_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(gen.get_min_size_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(gen.get_max_size_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

// Evaporation model:
  EvaporatingParticle::CEvaporationModel* pMdl = pObj->get_evapor_model();

  pProp = m_wndPropList.FindItemByData(pMdl->get_enable_surf_tens_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  BOOL bEnableEvapor = FALSE;
  pProp = m_wndPropList.FindItemByData(pObj->get_evapor_model_type_ptr());
  if(pProp != NULL)
  {
    pProp->Enable(bEnable);
    CString cMod = (CString)pProp->GetValue();
    bEnableEvapor = bEnable && (BOOL)(cMod != "None");
  }

  pProp = m_wndPropList.FindItemByData(pMdl->get_env_humidity_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableEvapor);

  pProp = m_wndPropList.FindItemByData(pMdl->get_mass_trans_model_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnableEvapor);
}
