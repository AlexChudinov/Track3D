
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "ResponseProperty.h"
#include "Button.h"

using namespace EvaporatingParticle;
//---------------------------------------------------------------------------------------
// Source of particles
//---------------------------------------------------------------------------------------
void CPropertiesWnd::add_source_ctrls()
{
  EvaporatingParticle::CSource* pSrc = CParticleTrackingApp::Get()->GetTracker()->get_src();
  CMFCPropertyGridProperty* pSrcGroup = new CMFCPropertyGridProperty(_T("Particle Source"));

// Type of the source:
  int nSrcType = pSrc->get_src_type();
  COleVariant var(pSrc->get_src_type_name(nSrcType));
  CGeneralResponseProperty* pType = new CGeneralResponseProperty(this, _T("Source Type"), var, _T("Specify the type of particle's source: Cone, Round Spot, Cylinder (for 2D), Rectangular Spot, Selected Regions and Positions from File."), pSrc->get_src_type_ptr());
  for(int i = 0; i < EvaporatingParticle::CSource::stCount; i++)
    pType->AddOption(pSrc->get_src_type_name(i));

  pType->AllowEdit(FALSE);
  pSrcGroup->AddSubItem(pType);

  int nSrcDistrType = pSrc->get_src_inject_type();

  bool bFromFile = (nSrcType == CSource::stPosFromFile);
  if(!bFromFile)
  {
    COleVariant vur(pSrc->get_inj_type_name(nSrcDistrType));
    CGeneralResponseProperty* pInjType = new CGeneralResponseProperty(this, _T("Distribution Type"), vur, _T("Distribution can be of two types: random or homogenious."), pSrc->get_src_inject_type_ptr());
    for(int j = 0; j < EvaporatingParticle::CSource::itCount; j++)
      pInjType->AddOption(pSrc->get_inj_type_name(j));

    pInjType->AllowEdit(FALSE);
    pSrcGroup->AddSubItem(pInjType);
  }

  CMFCPropertyGridProperty* pProp = NULL;
  

  bool bRandom = nSrcDistrType == EvaporatingParticle::CSource::itRandom;
  if(bRandom && !bFromFile)
  {
    pProp = new CMFCPropertyGridProperty(_T("Random Seed"), COleVariant((long)pSrc->get_random_seed()), _T("Initial number in a pseudo-random sequence of numbers. Two sequences are equal to each other if their random seeds are equal."), pSrc->get_random_seed_ptr());
    pSrcGroup->AddSubItem(pProp);
  }

// Position of the source:
  bool bNeedPos = (nSrcType != CSource::stSelReg) && (!bFromFile);
  if(bNeedPos)
  {
    CMFCPropertyGridProperty* pPosGroup = new CMFCPropertyGridProperty(_T("Position, mm"), pSrc->get_src_pos_ptr());
    pProp = new CMFCPropertyGridProperty(_T("X"), COleVariant(10 * pSrc->get_src_pos().x), _T("X coordinate of particle source."));
    pPosGroup->AddSubItem(pProp);
    pProp = new CMFCPropertyGridProperty(_T("Y"), COleVariant(10 * pSrc->get_src_pos().y), _T("Y coordinate of particle source."));
    pPosGroup->AddSubItem(pProp);
    pProp = new CMFCPropertyGridProperty(_T("Z"), COleVariant(10 * pSrc->get_src_pos().z), _T("Z coordinate of particle source."));
    pPosGroup->AddSubItem(pProp);
    pSrcGroup->AddSubItem(pPosGroup);
  }

  CMFCPropertyGridProperty* pDirGroup = new CMFCPropertyGridProperty(_T("Normal Direction"), pSrc->get_inject_dir_ptr());
  pProp = new CMFCPropertyGridProperty(_T("X"), COleVariant(pSrc->get_inject_dir().x), _T("X component of normal vector to the Spot or Rectangle plane."));
  pDirGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Y"), COleVariant(pSrc->get_inject_dir().y), _T("Y component of normal vector to the Spot or Rectangle plane."));
  pDirGroup->AddSubItem(pProp);
  pProp = new CMFCPropertyGridProperty(_T("Z"), COleVariant(pSrc->get_inject_dir().z), _T("Z component of normal vector to the Spot or Rectangle plane."));
  pDirGroup->AddSubItem(pProp);
  pSrcGroup->AddSubItem(pDirGroup);

// Cone parameters:
  bool bNeedCone = (nSrcType == CSource::stCone) || (nSrcType == CSource::stCylinder);
  if(bNeedCone)
  {
    CMFCPropertyGridProperty* pConeGroup = new CMFCPropertyGridProperty(_T("Cone and Cylinder"));
    pProp = new CMFCPropertyGridProperty(_T("Cone angle, deg"), COleVariant(pSrc->get_cone_angle()), _T("Cone aperture angle."), pSrc->get_cone_angle_ptr());
    pConeGroup->AddSubItem(pProp);
    pSrcGroup->AddSubItem(pConeGroup);
  }

// Spot and sphere parameters:
  bool bNeedRadius = (nSrcType == CSource::stSpot) || (nSrcType == CSource::stCylinder);
  if(bNeedRadius)
  {
    CMFCPropertyGridProperty* pSpotGroup = new CMFCPropertyGridProperty(_T("Spot and Cylinder"));
    pProp = new CMFCPropertyGridProperty(_T("Radius, mm"), COleVariant(10 * pSrc->get_radius()), _T("Radius of the circular spot."), pSrc->get_radius_ptr());
    pSpotGroup->AddSubItem(pProp);
    pSrcGroup->AddSubItem(pSpotGroup);
  }

// Rectangular spot:
  bool bRectSpot = (nSrcType == CSource::stRect);
  if(bRectSpot)
  {
    CMFCPropertyGridProperty* pRectGroup = new CMFCPropertyGridProperty(_T("Rectangular Spot"));
    pProp = new CMFCPropertyGridProperty(_T("Width, mm"), COleVariant(10 * pSrc->get_rect_width()), _T("Width of the rectangular spot."), pSrc->get_rect_width_ptr());
    pRectGroup->AddSubItem(pProp);
    pProp = new CMFCPropertyGridProperty(_T("Height, mm"), COleVariant(10 * pSrc->get_rect_height()), _T("Height of the rectangular spot or of the cylinder."), pSrc->get_rect_height_ptr());
    pRectGroup->AddSubItem(pProp);
    pSrcGroup->AddSubItem(pRectGroup);
  }

// Selected region source:
  bool bSelRegSrc = (nSrcType == CSource::stSelReg);
  if(bSelRegSrc)
  {
    CMFCPropertyGridProperty* pSelRegGroup = new CMFCPropertyGridProperty(_T("Start from Selected 2D-Region"));
    CString cSelRegNames = pSrc->get_selected_rgn_names();
    CSelectRegionButton* pSelRegButton = new CSelectRegionButton(this, _T("2D Regions"), cSelRegNames, _T("Click to select 2D regions from which the particles will start."), pSrc->get_selected_rgn_names_ptr());
    pSelRegGroup->AddSubItem(pSelRegButton);
    pSrcGroup->AddSubItem(pSelRegGroup);
  }

  CMFCPropertyGridProperty* pVelGroup = new CMFCPropertyGridProperty(_T("Initial Velocity, m/s"));
  CSourceCheckBox* pCheckBox = new CSourceCheckBox(this, _T("Use Gas Velocity"), (_variant_t)pSrc->get_use_initial_gas_vel(), _T("If ON the initial velocity is taken from the gas-dynamic data. Otherwise, the velocity direction is defined by normal to Spot and Rectangle planes."), pSrc->get_use_initial_gas_vel_ptr());
  pVelGroup->AddSubItem(pCheckBox);
  pProp = new CMFCPropertyGridProperty(_T("Abs. Velocity"), COleVariant(0.01 * pSrc->get_abs_vel()), _T("Absolute value of the starting velocity of a particle."), pSrc->get_abs_vel_ptr());
  pVelGroup->AddSubItem(pProp);
  pSrcGroup->AddSubItem(pVelGroup);

  CMFCPropertyGridProperty* pCountGroup = new CMFCPropertyGridProperty(_T("Count of Particles"));
  if(nSrcType != CSource::stPosFromFile)
  {
    pProp = new CMFCPropertyGridProperty(_T("Count of Particles"), COleVariant((long)pSrc->get_particles_count()), _T("Total number of particles."), pSrc->get_particles_count_ptr());
    pCountGroup->AddSubItem(pProp);
  }

  bool bPitch = ((nSrcType == CSource::stCone) || (nSrcType == CSource::stCylinder)) && (nSrcDistrType == CSource::itHomogen);
  if(bPitch)
  {
    pProp = new CMFCPropertyGridProperty(_T("Pitch Number"), COleVariant((long)pSrc->get_pitch_count()), _T("Number of pitch subdivisions."), pSrc->get_pitch_count_ptr());
    pCountGroup->AddSubItem(pProp);
    pProp = new CMFCPropertyGridProperty(_T("Azimuth Number"), COleVariant((long)pSrc->get_azim_count()), _T("Number of azimuthal subdivisions."), pSrc->get_azim_count_ptr());
    pCountGroup->AddSubItem(pProp);
  }

  pProp = new CMFCPropertyGridProperty(_T("Ensemble Size"), COleVariant((long)pSrc->get_ensemble_size()), _T("Count of equvalent ions in the ensemble."), pSrc->get_ensemble_size_ptr());
  pCountGroup->AddSubItem(pProp);
  pSrcGroup->AddSubItem(pCountGroup);

  if(bFromFile)
  {
    CMFCPropertyGridProperty* pFileGroup = new CMFCPropertyGridProperty(_T("Input File"));
    static TCHAR BASED_CODE szFilter[] = _T("CSV Files(*.csv)|*.csv|All Files(*.*)|*.*||");
    CMFCPropertyGridFileProperty* pFileProp = new CMFCPropertyGridFileProperty(_T("File Name"), TRUE, pSrc->get_filename(), _T("csv"),
      OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, _T("Specify location of file containing initial positions of particles."), pSrc->get_filename_ptr());

    pFileGroup->AddSubItem(pFileProp);
    pSrcGroup->AddSubItem(pFileGroup);
  }

  m_wndPropList.AddProperty(pSrcGroup);
}

void CPropertiesWnd::set_source_data()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  EvaporatingParticle::CSource* pSrc = pObj->get_src();

  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pSrc->get_src_type_ptr());
  if(pProp != NULL)
  {
    CString cType = (CString)pProp->GetValue();
    for(int i = 0; i < EvaporatingParticle::CSource::stCount; i++)
    {
      if(cType == pSrc->get_src_type_name(i))
      {
        pSrc->set_src_type(i);
        break;
      }
    }
  }

  pProp = m_wndPropList.FindItemByData(pSrc->get_src_inject_type_ptr());
  if(pProp != NULL)
  {
    CString cType = (CString)pProp->GetValue();
    for(int j = 0; j < EvaporatingParticle::CSource::itCount; j++)
    {
      if(cType == pSrc->get_inj_type_name(j))
      {
        pSrc->set_src_inject_type(j);
        break;
      }
    }
  }

  pProp = m_wndPropList.FindItemByData(pSrc->get_random_seed_ptr());
  if(pProp != NULL)
    pSrc->set_random_seed(pProp->GetValue().lVal);

  pProp = m_wndPropList.FindItemByData(pSrc->get_use_initial_gas_vel_ptr());
  if(pProp != NULL)
    pSrc->set_use_initial_gas_vel(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(pSrc->get_abs_vel_ptr());
  if(pProp != NULL)
    pSrc->set_abs_vel(100 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pSrc->get_cone_angle_ptr());
  if(pProp != NULL)
    pSrc->set_cone_angle(pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pSrc->get_radius_ptr());
  if(pProp != NULL)
    pSrc->set_radius(0.1 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pSrc->get_rect_width_ptr());
  if(pProp != NULL)
    pSrc->set_rect_width(0.1 * pProp->GetValue().dblVal);

  pProp = m_wndPropList.FindItemByData(pSrc->get_rect_height_ptr());
  if(pProp != NULL)
    pSrc->set_rect_height(0.1 * pProp->GetValue().dblVal);

  EvaporatingParticle::Vector3D vPos;
  pProp = m_wndPropList.FindItemByData(pSrc->get_src_pos_ptr());
  if(pProp != NULL)
  {
    vPos.x = 0.1 * pProp->GetSubItem(0)->GetValue().dblVal;
    vPos.y = 0.1 * pProp->GetSubItem(1)->GetValue().dblVal;
    vPos.z = 0.1 * pProp->GetSubItem(2)->GetValue().dblVal;
    pSrc->set_src_pos(vPos);
  }

  EvaporatingParticle::Vector3D vDir;
  pProp = m_wndPropList.FindItemByData(pSrc->get_inject_dir_ptr());
  if(pProp != NULL)
  {
    vDir.x = pProp->GetSubItem(0)->GetValue().dblVal;
    vDir.y = pProp->GetSubItem(1)->GetValue().dblVal;
    vDir.z = pProp->GetSubItem(2)->GetValue().dblVal;
    pSrc->set_inject_dir(vDir);
  }

  pProp = m_wndPropList.FindItemByData(pSrc->get_particles_count_ptr());
  if(pProp != NULL)
    pSrc->set_particles_count(pProp->GetValue().lVal);

  pProp = m_wndPropList.FindItemByData(pSrc->get_pitch_count_ptr());
  if(pProp != NULL)
    pSrc->set_pitch_count(pProp->GetValue().lVal);

  pProp = m_wndPropList.FindItemByData(pSrc->get_azim_count_ptr());
  if(pProp != NULL)
    pSrc->set_azim_count(pProp->GetValue().lVal);

  pProp = m_wndPropList.FindItemByData(pSrc->get_ensemble_size_ptr());
  if(pProp != NULL)
    pSrc->set_ensemble_size(pProp->GetValue().lVal);

  pProp = m_wndPropList.FindItemByData(pSrc->get_filename_ptr());
  if(pProp != NULL)
  {
    CString sFileName = (CString)pProp->GetValue();
    pSrc->set_filename((const char*)sFileName);
  }
}

void CPropertiesWnd::update_source_ctrls()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  EvaporatingParticle::CSource* pSrc = pObj->get_src();

  CMFCPropertyGridProperty* pTypeProp = m_wndPropList.FindItemByData(pSrc->get_src_type_ptr());
  bool bCone = false, bRadius = false, bRect = false, bSelReg = false, bCyl = false;
  if(pTypeProp != NULL)
  {
    CString cName = (CString)pTypeProp->GetValue();
    bCone = cName == pSrc->get_src_type_name(EvaporatingParticle::CSource::stCone);
    bRadius = cName == pSrc->get_src_type_name(EvaporatingParticle::CSource::stSpot);
    bRect = cName == pSrc->get_src_type_name(EvaporatingParticle::CSource::stRect);
    bSelReg = cName == pSrc->get_src_type_name(EvaporatingParticle::CSource::stSelReg);
    bCyl = cName == pSrc->get_src_type_name(EvaporatingParticle::CSource::stCylinder);
  }

  CMFCPropertyGridProperty* pInjProp = m_wndPropList.FindItemByData(pSrc->get_src_inject_type_ptr());
  bool bHomogen = (pInjProp != NULL) && ((CString)pInjProp->GetValue() == pSrc->get_inj_type_name(EvaporatingParticle::CSource::itHomogen));

  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pSrc->get_random_seed_ptr());
  if(pProp != NULL)
    pProp->Enable(!bHomogen);

// Absolute velocity controls:
  CMFCPropertyGridProperty* pUseGasVelProp = m_wndPropList.FindItemByData(pSrc->get_use_initial_gas_vel_ptr());
  bool bUseGasVel = (pUseGasVelProp != NULL) && (pUseGasVelProp->GetValue().boolVal);
  pProp = m_wndPropList.FindItemByData(pSrc->get_abs_vel_ptr());
  if(pProp != NULL)
    pProp->Enable(!bUseGasVel);

// Count controls:
  long nAzimCount = 1, nPitchCount = 1, nCount = 1;
  bool bEnablePitchAzim = bHomogen && !bRect && !bSelReg;
  pProp = m_wndPropList.FindItemByData(pSrc->get_pitch_count_ptr());
  if(pProp != NULL)
  {
    pProp->Enable(bEnablePitchAzim);
    if(bEnablePitchAzim)
      nPitchCount = pProp->GetValue().lVal;
  }

  pProp = m_wndPropList.FindItemByData(pSrc->get_azim_count_ptr());
  if(pProp != NULL)
  {
    pProp->Enable(bEnablePitchAzim);
    if(bEnablePitchAzim)
      nAzimCount = pProp->GetValue().lVal;
  }

  pProp = m_wndPropList.FindItemByData(pSrc->get_particles_count_ptr());
  if(pProp != NULL)
  {
    pProp->Enable(!bEnablePitchAzim);
    if(bEnablePitchAzim)
      pProp->SetValue(COleVariant(long(pSrc->get_particles_count())));
  }

// Cone
  pProp = m_wndPropList.FindItemByData(pSrc->get_cone_angle_ptr());
  if(pProp != NULL)
    pProp->Enable(bCone || bCyl);

// Circular spot and sphere
  pProp = m_wndPropList.FindItemByData(pSrc->get_radius_ptr());
  if(pProp != NULL)
    pProp->Enable(bRadius || bCyl);

// Rectangular spot
  pProp = m_wndPropList.FindItemByData(pSrc->get_rect_width_ptr());
  if(pProp != NULL)
    pProp->Enable(bRect);

  pProp = m_wndPropList.FindItemByData(pSrc->get_rect_height_ptr());
  if(pProp != NULL)
    pProp->Enable(bRect || bCyl);

// Ensemble size
  pTypeProp = m_wndPropList.FindItemByData(pObj->get_particle_type_ptr());
  bool bIonType = pObj->get_particle_type() == EvaporatingParticle::CTrack::ptIon; //(pTypeProp != NULL) && ((CString)pTypeProp->GetValue() == "Ions");
  pProp = m_wndPropList.FindItemByData(pSrc->get_ensemble_size_ptr());
  if(pProp != NULL)
    pProp->Enable(bIonType && !bSelReg);

// Selected region source:
  pProp = m_wndPropList.FindItemByData(pSrc->get_selected_rgn_names_ptr());
  if(pProp != NULL)
    pProp->Enable(bSelReg && pObj->is_ready());

  pProp = m_wndPropList.FindItemByData(pSrc->get_src_pos_ptr());
  if(pProp != NULL)
  {
    pProp->GetSubItem(0)->Enable(!bSelReg);
    pProp->GetSubItem(1)->Enable(!bSelReg);
    pProp->GetSubItem(2)->Enable(!bSelReg);
  }
}
