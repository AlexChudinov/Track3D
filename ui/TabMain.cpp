
#include "stdafx.h"

#include "PropertiesWnd.h"
#include "ParticleTracking.h"
#include "matrix3d.hpp"
#include "Button.h"


//---------------------------------------------------------------------------------------
// Particle type and data file name
//---------------------------------------------------------------------------------------
static const char* cSymPlanes[4] = { "None", "XY Only", "XZ Only", "Both XY and YZ" };

void CPropertiesWnd::add_type_ctrls()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

// Type of particles:
  CMFCPropertyGridProperty* pTypeGroup = new CMFCPropertyGridProperty(_T("Type of Particles"));
  COleVariant var(pObj->get_particle_type() == EvaporatingParticle::CTrack::ptDroplet ? _T("Droplets") : _T("Ions"));
  CMFCPropertyGridProperty* pType = new CMFCPropertyGridProperty(_T("Type of Particles"), var, _T("Specify type of particles, either droplets or ions."), pObj->get_particle_type_ptr());
  pType->AddOption(_T("Droplets"));
  pType->AddOption(_T("Ions"));
  pType->AllowEdit(FALSE);
  pTypeGroup->AddSubItem(pType);
  m_wndPropList.AddProperty(pTypeGroup);

// Input data filename:
  CMFCPropertyGridProperty* pDataGroup = new CMFCPropertyGridProperty(_T("Input"));
  static TCHAR BASED_CODE szFilter[] = _T("Geometry Data Files(*.geom)|*.geom|All Files(*.*)|*.*||");
// 1. Gas-dynamic data:
  CMFCPropertyGridFileProperty* pFileProp = new CMFCPropertyGridFileProperty(_T("Gas-Dynamics Data"), TRUE, pObj->get_filename(), _T("csv"),
    OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilter, _T("Specify location of the gas-dynamic data file."), pObj->get_filename_ptr());
  pDataGroup->AddSubItem(pFileProp);

  // AC 27.06.2016 Gas parameters
  CMFCPropertyGridProperty* pGasProp = new CMFCPropertyGridProperty(_T("Gas Molecular Mass, g/mol"), COleVariant(pObj->get_molar_mass() / Const_AMU_CGS), _T("Specify molecular mass of the buffer gas."), pObj->get_molar_mass_ptr());
  pDataGroup->AddSubItem(pGasProp);

  m_wndPropList.AddProperty(pDataGroup);

// File output:
  EvaporatingParticle::COutputEngine& outEng = pObj->get_output_engine();

  CMFCPropertyGridProperty* pOutGroup = new CMFCPropertyGridProperty(_T("Output to a File"));
  CMFCPropertyGridProperty* pOutProp = new CMFCPropertyGridProperty(_T("Enable File Output"), (_variant_t)outEng.get_enable_file_output(), _T("Turns ON/OFF the file output."), outEng.get_enable_file_output_ptr());
  pOutGroup->AddSubItem(pOutProp);

  CMFCPropertyGridProperty* pOutTypeProp = new CMFCPropertyGridProperty(_T("Collect tracks by initial radius"), (_variant_t)outEng.get_enable_ens_by_radius(), _T("If this is 'true' the output tracks are grouped into bins by their initial radii."), outEng.get_enable_ens_by_radius_ptr());
  pOutGroup->AddSubItem(pOutTypeProp);

  CMFCPropertyGridProperty* pCountProp = new CMFCPropertyGridProperty(_T("Count of radial ranges"), COleVariant((long)outEng.get_ens_by_radius_count()), _T("Count of initial radius bins, by which the output ensembles will be built."), outEng.get_ens_by_radius_count_ptr());
  pOutGroup->AddSubItem(pCountProp);

  CMFCPropertyGridProperty* pQ00Prop = new CMFCPropertyGridProperty(_T("Output only Q00 tracks"), (_variant_t)outEng.get_restrict_output(), _T("If true, only those tracks will be output to a file which have passed to the Q00 region."), outEng.get_restrict_output_ptr());
  pOutGroup->AddSubItem(pQ00Prop);
  m_wndPropList.AddProperty(pOutGroup);

// Saving track data to the *.tsk file:
  CMFCPropertyGridProperty* pSaveTrackGroup = new CMFCPropertyGridProperty(_T("Save track data onto the disk"));
  CMFCPropertyGridProperty* pEnableSaveTrack = new CMFCPropertyGridProperty(_T("Enable"), (_variant_t)pObj->get_enable_save_tracks(), _T("If this flag is ON the File/Save... command saves the full tracks data on the disk."), pObj->get_enable_save_tracks_ptr());
  pSaveTrackGroup->AddSubItem(pEnableSaveTrack);
  m_wndPropList.AddProperty(pSaveTrackGroup);

// Symmetry:
  CMFCPropertyGridProperty* pSymGroup = new CMFCPropertyGridProperty(_T("Symmetry"));
  int nSymPlanes = pObj->get_sym_plane();
  if((nSymPlanes & EvaporatingParticle::CTracker::spXY) && (nSymPlanes & EvaporatingParticle::CTracker::spXZ))
    var = COleVariant(_T(cSymPlanes[3]));
  else if(nSymPlanes & EvaporatingParticle::CTracker::spXZ)
    var = COleVariant(_T(cSymPlanes[2]));
  else if(nSymPlanes & EvaporatingParticle::CTracker::spXY)
    var = COleVariant(_T(cSymPlanes[1]));
  else
    var = COleVariant(_T(cSymPlanes[0]));

  CMFCPropertyGridProperty* pSym = new CMFCPropertyGridProperty(_T("Type of Symmetry"), var, _T("Specify symmetry planes."), pObj->get_sym_plane_ptr());
  for(UINT i = 0; i < 4; i++)
    pSym->AddOption(_T(cSymPlanes[i]));

  pSymGroup->AddSubItem(pSym);
  m_wndPropList.AddProperty(pSymGroup);

// Mesh transformation:
  CMFCPropertyGridProperty* pTransGroup = new CMFCPropertyGridProperty(_T("Mesh Transformation"));
  EvaporatingParticle::CTransform& trans = pObj->get_transform();
// Enable
  CMFCPropertyGridProperty* pProp = new CMFCPropertyGridProperty(_T("Enable"), (_variant_t)trans.get_enable(), _T("Turns ON/OFF the mesh transformation."), trans.get_enable_ptr());
  pTransGroup->AddSubItem(pProp);
// Rotation axis:
  var = COleVariant(_T(trans.axis_name(trans.get_rot_axis())));
  pProp = new CMFCPropertyGridProperty(_T("Rotation Axis"), var, _T("Specify rotation axis."), trans.get_rot_axis_ptr());
  for(int j = 0; j < EvaporatingParticle::CTransform::axesCount; j++)
    pProp->AddOption(_T(trans.axis_name(j)));

  pTransGroup->AddSubItem(pProp);
// Rotation angle:
  pProp = new CMFCPropertyGridProperty(_T("Rotation angle, degree"), COleVariant(trans.get_rot_angle() * Const_RadianToDegree), _T("Mesh rotation angle in degrees."), trans.get_rot_angle_ptr());
  pTransGroup->AddSubItem(pProp);

  m_wndPropList.AddProperty(pTransGroup);
}

void CPropertiesWnd::set_type_data()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

// Type of particles:
  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(pObj->get_particle_type_ptr());
  if(pProp != NULL)
  {
    CString cType = (CString)pProp->GetValue();
    if(cType == "Droplets")
      pObj->set_particle_type(EvaporatingParticle::CTrack::ptDroplet);
    else if(cType == "Ions")
      pObj->set_particle_type(EvaporatingParticle::CTrack::ptIon);
  }

// Input data filename, gas-dynamic data:
  pProp = m_wndPropList.FindItemByData(pObj->get_filename_ptr());
  if(pProp != NULL)
  {
    CString cFile = (CString)pProp->GetValue();
    std::string str = std::string(CT2CA(cFile));
    pObj->set_filename(str.c_str());
  }

// Gas molecular mass AC 27/06/2016
  pProp = m_wndPropList.FindItemByData(pObj->get_molar_mass_ptr());
	if(pProp != NULL)
    pObj->set_molar_mass(Const_AMU_CGS * pProp->GetValue().dblVal);

// File output:
  EvaporatingParticle::COutputEngine& outEng = pObj->get_output_engine();

  pProp = m_wndPropList.FindItemByData(outEng.get_enable_file_output_ptr());
  if(pProp != NULL)
    outEng.set_enable_file_output(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(outEng.get_enable_ens_by_radius_ptr());
  if(pProp != NULL)
    outEng.set_enable_ens_by_radius(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(outEng.get_ens_by_radius_count_ptr());
  if(pProp != NULL)
    outEng.set_ens_by_radius_count(pProp->GetValue().lVal);

  pProp = m_wndPropList.FindItemByData(outEng.get_restrict_output_ptr());
  if(pProp != NULL)
    outEng.set_restrict_output(pProp->GetValue().boolVal);

// Saving track data to the *.tsk file:
  pProp = m_wndPropList.FindItemByData(pObj->get_enable_save_tracks_ptr());
  if(pProp != NULL)
    pObj->set_enable_save_tracks(pProp->GetValue().boolVal);

// Symmetry:
  pProp = m_wndPropList.FindItemByData(pObj->get_sym_plane_ptr());
  if(pProp != NULL)
  {
    CString cType = (CString)pProp->GetValue();
    int nSymPlanes = 0;
    if(cType == cSymPlanes[1])
      nSymPlanes = EvaporatingParticle::CTracker::spXY;
    else if(cType == cSymPlanes[2])
      nSymPlanes = EvaporatingParticle::CTracker::spXZ;
    else if(cType == cSymPlanes[3])
      nSymPlanes = EvaporatingParticle::CTracker::spXY | EvaporatingParticle::CTracker::spXZ;

    pObj->set_sym_plane(nSymPlanes);
  }

// Mesh transformation:
  EvaporatingParticle::CTransform& trans = pObj->get_transform();
  pProp = m_wndPropList.FindItemByData(trans.get_enable_ptr());
  if(pProp != NULL)
    trans.set_enable(pProp->GetValue().boolVal);

  pProp = m_wndPropList.FindItemByData(trans.get_rot_axis_ptr());
  if(pProp != NULL)
  {
    CString cAxis = (CString)pProp->GetValue();
    for(int i = 0; i < EvaporatingParticle::CTransform::axesCount; i++)
    {
      if(cAxis == trans.axis_name(i))
      {
        trans.set_rot_axis(i);
        break;
      }
    }
  }

  pProp = m_wndPropList.FindItemByData(trans.get_rot_angle_ptr());
  if(pProp != NULL)
    trans.set_rot_angle(pProp->GetValue().dblVal * Const_DegreeToRadian);
}

void CPropertiesWnd::update_type_ctrls()
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  EvaporatingParticle::COutputEngine& outEng = pObj->get_output_engine();

  bool bEnable = false;
  CMFCPropertyGridProperty* pProp = m_wndPropList.FindItemByData(outEng.get_enable_file_output_ptr());
  if(pProp != NULL)
    bEnable = pProp->GetValue().boolVal;

  pProp = m_wndPropList.FindItemByData(outEng.get_enable_ens_by_radius_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(outEng.get_ens_by_radius_count_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  pProp = m_wndPropList.FindItemByData(outEng.get_restrict_output_ptr());
  if(pProp != NULL)
    pProp->Enable(bEnable);

  bool bEnableButton = pObj->is_ready();
  pProp = m_wndPropList.FindItemByData((DWORD_PTR)&outEng);
  if(pProp != NULL)
    pProp->Enable(bEnableButton);
}

