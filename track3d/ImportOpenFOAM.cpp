#include "stdafx.h"

#include <stdio.h>

#include "math.h"
#include "float.h"

#include "ImportOpenFOAM.h"

#include "ParticleTracking.h"
#include "constant.hpp"


namespace EvaporatingParticle
{

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
CImportOpenFOAM::CImportOpenFOAM()
{
  set_default();
}

CImportOpenFOAM::~CImportOpenFOAM()
{
}

void CImportOpenFOAM::set_default()
{
  m_nStep = opFirstStep;
  m_fShiftX = 5.85; // cm, the shift in S-Lens.
}

bool CImportOpenFOAM::do_import()
{
  terminate(false);
  switch(m_nStep)
  {
    case opFirstStep: return first_step();
    case opSecondStep: return second_step();
    case opRestore: return restore_ansys_data();
  }

  return true;
}

//-------------------------------------------------------------------------------------------------
// First (preliminary) step. Reading OpenFOAM file, searching for proper nodes in a short mesh
// used for OpenFOAM export and changing corresponding *.var file on disk. Input: *.geom file of a
// short mesh as the "Gas-Dynamic data file" and a *.csv OpenFOAM file as the "Import data file".
// At this step the proper nodes are searched by a simple comparison of the distance with a given
// tolerance eps ~ 0.001 mm.
//-------------------------------------------------------------------------------------------------
bool CImportOpenFOAM::first_step()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
// There must be nodes of the "short" mesh used by the ExportOpenFOAM class to create points, faces, owner, neighbour and boundary files.
  CNodesCollection& vNodes = pObj->get_nodes();
  double fMolarMass = pObj->get_molar_mass();   // g.

  size_t nNodeCount = vNodes.size();
  if(nNodeCount == 0)
    return false;

  set_job_name("OpenFOAM import, first step..."); // the handlers are supposed to be set in the calling CTracker class.
  set_progress(0);

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, m_sFile.c_str(), (const char*)("r"));
  if(nErr != 0 || pStream == 0)
    return false;

  char sHeader[256];
  int nRes = fscanf(pStream, "%s", sHeader);

  bool bOK = true;
  UINT nFailsCount = 0, nSuccessCount = 0, i;
  float fT, fN, fP, fVx, fVy, fVz, fX, fY, fZ;
  while(true)
  {
    nRes = fscanf_s(pStream, "%f, %e, %f, %f, %f, %f, %f, %f, %f", &fT, &fN, &fP, &fVx, &fVy, &fVz, &fX, &fY, &fZ);
    if(nRes == EOF)
      break;

    convert_to_cgs(fN, fVx, fVy, fVz, fX, fY, fZ);

    CNode3D* pNode = find_node(fX, fY, fZ);
    if(pNode == NULL)
    {
      nFailsCount++;
      continue;
    }

// Substitute node parameters:
    pNode->dens = fN * fMolarMass;
    pNode->press = fN * Const_Boltzmann * fT;         // P = nkT.
    pNode->vel = Vector3D(fVx, fVy, fVz);
    pNode->temp = fT;
    nSuccessCount++;

    i = nSuccessCount + nFailsCount;
    if(i % 100 == 0)
      set_progress(int(0.5 + 100. * i / nNodeCount));
    if(get_terminate_flag())
    {
      bOK = false;
      break;
    }
  }

  fclose(pStream);

  if(bOK)
    update_var_file();

  print_statistics(nSuccessCount, nFailsCount);
  return bOK;
}

void CImportOpenFOAM::convert_to_cgs(float& fN, float& fVx, float& fVy, float& fVz, float& fX, float& fY, float& fZ) const
{
  fN *= (float)SI_to_CGS_NumDens;

  fVx *= (float)SI_to_CGS_Vel;
  fVy *= (float)SI_to_CGS_Vel;
  fVz *= (float)SI_to_CGS_Vel;

  fX *= (float)SI_to_CGS_Len;
  fX += (float)m_fShiftX;

  fY *= (float)SI_to_CGS_Len;
  fZ *= (float)SI_to_CGS_Len;
}

void CImportOpenFOAM::convert_to_si(float& fPress, float& fDens, float& fVx, float& fVy, float& fVz) const
{
  fPress *= (float)CGS_to_SI_Press;
  fDens *= (float)CGS_to_SI_Dens;

  fVx *= (float)CGS_to_SI_Vel;
  fVy *= (float)CGS_to_SI_Vel;
  fVz *= (float)CGS_to_SI_Vel;
}

static const double scfTol = 0.0001; // tolerance, cm.

CNode3D* CImportOpenFOAM::find_node(float fX, float fY, float fZ) const
{
  CNodesCollection& vNodes = CParticleTrackingApp::Get()->GetTracker()->get_nodes();

  Vector3D vPos;
  CNode3D* pNode = NULL;
  size_t nNodeCount = vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    pNode = vNodes.at(i);
    vPos = pNode->pos;
    if(fabs(vPos.x - fX) > scfTol || fabs(vPos.y - fY) > scfTol || fabs(vPos.z - fZ) > scfTol)
      continue;

    return pNode;
  }

  return NULL;
}

static const char scHeader[] = "   Abs Press      Density         Temp     Dyn Visc   Therm Cond           Cp           Vx           Vy           Vz        DC Ex        DC Ey        DC Ez        RF Ex        RF Ey        RF Ez";

bool CImportOpenFOAM::update_var_file()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
// There must be nodes of the "short" mesh used by the ExportOpenFOAM class to create points, faces, owner, neighbour and boundary files:
  CNodesCollection& vNodes = pObj->get_nodes();

  set_job_name("OpenFOAM import: Updating gas-dynamic data on disk...");
  set_progress(0);

  std::string cBase = COutputEngine::get_base_name(pObj->get_filename());

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (cBase + "var").c_str(), (const char*)("w+"));
  if(nErr != 0 || pStream == 0)
    return false;

  fprintf(pStream, "%s\n", scHeader);

  Vector3D vE(0, 0, 0), vErf(0, 0, 0);
  float fPress, fDens, fTemp, fDynVisc = 0.0f, fThermCond = 0.0f, fCp = 0.0f, fVx, fVy, fVz;

  CNode3D* pNode = NULL;
  size_t nNodeCount = vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    pNode = vNodes.at(i);
    fPress = pNode->press;
    fDens = pNode->dens;
    fTemp = pNode->temp;
    fVx = pNode->vel.x;
    fVy = pNode->vel.y;
    fVz = pNode->vel.z;

    convert_to_si(fPress, fDens, fVx, fVy, fVz);  // the parameters are stored on the disk in SI units.

    fprintf(pStream, "%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
      fPress, fDens, fTemp, fDynVisc, fThermCond, fCp, fVx, fVy, fVz, vE.x, vE.y, vE.z, vErf.x, vErf.y, vErf.z);

    if(i % 100 == 0)
      set_progress(int(0.5 + 100. * i / nNodeCount));
  }

  fclose(pStream);
  return true;
}

void CImportOpenFOAM::print_statistics(UINT nSuccess, UINT nFail)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  std::string cPath = COutputEngine::get_full_path(pObj->get_filename());

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (cPath + "import_statistics.txt").c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  fprintf(pStream, "Count of successes: %d, count of failures: %d, total: %d", nSuccess, nFail, nSuccess + nFail);

  fclose(pStream);
}

//-------------------------------------------------------------------------------------------------
// Second (run-time) step. At this step pressure, temperature and velocity at the nodes of the full
// mesh are substituted by corresponding values interpolated from the short mesh during run-time.
//-------------------------------------------------------------------------------------------------
bool CImportOpenFOAM::second_step()
{
// There must be nodes of the "full" mesh suitable for tracking particles.
  CNodesCollection& vNodes = CParticleTrackingApp::Get()->GetTracker()->get_nodes();
  size_t nNodeCount = vNodes.size();
  if(nNodeCount == 0)
    return false;

  CTracker AuxObj(true);
// At this step m_sFile must point to the short geometry, the "*.var" file of which has been earlier modified by
// CImportOpenFOAM::first_step(). Since that time it contains pressure, temperature and velocity from OpenFOAM.
  AuxObj.set_filename(m_sFile.c_str());
  AuxObj.set_handlers(m_hJobNameHandle, m_hProgressBarHandle);

  if(!AuxObj.read_geometry())
    return false;

  if(!AuxObj.read_gasdyn_data())
    return false;

// Now AuxObj contains OpenFOAM data in the vertices of the "short" geometry. The task of the second step is to
// interpolate these data into the vertices of the "full" mesh. Do not forget about the dummy zero node.

  set_job_name("OpenFOAM import: Interpolating data into mesh vertices...");
  set_progress(0);

  Vector3D vPos;
  CElem3D* pElem = NULL;
  CNode3D* pNode;
  CNode3D node;   // this is just a container for interpolated data.

  CBox box = AuxObj.get_box();
  for(size_t i = 1; i < nNodeCount; i++)
  {
    pNode = vNodes.at(i);
    vPos = pNode->pos;
    if(!box.inside(vPos) || !AuxObj.interpolate(vPos, 0, 0, node, pElem))
      continue;

// Modification of the main gas-dynamic parameters by the imported data:
    pNode->press = node.press;
    pNode->dens = node.dens;
    pNode->temp = node.temp;
    pNode->vel = node.vel;

    if(i % 100 == 0)
      set_progress(int(0.5 + 100. * i / nNodeCount));
    if(get_terminate_flag())
      return false;
  }

  return true;
}

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
bool CImportOpenFOAM::restore_ansys_data()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  pObj->terminate(false);
  pObj->set_handlers(m_hJobNameHandle, m_hProgressBarHandle);
  bool bRes = pObj->read_gasdyn_data();
  pObj->set_handlers(NULL, NULL);
  return bRes;
}

//-------------------------------------------------------------------------------------------------
// Streamability.
//-------------------------------------------------------------------------------------------------
void CImportOpenFOAM::save(CArchive& ar)
{
  const UINT nVersion = 3;
  ar << nVersion;

// Data file:
  CString cFileName(m_sFile.c_str());
  ar << cFileName;

  ar << m_fShiftX;
  ar << m_nStep;
}

void CImportOpenFOAM::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

// Data file:
  CString cFileName;
  ar >> cFileName;
  set_filename((const char*)cFileName);

  ar >> m_fShiftX;

  if(nVersion >= 3)
    ar >> m_nStep;

  if(nVersion == 1) // since version 2 the environment gas molar mass was moved to the CTracker class.
  {
    double fMolarMass;
    ar >> fMolarMass;
  }
}

};  // namespace EvaporatingParticle
