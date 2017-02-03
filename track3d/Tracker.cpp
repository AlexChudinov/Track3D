
#include "stdafx.h"

#include "tchar.h"
#include <stdio.h>
#include <direct.h>
#include "float.h"

#include "Tracker.hpp"
#include "BarnesHut.h"
#include "Symmetry.hpp"

#include "EvaporationModel.h"
#include "CalcThread.h"
#include <algorithm>

#include "mathematics.h"

#include "ParticleTracking.h"
#include "PropertiesWnd.h"
#include "MainFrm.h"

#include "ExecutionDialog.h"

#include <libIntegrators.h>


namespace EvaporatingParticle
{
//-------------------------------------------------------------------------------------------------
// CTracker - the main class for data reading from ANSYS data file and tracking particles taking
//            into account evaporation and heat exchange with the environment.
//-------------------------------------------------------------------------------------------------
CTracker::CTracker()
  : m_pEvaporModel(NULL), m_pBarnesHut(NULL)
{
  set_default();
}

CTracker::~CTracker()
{
  if(m_pEvaporModel != NULL)
    delete m_pEvaporModel;

  if(m_pBarnesHut != NULL)
    delete m_pBarnesHut;

  clear();
}

void CTracker::clear()
{
  size_t nElemCount = m_vElems.size();
  for(size_t j = 0; j < nElemCount; j++)
  {
    set_status("Deleting elements", 100 * j / nElemCount);
    delete m_vElems.at(j);
  }

  size_t nRegCount = m_vRegions.size();
  for(size_t k = 0; k < nRegCount; k++)
  {
    CRegion* pReg = m_vRegions.at(k);
    size_t nFaceCount = pReg->vFaces.size();
    for(size_t l = 0; l < nFaceCount; l++)
      delete pReg->vFaces.at(l);

    set_status("Deleting regions", 100 * k / nRegCount);
    delete pReg;
  }

  size_t nNodeCount = m_vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    set_status("Deleting nodes", 100 * i / nNodeCount);
    delete m_vNodes.at(i);
  }

  m_vNodes.clear();
  m_vElems.clear();
  m_vRegions.clear();
  m_vExtRegions.clear();

  m_Transform.set_default();

  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->clear();
  pDrawObj->invalidate_all();

  CCrossSectColl* pColl = CParticleTrackingApp::Get()->GetPlanes();
  size_t nPlanesCount = pColl->size();
  for(size_t l = 0; l < nPlanesCount; l++)
    delete pColl->at(l);

  pColl->clear();

  CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
  pFields->clear_fields();

  set_status("Ready", -1);
}

void CTracker::set_default()
{
  m_nType = CTrack::ptDroplet;
  m_bUseMultiThread = true;
  m_nIntegrType = intModMidpnt;
  m_bOldIntegrator = false;
  m_bUseRadialCoulomb = true;
  m_bSaveTracks = false;

  m_bReady = false;
  m_bTerminate = false;

  m_bConv2CGS = true;

  double fTimeStep = 2.e-8;       // sec
  set_time_step(fTimeStep);

  m_fInitD = 1.0e-4;      // cm, 1 mcm particles by default.
  m_fPartDens = 1.0;      // g/cm3, water by default.

  m_fMolarMass = Const_Molar_Mass_Air * Const_AMU_CGS;
  m_fInitMass = get_particle_mass(m_fInitD);  // g.
  m_nSymPlanes = spXY | spXZ;

// Electrostatics:
  m_bEnableField = true;
  m_fCharge = 1.e+4 * Const_Charge_CGS;
  m_fAmplDC = 1.; // this is a multiplier, as ANSYS data provide DC field for 2500 V at the emitter.

// Restrictions:
  m_fMinD = 1.0e-6;       // cm.
  m_fMinMass = get_particle_mass(m_fMinD);  // g.
  m_fMaxIntegrTime = 1.0; // sec.

// Evaporation:
  m_nEvaporModelType = 0; // none.
  create_evapor_model();

// Ion type of particles:
  m_fIonMobility = 1. / SI_to_CGS_Voltage;  // 1 cm2/(s*V) in CGSE.
  const double cfA2 = Const_Angstrem_CGS * Const_Angstrem_CGS;
  set_ion_cross_section(160 * cfA2);        // collision cross-section, 160 squared angstrem, 1.6e-14 cm^2.
  m_bVelDependent = false;                  // if true, the collision cross-section is inversely proportional to the relative velocity.
  set_ion_mass(Const_AMU_CGS * 16000);      // protein molecule.
  m_fActEnergy = 1.0;       // eV.

// Coulomb effects:
  m_bEnableCoulomb = false;
  m_fInitBunchRadius = 0.029;  // cm, corresponds to 0.29 mm - full radius of the round capillary.
  set_full_current(100 * Const_nA_to_CGSE); // 100 nA by default;

  m_bAxialSymm = true;

// Iterational calculations of the Coulomb effects:
  m_nIterCount = 10;
  m_fRadialCoulombX = 10000;
// Barnes-Hut parameters:
  m_fBHDist = 1.5;
  m_fBHCritR = 0.001; // cm.
  m_nMaxRecDepth = 12;
  m_bEnableQuadTerms = false;

// Radio-frequency field:
  m_bEnableRF = false;
  m_fAmplRF = 80.;          // V, this is only a scaling factor.
  set_rf_frequency(6.3e+5); // 630 kHz by default.
// RF in Q00:
  m_fAmplRF_Q00 = 200.;
  set_rf_Q00_freq(3.0e+6);  // 3 MHz by default.
  m_fX_Q00 = 10000.; // Q00 region does not exists by default.
// RF in flatapole (Q0):
  m_fAmplRF_Q0 = 350.;
  set_rf_flatapole_freq(3.0e+6);  // 3 MHz by default.
  m_fX_Q0 = 10000.; // Q0 region does not exists by default.
}

//-------------------------------------------------------------------------------------------------
// Particle tracking
//-------------------------------------------------------------------------------------------------
static STEPPER_TYPE_ID GetIntegrType(int nType)
{
  switch(nType)
  {
    case CTracker::intExplEuler: return BOOST_EXPLICIT_EULER;
    case CTracker::intModMidpnt: return BOOST_MODIFIED_MIDPOINT;
    case CTracker::intRK2: return RUNGE_KUTTA2;
    case CTracker::intRK4: return BOOST_RUNGE_KUTTA4;
  }

  return UNSUPPORTED_TYPE;
}

static void SetZeroDeriv(double* pTimeDeriv, ULONG nSize)
{
  for(UINT i = 0; i < nSize; i++)
    pTimeDeriv[i] = 0;
}

void CTracker::GetTimeDeriv(const void* pData, const double* pItemState, double* pTimeDeriv, const double* pTime)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  ULONG nStateSize = pObj->get_particle_type() == CTrack::ptDroplet ? DROPLET_STATE_SIZE : ION_STATE_SIZE;

  CIntegrInterface* pI = (CIntegrInterface*)pData;
  if(pI->nElemId < 0)
  {
    SetZeroDeriv(pTimeDeriv, nStateSize);
    return;
  }

  CElem3D* pElem = pObj->get_elems().at(pI->nElemId);

  Vector3D vAccel;
  Vector3D vPos(pItemState[0], pItemState[1], pItemState[2]);
  Vector3D vVel(pItemState[3], pItemState[4], pItemState[5]);
  double fIonTemp = pItemState[6];

// After this call vPos and vVel are in the gas-dynamic domain. If the reflection is really
// done, the bReflDone flag will be "true". Note that vAccel is a dummy parameter here.
  bool bReflDone = pObj->sym_corr(vPos, vVel, vAccel);

  CNode3D node;
// Compute all ANSYS fields at pItem->pos and place them into the node structure.
  if(!pObj->interpolate(vPos, node, pElem))
  {
    pI->nElemId = -1;  // if interpolate(...) fails terminate integration of this track.
    pI->bOk = false;
    SetZeroDeriv(pTimeDeriv, nStateSize);
    return;
  }

  pI->nElemId = pElem->nInd;
  double fTime = *pTime;

  if(pObj->get_particle_type() == CTrack::ptIon)
  {
    double fExpCoeff;
    vAccel = pObj->get_ion_accel(node, vVel, fTime, pObj->m_fTimeStep, pI->fPhase, pI->fCurr, fExpCoeff);
// vAccel is the reflected particle acceleration. Turn back the position, velocity and acceleration:
    if(bReflDone)
      pObj->sym_corr(vPos, vVel, vAccel, true);

    double fTionInf;
    pTimeDeriv[6] = pObj->get_dTi(node, vVel, fIonTemp, fExpCoeff, fTionInf);
    pI->fTionInf = fTionInf;

    pTimeDeriv[7] = 0;
  }
  else
  {
    double fRe;
    double fTemp = pItemState[6];
    double fMass = pItemState[7];
    if(fMass < pObj->m_fMinMass)
    {
      pI->bOk = false;
      SetZeroDeriv(pTimeDeriv, nStateSize);
      return;
    }

    double fD = pObj->get_particle_diameter(fMass);
    vAccel = pObj->get_accel(node, vVel, fMass, fD, fTime, fRe);  // the Reynolds number is computed here.
// vAccel is the reflected particle acceleration. Turn back the position, velocity and acceleration:
    if(bReflDone)
      pObj->sym_corr(vPos, vVel, vAccel, true);

    CEvaporationModel* pEvaporModel = pObj->get_evapor_model();
    pTimeDeriv[6] = pEvaporModel->get_cooling_rate(node, fTemp, fD, fRe);
    pTimeDeriv[7] = pEvaporModel->get_evaporation_rate(node, fTemp, fD, fRe);
  }

  pTimeDeriv[0] = vVel.x;
  pTimeDeriv[1] = vVel.y;
  pTimeDeriv[2] = vVel.z;

  pTimeDeriv[3] = vAccel.x;
  pTimeDeriv[4] = vAccel.y;
  pTimeDeriv[5] = vAccel.z;
}

void CTracker::get_output_freq(UINT& nOutFreq) const
{
  nOutFreq = UINT(m_OutputEngine.get_output_time_step() / m_fTimeStep);
  if(nOutFreq == 0)
    nOutFreq = 1;
}

void CTracker::single_thread_calculate()
{
  if(!prepare())
  {
    m_bResult = false;
    return;
  }

  bool bIter = (m_nType == CTrack::ptIon) && m_bEnableCoulomb && !m_bAxialSymm;
  if(bIter)
    do_iterations();
  else
    do_track();

  relax();
}

void CTracker::do_track()
{
  set_job_name("Calculating particles...");
  set_progress(0);

  UINT nOutFreq;
  get_output_freq(nOutFreq);
  init_currents();

  const ULONG nStateSize = m_nType == CTrack::ptDroplet ? DROPLET_STATE_SIZE : ION_STATE_SIZE;
  STEPPER_TYPE_ID nIntegrType = GetIntegrType(get_integr_type());

  size_t nPartCount = m_Tracks.size();
  for(size_t i = 0; i < nPartCount; i++)
  {
    size_t nStep = 1; // to prevent from writing the first item to the track twice.
    CTrack& track = m_Tracks.at(i);
    if(m_bOldIntegrator)
    {
      CBaseTrackItem* pItem = track.at(0)->copy();
      while(true)
      {
        bool bOK = m_nType == CTrack::ptDroplet ? 
          do_time_step(pItem) : do_ion_time_step(pItem, track.get_phase(), track.get_current());

        if(!bOK)
        {
          track.push_back(pItem->copy());  // insert the last item anyway.
          break;  // the track went out of the mesh.
        }

        if(nStep % nOutFreq == 0)
          track.push_back(pItem->copy());

        if(track_is_over(pItem->time, pItem->get_mass(), pItem->get_temp()))
          break;  // the integration time exceeded the limit or the Rayleigh criterion was met.

        nStep++;
      }

      delete pItem;
    }
    else
    {
      double fMass, fTemp;
      CBaseTrackItem* pItem = track.at(0);
      double fTime = pItem->time;
      double pState[max(ION_STATE_SIZE, DROPLET_STATE_SIZE)];
      pItem->state(pState);
      CIntegrInterface data(pItem->nElemId, track.get_index(), track.get_phase(), track.get_current());
      void* pI = create_integrator_interface(nStateSize, nIntegrType, (const void*)&data, CTracker::GetTimeDeriv);
      while(true)
      {
        do_integrator_step(pI, pState, &fTime, &m_fTimeStep);

        if(!data.bOk)
        {
          CBaseTrackItem* pCurrItem = CBaseTrackItem::create(m_nType, -1, fTime, pState);
          track.push_back(pCurrItem);  // insert the last item anyway.
          break;  // the track went out of the mesh.
        }

        if(nStep % nOutFreq == 0)
        {
          CBaseTrackItem* pCurrItem = CBaseTrackItem::create(m_nType, data.nElemId, fTime, pState);
          if(m_nType == CTrack::ptIon)
          {
            CIonTrackItem* pIonItem = (CIonTrackItem*)pCurrItem;
            pIonItem->tempinf = data.fTionInf;
          }

          track.push_back(pCurrItem);
        }

        fTemp = pState[6];
        fMass = m_nType == CTrack::ptDroplet ? pState[7] : 0;
        if(track_is_over(fTime, fMass, fTemp))
          break;  // the integration time exceeded the limit or the Rayleigh criterion was met.

        fTime += m_fTimeStep;
        nStep++;
      }

      delete_integrator_interface(pI);
    }

    set_progress(int(0.5 + 100. * (i + 1) / nPartCount));
    if(m_bTerminate)
      return;
  }
}

void CTracker::do_iterations()
{
  for(UINT i = 0; i < m_nIterCount; i++)
  {
    do_track();
    if(m_bTerminate)
      return;

    if(i == m_nIterCount - 1)
      break;
    if(i == 0)
      m_OutputEngine.prepare_current_output();

    if((m_nIterCount > 4) && (i >= m_nIterCount - 4))
      m_OutputEngine.add_current();    // average output currents over 4 latest iterations.

    CalcThreadVector vThreads;
    create_BH_object(vThreads, i + 1);
    clear_tracks(false);    // keep the initial positions only.
  }
}

bool CTracker::create_BH_object(CalcThreadVector& vThreads, UINT nIter)
{
  set_job_name("Building Space Charge distribution volume...");
  set_progress(0);

  if(m_pBarnesHut != NULL)
    delete m_pBarnesHut;

  size_t nTrackCount = m_Tracks.size();
  if(nTrackCount < 1)
    return false;

  m_pBarnesHut = new CBarnesHut();
  m_pBarnesHut->set_dist_coeff(m_fBHDist);
  m_pBarnesHut->set_crit_radius(m_fBHCritR);
  m_pBarnesHut->set_max_rec_depth(m_nMaxRecDepth);
  m_pBarnesHut->set_enable_quad_terms(m_bEnableQuadTerms);
  m_pBarnesHut->set_sym_type(get_symmetry_type());

  Vector3D vCenter;
  double fEdge, fMinX, fMaxX;
  get_BH_cube(vCenter, fEdge, fMinX, fMaxX);
  m_pBarnesHut->create_main_cell(vCenter, fEdge); 

  double fCurrPerTrack = get_full_current_at(nIter) / nTrackCount;   // CGSE.

  if(m_bTerminate || m_hJobNameHandle == NULL || m_hProgressBarHandle == NULL)
    return false;

  m_SpaceChargeDist.clear();
  m_SpaceChargeDist.set_run_time_data(fMinX, fMaxX, fCurrPerTrack);
  m_SpaceChargeDist.set_handlers(m_hJobNameHandle, m_hProgressBarHandle);

  if(!m_SpaceChargeDist.set_BH_object(m_pBarnesHut))
  {
    m_bTerminate = true;
    return false;
  }
  
  m_pBarnesHut->prepare(vThreads);

  m_SpaceChargeDist.set_handlers(NULL, NULL);
  return true;
}

int CTracker::get_symmetry_type() const
{
  int nSymXY = m_nSymPlanes & spXY;
  int nSymXZ = m_nSymPlanes & spXZ;

  if((nSymXY == 0) && (nSymXZ == 0))
    return CBarnesHut::symNone;
  if((nSymXY != 0) && (nSymXZ == 0))
    return CBarnesHut::symXYonly;
  if((nSymXY == 0) && (nSymXZ != 0))
    return CBarnesHut::symXZonly;
  if((nSymXY != 0) && (nSymXZ != 0))
    return CBarnesHut::symBoth;

  return CBarnesHut::symNone;
}

void CTracker::get_BH_cube(Vector3D& c, double& edge, double& fMinX, double& fMaxX) const
{
  c = m_Box.get_center();
  Vector3D vMin = m_Box.vMin;
// Temporary limitation, which is, in general, incorrect:
  Vector3D vSrcPos = m_Src.get_src_pos();
  if(vSrcPos.x > vMin.x)
    vMin.x = vSrcPos.x;

  Vector3D vMax = m_Box.vMax;
  double vEdges[3] = { vMax.x - vMin.x, vMax.y - vMin.y, vMax.z - vMin.z };
  edge = vEdges[0];
  for(UINT i = 1; i < 3; i++)
  {
    if(edge < vEdges[i])
      edge = vEdges[i];
  }

  fMinX = vMin.x;
  fMaxX = vMax.x;
}

bool CTracker::capture_save_image(UINT nIter)
{
  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->invalidate_tracks();

  CMainFrame* pMainFrame = (CMainFrame*)(CParticleTrackingApp::Get()->GetMainWnd());
  CView* pView = pMainFrame->GetActiveView();
  if(pView == NULL)
    return false;

  show_dlg(SW_HIDE);
  Sleep(500);
  pView->SendMessage(WM_PAINT);
  pDrawObj->capture_image();
  show_dlg(SW_SHOW);

  std::string cPath = COutputEngine::get_full_path(get_filename());
  std::string cName("iter_#");
  char buff[4];
  std::string cIter(itoa(nIter, buff, 10));
  std::string cExt(".png");
  std::string cFileName = cPath + cName + cIter + cExt;

  pDrawObj->save_image(cFileName.c_str());
  return true;
}

//-------------------------------------------------------------------------------------------------
// Multi-threading support
//-------------------------------------------------------------------------------------------------
bool CTracker::prepare()
{
  m_bResult = true;
  m_bTerminate = false;

  if(!m_bReady && !read_data())
    return false;

  initial_conditions();

  if(m_pBarnesHut != NULL)
  {
    delete m_pBarnesHut;
    m_pBarnesHut = NULL;
  }

  return true;
}

UINT CTracker::main_thread_func(LPVOID pData)
{
  CalcThread* pThread = (CalcThread*)pData;
  CTracker* pObj = (CTracker*)pThread->m_pData;

  UINT nOutFreq;
  pObj->get_output_freq(nOutFreq);
  pObj->init_currents();

  const ULONG nStateSize = pObj->m_nType == CTrack::ptDroplet ? DROPLET_STATE_SIZE : ION_STATE_SIZE;
  STEPPER_TYPE_ID nIntegrType = GetIntegrType(pObj->get_integr_type());

  for(size_t i = pThread->get_first_job(); i <= pThread->get_last_job(); i++)
  {
    size_t nStep = 1;
    CTrack& track = pObj->m_Tracks.at(i);

    if(pObj->m_bOldIntegrator)
    {
      CBaseTrackItem* pItem = track.at(0)->copy();
      while(true)
      {
        if(pObj->m_bTerminate)
          break;

        bool bOK = pObj->m_nType == CTrack::ptDroplet ? 
          pObj->do_time_step(pItem) : pObj->do_ion_time_step(pItem, track.get_phase(), track.get_current());

        if(!bOK)
        {
          track.push_back(pItem->copy());  // insert the last item anyway.
          break;  // the track went out of the mesh.
        }

        if(nStep % nOutFreq == 0)
          track.push_back(pItem->copy());

        if(pObj->track_is_over(pItem->time, pItem->get_mass(), pItem->get_temp()))
          break;  // the integration time exceeded the limit or the Rayleigh criterion was met.

        nStep++;
      }

      delete pItem;
    }
    else
    {
      double fMass, fTemp;
      CBaseTrackItem* pItem = track.at(0);
      double fTime = pItem->time;
      double pState[max(ION_STATE_SIZE, DROPLET_STATE_SIZE)];
      pItem->state(pState);
      CIntegrInterface data(pItem->nElemId, track.get_index(), track.get_phase(), track.get_current());
      void* pI = create_integrator_interface(nStateSize, nIntegrType, (const void*)&data, CTracker::GetTimeDeriv);
      while(true)
      {
        do_integrator_step(pI, pState, &fTime, &pObj->m_fTimeStep);

        if(!data.bOk)
        {
          CBaseTrackItem* pCurrItem = CBaseTrackItem::create(pObj->m_nType, -1, fTime, pState);
          track.push_back(pCurrItem);  // insert the last item anyway.
          break;  // the track went out of the mesh.
        }

        if(nStep % nOutFreq == 0)
        {
          CBaseTrackItem* pCurrItem = CBaseTrackItem::create(pObj->m_nType, data.nElemId, fTime, pState);
          track.push_back(pCurrItem);
        }

        fTemp = pState[6];
        fMass = pObj->m_nType == CTrack::ptDroplet ? pState[7] : 0;
        if(pObj->track_is_over(fTime, fMass, fTemp))
          break;  // the integration time exceeded the limit or the Rayleigh criterion was met.

        fTime += pObj->m_fTimeStep;
        nStep++;
      }

      delete_integrator_interface(pI);
    }

    pThread->done_job();
    pObj->set_tracking_progress();
    if(pObj->m_bTerminate)
      break;
  }

  return 0;
}

void CTracker::multi_thread_calculate()
{
  if(!prepare())
  {
    m_bResult = false;
    return;
  }

  UINT nPartCount = m_Tracks.size();
  if(nPartCount == 0)
    return;

  bool bIter = (m_nType == CTrack::ptIon) && m_bEnableCoulomb && !m_bAxialSymm;
  UINT nIterCount = bIter ? m_nIterCount : 1;

  char buff[8];
  for(UINT i = 0; i < nIterCount; i++)
  {
    set_progress(0);
    CString sJobName(_T("Calculating particles")); 
    if(bIter)
      sJobName += CString(_T(", iteration # ")) + CString(itoa(i + 1, buff, 10));
    else
      sJobName += CString(_T("..."));

    set_job_name((const char*)sJobName);

    CalcThreadVector vThreads;
    m_pCalcThreads = &vThreads;
    vThreads.distribute_jobs(0, nPartCount - 1, &main_thread_func, (void*)this);
    vThreads.start_execution();
    vThreads.wait();

    if(m_bTerminate)
      return;

    if(bIter)
      capture_save_image(i);  // capture and save the screen image at every iteration.

    if(i == nIterCount - 1)
      break;
    if(i == 0)
      m_OutputEngine.prepare_current_output();

    if((nIterCount > 4) && (i >= nIterCount - 4))
      m_OutputEngine.add_current();    // average output currents over 4 latest iterations.

    create_BH_object(vThreads, i + 1);
    clear_tracks(false);    // keep the initial positions only.
  }

  relax();
}

void CTracker::relax()
{
  if(m_nType == CTrack::ptIon)
    m_OutputEngine.output_ion_current();

  if(!m_OutputEngine.get_enable_file_output())
    return;

  size_t nPartCount = m_Tracks.size();
  if(m_nType == CTrack::ptDroplet)
    for(size_t i = 0; i < nPartCount; i++)
      m_OutputEngine.output_droplet_track(i);
  else
    m_OutputEngine.output_ion_tracks();
}

void CTracker::set_tracking_progress()
{
  set_progress(m_pCalcThreads->get_progress());
}

//-------------------------------------------------------------------------------------------------
// Droplets type of particles.
//-------------------------------------------------------------------------------------------------
bool CTracker::do_time_step(CBaseTrackItem* pBaseItem)
{
  if(pBaseItem->nElemId < 0)
    return false;

  CDropletTrackItem* pItem = (CDropletTrackItem*)pBaseItem;

  Vector3D vAccel;
  bool bReflDone = sym_corr(pItem->pos, pItem->vel, vAccel);

  CNode3D node;
  CElem3D* pElem = m_vElems.at(pItem->nElemId);
// Predictor:
  if(!interpolate(pItem->pos, node, pElem))
  {
    if(bReflDone)
      sym_corr(pItem->pos, pItem->vel, vAccel, true);

    return false;
  }

  double fRe;
  double fD = get_particle_diameter(pItem->mass);
  vAccel = get_accel(node, pItem->vel, pItem->mass, fD, pItem->time, fRe);  // the Reynolds number is computed here.
// vAccel is the reflected particle acceleration. Turn back the position, velocity and acceleration and do a half time step:
  if(bReflDone)
    sym_corr(pItem->pos, pItem->vel, vAccel, true);

  double fCoolingRate = m_pEvaporModel->get_cooling_rate(node, pItem->temp, fD, fRe);
  double fEvaporRate = m_pEvaporModel->get_evaporation_rate(node, pItem->temp, fD, fRe);

  double fTemp05 = pItem->temp - fCoolingRate * m_fHalfTimeStep;
  double fMass05 = pItem->mass - fEvaporRate * m_fHalfTimeStep;
  if(fMass05 < m_fMinMass)
    return false;

  Vector3D vVel05 = pItem->vel + vAccel * m_fHalfTimeStep;
  Vector3D vPos05 = pItem->pos + (pItem->vel + 0.5 * vAccel * m_fHalfTimeStep) * m_fHalfTimeStep;

// Corrector: getting more precise estimates of average acceleration during all dt.
// If vPos05 is out of the gas-dynamic domain, reflect both position and velocity (acceleration is a dummy parameter here)
// and remember the result of this action:
  bReflDone = sym_corr(vPos05, vVel05, vAccel);

  if(!interpolate(vPos05, node, pElem))
    return false;

  double fRe05;
  double fD05 = get_particle_diameter(fMass05);
  double fTime05 = pItem->time + m_fHalfTimeStep;
  vAccel = get_accel(node, vVel05, fMass05, fD05, fTime05, fRe05); // the Reynolds number is computed here.
// vAccel is the reflected particle acceleration. Turn back the velocity and acceleration and do the full time step:
  if(bReflDone)
    sym_corr(vPos05, vVel05, vAccel, true); // vPos05 is a dummy parameter here.

  fCoolingRate = m_pEvaporModel->get_cooling_rate(node, fTemp05, fD05, fRe05);
  fEvaporRate = m_pEvaporModel->get_evaporation_rate(node, fTemp05, fD05, fRe05);

  pItem->temp -= fCoolingRate * m_fTimeStep;
  pItem->mass -= fEvaporRate * m_fTimeStep;
  if(pItem->mass < m_fMinMass)
    return false;

  pItem->nElemId = pElem->nInd;
  pItem->pos += vVel05 * m_fTimeStep;
  pItem->vel += vAccel * m_fTimeStep;
  pItem->time += m_fTimeStep;
  return true;
}

Vector3D CTracker::get_accel(const CNode3D& node, const Vector3D& vVel, double fMass, double fD, double fTime, double& fRe) const
{
// Dynamics of a body with a variable mass: dV/dt = F/m + (U/m)dm/dt, where U is a relative velocity of the
// detached mass dm. In the case of the evaporating particle U = 0, so that the 2-nd Newton's law will not change.
  Vector3D accel(0, 0, 0);
// Electrostatics:
  if(m_bEnableField || m_bEnableRF)
  {
    double fEoverM = m_fCharge / fMass;

    if(m_bEnableField)
      accel += fEoverM * (node.field * m_fAmplDC + m_vFieldPtbColl.apply(node.pos));

    if(m_bEnableRF)
    {
      double fAmpl, fOmega;
      get_ampl_freq(node.pos, fAmpl, fOmega);
      accel += fEoverM * node.rf * fAmpl * sin(fOmega * fTime);
    }
  }
// Gas-dynamics:
  Vector3D vDiff = node.vel - vVel;
  fRe = get_Re(vDiff, node.dens, node.visc, fD);  // the character length fD varies with time.
  double fCd = get_Cd(fRe);
  accel += 0.75 * (node.dens / m_fPartDens) * (fCd / fD) * vDiff.length() * vDiff;

  return accel;
}

//-------------------------------------------------------------------------------------------------
// Track limitations.
//-------------------------------------------------------------------------------------------------
bool CTracker::track_is_over(double fTime, double fMass, double fTemp) const
{
  if(fTime >= m_fMaxIntegrTime)
    return true;

  if(m_nType == CTrack::ptIon)
    return false;

  double fD = get_particle_diameter(fMass);

  return fD <= m_fMinD || limit_of_Rayleigh(fTemp, fD);
}

double CTracker::get_max_charge(double fT, double fD) const
{
  if(fT < 273)
    return FLT_MAX;

  double fSigma = m_pEvaporModel->get_surface_tension(fT);
  return sqrt(Const_2PI * fSigma * fD * fD * fD);
}

bool CTracker::limit_of_Rayleigh(double fT, double fD) const
{
  if(fT < 273 || !m_pEvaporModel->get_enable_surf_tens())
    return false; // the Rayleigh criterion does not work at temperatures lower than the freezing point.

  double fLimCharge = get_max_charge(fT, fD);
  if(m_fCharge >= fLimCharge)
    return true;

  return false;
}

//-------------------------------------------------------------------------------------------------
// Ion type of particles.
//-------------------------------------------------------------------------------------------------
bool CTracker::do_ion_time_step(CBaseTrackItem* pBaseItem, double fPhase, double fCurr)
{
  if(pBaseItem->nElemId < 0)
    return false;

  CIonTrackItem* pItem = (CIonTrackItem*)pBaseItem;

  Vector3D vAccel;
// After this call item.pos and item.vel are in the gas-dynamic domain. If the reflection is really
// done, the bReflDone flag will be "true". Note that vAccel is a dummy parameter here.
  bool bReflDone = sym_corr(pItem->pos, pItem->vel, vAccel);

  CNode3D node;
  CElem3D* pElem = m_vElems.at(pItem->nElemId);
// Predictor:
// Compute all ANSYS fields at item.pos and place them into the node structure.
  if(!interpolate(pItem->pos, node, pElem))
  {
    if(bReflDone) // if interpolate(...) fails set the reflected vectors back and terminate the integration.
      sym_corr(pItem->pos, pItem->vel, vAccel, true);

    return false;
  }

  double fExpCoeff;
  vAccel = get_ion_accel(node, pItem->vel, pItem->time, m_fHalfTimeStep, fPhase, fCurr, fExpCoeff);
// vAccel is the reflected particle acceleration. Turn back the position, velocity and acceleration and do a half time step:
  if(bReflDone)
    sym_corr(pItem->pos, pItem->vel, vAccel, true);

  Vector3D vVel05 = pItem->vel + vAccel * m_fHalfTimeStep;
  Vector3D vPos05 = pItem->pos + (pItem->vel + 0.5 * vAccel * m_fHalfTimeStep) * m_fHalfTimeStep;

  double fTionInf;
  double fTi05 = pItem->temp + get_dTi(node, pItem->vel, pItem->temp, fExpCoeff, fTionInf) * m_fHalfTimeStep;

// Corrector: getting more precise estimates of average acceleration during all dt.
// If vPos05 is out of the gas-dynamic domain, reflect both position and velocity (acceleration is a dummy parameter here)
// and remember the result of this action:
  bReflDone = sym_corr(vPos05, vVel05, vAccel);

  if(!interpolate(vPos05, node, pElem))
    return false;

  double fTime05 = pItem->time + m_fHalfTimeStep;
  vAccel = get_ion_accel(node, vVel05, fTime05, m_fTimeStep, fPhase, fCurr, fExpCoeff);
// vAccel is the reflected particle acceleration. Turn back the velocity and acceleration and do the full time step:
  if(bReflDone)
    sym_corr(vPos05, vVel05, vAccel, true); // vPos05 is a dummy parameter here.

  pItem->nElemId = pElem->nInd;
// Correct position and velocity. Whether or not they are in the gas-dynamic domain will be checked at the next time step.
  pItem->pos += vVel05 * m_fTimeStep;
  pItem->vel += vAccel * m_fTimeStep;
  pItem->temp += get_dTi(node, vVel05, fTi05, fExpCoeff, fTionInf) * m_fTimeStep;
  pItem->tempinf = fTionInf;
  pItem->time += m_fTimeStep;
  return true;
}

Vector3D CTracker::get_ion_accel(const CNode3D&  node, 
                                 const Vector3D& vVel,
                                 double          fTime,
                                 double          fTimeStep,
                                 double          fPhase,
                                 double          fCurr,
                                 double&         fExpCoeff) const
{
// Entraining by the gas:
  double fMob = get_ion_mob(node.press, node.temp);
  double fOneOvrTau = m_fChargeMassRatio / fMob;  // 1 / tau.
  double fPow = fTimeStep * fOneOvrTau;           // in fact, fPow = dt / tau.
  fExpCoeff = fPow < 0.01 ? fOneOvrTau : (1. - exp(-fPow)) / fTimeStep;

  Vector3D vE(0, 0, 0);
// DC Field:
  if(m_bEnableField)
    vE += (node.field * m_fAmplDC + m_vFieldPtbColl.apply(node.pos));

// RF field:
  if(m_bEnableRF)
    vE += get_rf_field(node, fTime, fPhase);

// Coulomb:
  if(m_bEnableCoulomb)
  {
    if(m_bAxialSymm && (vVel.x > Const_Almost_Zero))
    {
      vE += CSpaceChargeDistrib::radial_coulomb_field(node.pos, vVel.x, fCurr);
    }
    else if(m_pBarnesHut != NULL) // in the case of iterations m_pBarnesHut is non-zero since the second iteration.
    {
      Vector3D vClmb(0, 0, 0);
      if(m_bUseRadialCoulomb && (node.pos.x > m_fRadialCoulombX) && (vVel.x > Const_Almost_Zero))
      {
        vClmb = m_SpaceChargeDist.radial_coulomb(node.pos, vVel.x);
      }
      else
      {
        vClmb = m_pBarnesHut->coulomb_force(node.pos);
// DEBUG
        if((node.pos.x > m_fX_Q0) && (vClmb.x < 0))
          vClmb.x = 0;
// END DEBUG
      }

      vE += vClmb;
    }
  }

  return (node.vel - vVel + fMob * vE) * fExpCoeff;
}

double CTracker::get_dTi(const CNode3D& node, const Vector3D& vVel, double fIonTemp, double fExpCoeff, double& fTinf) const
{
  fTinf = node.temp + 0.2 * m_fMolarMass * (node.vel - vVel).sqlength() / Const_Boltzmann;
  return (fTinf - fIonTemp) * fExpCoeff;
}

Vector3D CTracker::get_rf_field(const CNode3D& node, double fTime, double fPhase) const
{
  if(node.pos.x < m_fX_Q00) // in the funnel (S-Lens).
    return node.rf * (m_fAmplRF * sin(m_fOmega * fTime + fPhase));

  const double cfHalfWidth = 0.4; // cm, the half width of the transition zone.
  double fLowX = m_fX_Q0 - cfHalfWidth;
  double fHighX = m_fX_Q0 + cfHalfWidth;

  if(node.pos.x < fLowX)
  {
    return node.rf * (m_fAmplRF_Q00 * sin(m_fOmega_Q00 * fTime + fPhase));  // purely Q00 field.
  }
  else if(node.pos.x > fHighX)
  {
// DEBUG
    double fR0 = 0.21;  // 2.1 mm - an inscribed radius in the flatapole.
    double fCoeff = 2. * SI_to_CGS_Voltage / (fR0 * fR0);
    Vector3D vRF(0, -fCoeff * node.pos.y, fCoeff * node.pos.z);
/*
    Vector3D vRF = node.rf;
    vRF.x = 0;  // artificial move to get rid of the mirroring effect in the Q0 region.
*/
// END DEBUG
    return vRF * (m_fAmplRF_Q0 * sin(m_fOmega_Q0 * fTime + fPhase));  // purely Q0 field.
  }

  double fKsi = (node.pos.x - fLowX) / (fHighX - fLowX);
  double fSmoothStep = fKsi * fKsi * (3 - 2 * fKsi);

  double fAQ00 = (1 - fSmoothStep) * m_fAmplRF_Q00;
  double fAQ0 = fSmoothStep * m_fAmplRF_Q0;

  return node.rf * (fAQ00 * sin(m_fOmega_Q00 * fTime + fPhase) + fAQ0 * sin(m_fOmega_Q0 * fTime + fPhase));
}

void CTracker::init_currents()
{
  if(m_nType != CTrack::ptIon || !m_bEnableCoulomb || m_fInitBunchRadius < Const_Almost_Zero)
    return;

  size_t nIonCount = m_Tracks.size();
  for(size_t i = 0; i < nIonCount; i++)
  {
    CTrack& track = m_Tracks.at(i);
    if(track.size() == 0)
      continue;

    CBaseTrackItem* pItem = track.at(0);
    double fR2 = pItem->pos.y * pItem->pos.y + pItem->pos.z * pItem->pos.z;
    double fR02 = m_fInitBunchRadius * m_fInitBunchRadius;
    if(m_bAxialSymm)
      track.set_current(m_fFullCurrent * fR2 / fR02);
  }
}

//-------------------------------------------------------------------------------------------------
// Reflection of positions and velocities due to symmetry. If any of XY or XZ symmetry planes exist
// POSITIVE z and y are supposed to be inside the domain.
//-------------------------------------------------------------------------------------------------
bool CTracker::sym_corr(Vector3D& vPos, Vector3D& vVel, Vector3D& vAccel, bool bForceReflect) const
{
  bool bReflDone = false;
  if((m_nSymPlanes & spXY) && (vPos.z < 0 || bForceReflect))
  {
    vPos.z = -vPos.z;
    vVel.z = -vVel.z;
    vAccel.z = -vAccel.z;
    bReflDone = true;
  }

  if((m_nSymPlanes & spXZ) && (vPos.y < 0 || bForceReflect))
  {
    vPos.y = -vPos.y;
    vVel.y = -vVel.y;
    vAccel.y = -vAccel.y;
    bReflDone = true;
  }

  if((m_nSymPlanes & spYZ) && (vPos.x < 0 || bForceReflect))
  {
    vPos.x = -vPos.x;
    vVel.x = -vVel.x;
    vAccel.x = -vAccel.x;
    bReflDone = true;
  }

  return bReflDone;
}

typedef math::Reflector<Vector3D, double, CTracker> CSymmReflect;
//-------------------------------------------------------------------------------------------------
// Initial conditions
//-------------------------------------------------------------------------------------------------
void CTracker::initial_conditions()
{
  set_job_name("Initial conditions...");
  set_progress(0);

  clear_tracks(); // clear all items from the tracks, even the initial positions.

  Vector3D vPos, vVel;
  double fTime, fPhase, fTemp;
  UINT nEnsIndex;
  size_t nElemId;

  if(!m_Src.generate_initial_cond())
    return;

  UINT nCount = m_Src.get_particles_count() * m_Src.get_ensemble_size();  // maximal expected count of particles.
  for(size_t i = 0; i < nCount; i++)
  {
    m_Src.get(i, vPos, vVel, fTime, fPhase, fTemp, nEnsIndex, nElemId);  // this function has a built-in index check.
    
    CTrack track(m_nType, nEnsIndex, fPhase);
    track.reserve(100);

    CBaseTrackItem* pItem = create_track_item(nElemId, vPos, vVel, m_fInitMass, fTemp, fTime);
    track.push_back(pItem);

    m_Tracks.push_back(track);
  }
}

CBaseTrackItem* CTracker::create_track_item(size_t          nElemId,
                                            const Vector3D& vPos,
                                            const Vector3D& vVel,
                                            double          fMass,
                                            double          fTemp,
                                            double          fTime) const
{
  if(m_nType == CTrack::ptIon)
  {
    CIonTrackItem* pIonItem = new CIonTrackItem(nElemId, vPos, vVel, fTemp, 1., fTime);
    return (CBaseTrackItem*)pIonItem;
  }
  else if(m_nType == CTrack::ptDroplet)
  {
    CDropletTrackItem* pItem = new CDropletTrackItem(nElemId, vPos, vVel, fTemp, fMass, fTime);
    return (CBaseTrackItem*)pItem;
  }

  return NULL;
}

void CTracker::clear_tracks(bool bFinally)
{
  size_t nTrackSize = m_Tracks.size();
  for(size_t i = 0; i < nTrackSize; i++)
  {
#ifndef _DEBUG
    if(bFinally && (i % 20 == 0))
      set_status("Clearing tracks", 100 * i / nTrackSize);
#endif
    CTrack& track = m_Tracks.at(i);
    size_t nFirst = bFinally ? 0 : 1;
    size_t nLast = track.size();
    for(size_t j = nFirst; j < nLast; j++)
      delete track.at(j);

    if(track.size() > 1)
      track.erase(track.begin() + nFirst, track.end());
  }

  if(bFinally)
    m_Tracks.clear();
#ifndef _DEBUG
  set_status("Ready", -1);
#endif
}

//-------------------------------------------------------------------------------------------------
// Mesh specific interface (data reading, elements finding, interpolation)
//-------------------------------------------------------------------------------------------------
bool CTracker::interpolate(const Vector3D& vPos, CNode3D& node, CElem3D*& pElem) const
{
  pElem = find_elem(pElem, vPos);
  if(pElem == NULL)
    return false;

  pElem->interpolate(vPos, node);

  return true;
}

CElem3D* CTracker::find_elem(CElem3D* pPrevElem, const Vector3D& vPos) const
{
  CElem3D* pElem = NULL;
// First, try the nearest neighbors of the previous face (if it is not zero), including the previous face itself:
  if(pPrevElem != NULL)
  {
    pElem = try_neighbors(pPrevElem, vPos);
    if(pElem != NULL)
      return pElem;
  }

  size_t nCount = m_vElems.size();  // run over all the elements collection:
  for(size_t i = 0; i < nCount; i++)
  {
    pElem = m_vElems.at(i);
    if(pElem->inside(vPos))
      return pElem;
  }

  return NULL;
}

CElem3D* CTracker::try_neighbors(CElem3D* pElem, const Vector3D& vPos) const
{
// First, try the input element itself.
  if(pElem->inside(vPos))
    return pElem;

// Try the neighbors of pElem.
  std::vector<CElem3D*> vUsedElems;
  size_t nNodeCount = pElem->vNodes.size();
  for(size_t j = 0; j < nNodeCount; j++)
  {
    CNode3D* pNode = pElem->vNodes.at(j);
    size_t nNbrCount = pNode->vNbrElems.size();
    for(size_t i = 0; i < nNbrCount; i++)
    {
      CElem3D* pNbrElem = pNode->vNbrElems.at(i);
      if(std::find(vUsedElems.begin(), vUsedElems.end(), pNbrElem) != vUsedElems.end())
        continue;

      if(pNbrElem->inside(vPos))
        return pNbrElem;

      vUsedElems.push_back(pNbrElem);
    }
  }

  return NULL;
}

bool CTracker::abort(FILE* pStream)
{
  if(pStream != NULL)
    fclose(pStream);

  m_bReady = false;

  return false;
}

bool CTracker::read_data()
{
  clear();

  if(!read_geometry())
    return false;

  if(!read_gasdyn_data()) // read the gas-dynamic variables in the nodes of the mesh.
    return false; // wrong gas-dynamic data format.

  if(!read_2D_regions())
    return false;

  invalidate_calculators();
  m_bReady = true;
  return true;
}

bool CTracker::read_geometry()
{
  set_job_name("Reading points...");
  set_progress(0);

  std::string cBase = COutputEngine::get_base_name(get_filename());

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (cBase + "geom").c_str(), (const char*)("r"));
  if(nErr != 0 || pStream == 0)
    return false;

  char cStr[32];
  int nNode, nNodeCount, nRes = 0;
  float fX, fY, fZ;

  nRes = fscanf_s(pStream, "%s", cStr, 32);
  nRes = fscanf_s(pStream, "%d", &nNodeCount);
// Nodes:
  m_vNodes.clear();
  m_vNodes.reserve(nNodeCount);

  CNode3D* pNode = NULL;
// Important! Nodes in ANSYS output are enumerated from 1 to nNodeCount. In connectivity data, too!
// In our arrays we enumerate nodes from 0 to nNodeCount - 1, i.e. subtract unity from the original node index.
  for(UINT i = 0; i < nNodeCount; i++)
  {
    nRes = fscanf_s(pStream, "%d %e %e %e", &nNode, &fX, &fY, &fZ);
    if(nRes == EOF)
      return false;

    fX *= SI_to_CGS_Len;
    fY *= SI_to_CGS_Len;
    fZ *= SI_to_CGS_Len;

    Vector3D vPos(fX, fY, fZ);
// Mesh transformation:
    if(m_Transform.get_enable())
      m_Transform.transform(vPos);

    pNode = new CNode3D(vPos);
    pNode->nInd = i;
    m_vNodes.push_back(pNode);

// Progress bar and termination support:
    if(i % 100 == 0)
      set_progress(int(0.5 + 100. * i / nNodeCount));
    if(m_bTerminate)
      return abort(pStream);
  }

// Connectivity information.
  m_vElems.clear();

// Tetrahedra (4 nodes, 4 planes):
  int nTetraCount, n0, n1, n2, n3;
  nRes = fscanf_s(pStream, "%s", cStr, 32);
  nRes = fscanf_s(pStream, "%d", &nTetraCount);
  m_vElems.reserve(nTetraCount);

  if(nTetraCount > 0)
  {
    set_job_name("Adding tetra...");

    for(UINT i = 0; i < nTetraCount; i++)
    {
      nRes = fscanf_s(pStream, "%d %d %d %d", &n0, &n1, &n2, &n3);
      if(nRes == EOF)
        return false;

      add_tetra(m_vNodes.at(n0 - 1), m_vNodes.at(n1 - 1), m_vNodes.at(n2 - 1), m_vNodes.at(n3 - 1));

// Progress bar and termination support:
      if(i % 100 == 0)
        set_progress(int(0.5 + 100. * i / nTetraCount));
      if(m_bTerminate)
        return abort(pStream);
    }
  }

// Pyramids (5 nodes, 5 planes):
  int nPyrCount, n4;
  nRes = fscanf_s(pStream, "%s", cStr, 32);
  nRes = fscanf_s(pStream, "%d", &nPyrCount);
  m_vElems.reserve(nTetraCount + nPyrCount);

  if(nPyrCount > 0)
  {
    set_job_name("Adding pyramids...");

    for(UINT i = 0; i < nPyrCount; i++)
    {
      nRes = fscanf_s(pStream, "%d %d %d %d %d", &n0, &n1, &n2, &n3, &n4);
      if(nRes == EOF)
        return false;

      add_pyramid(m_vNodes.at(n0 - 1), m_vNodes.at(n1 - 1), m_vNodes.at(n2 - 1), m_vNodes.at(n3 - 1), m_vNodes.at(n4 - 1));

// Progress bar and termination support:
      if(i % 100 == 0)
        set_progress(int(0.5 + 100. * i / nPyrCount));
      if(m_bTerminate)
        return abort(pStream);
    }
  }

// Wedges (or prisms, 6 nodes, 5 planes):
  int nWedgeCount, n5;
  nRes = fscanf_s(pStream, "%s", cStr, 32);
  nRes = fscanf_s(pStream, "%d", &nWedgeCount);
  m_vElems.reserve(nTetraCount + nPyrCount + nWedgeCount);

  if(nWedgeCount > 0)
  {
    set_job_name("Adding wedges...");

    for(UINT i = 0; i < nWedgeCount; i++)
    {
      nRes = fscanf_s(pStream, "%d %d %d %d %d %d", &n0, &n1, &n2, &n3, &n4, &n5);
      if(nRes == EOF)
        return false;

      add_wedge(m_vNodes.at(n0 - 1), m_vNodes.at(n1 - 1), m_vNodes.at(n2 - 1), m_vNodes.at(n3 - 1), m_vNodes.at(n4 - 1), m_vNodes.at(n5 - 1));

// Progress bar and termination support:
      if(i % 100 == 0)
        set_progress(int(0.5 + 100. * i / nWedgeCount));
      if(m_bTerminate)
        return abort(pStream);
    }
  }

// Hexahedra (8 nodes, 6 planes):
  int nHexaCount, n6, n7;
  nRes = fscanf_s(pStream, "%s", cStr, 32);
  nRes = fscanf_s(pStream, "%d", &nHexaCount);
  m_vElems.reserve(nTetraCount + nPyrCount + nWedgeCount + nHexaCount);

  if(nHexaCount > 0)
  {
    set_job_name("Adding hexa...");

    for(UINT i = 0; i < nHexaCount; i++)
    {
      nRes = fscanf_s(pStream, "%d %d %d %d %d %d %d %d", &n0, &n1, &n2, &n3, &n4, &n5, &n6, &n7);
      if(nRes == EOF)
        return false;

      add_hexa(m_vNodes.at(n0 - 1), m_vNodes.at(n1 - 1), m_vNodes.at(n3 - 1), m_vNodes.at(n2 - 1), 
               m_vNodes.at(n4 - 1), m_vNodes.at(n5 - 1), m_vNodes.at(n7 - 1), m_vNodes.at(n6 - 1));

// Progress bar and termination support:
      if(i % 100 == 0)
        set_progress(int(0.5 + 100. * i / nHexaCount));
      if(m_bTerminate)
        return abort(pStream);
    }
  }

  fclose(pStream);
  bounding_box();
  return true;
}

bool CTracker::read_2D_regions()
{
  set_job_name("Reading 2D regions...");
  set_progress(0);

  std::string cBase = COutputEngine::get_base_name(get_filename());

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (cBase + "rgn").c_str(), (const char*)("r"));
  if(nErr != 0 || pStream == 0)
    return false;

  char cName[32];
  int nRes, n0, n1, n2, n3;
  UINT nRegCount, nTriCount, nQuadCount, i, j, nFace;
  CFace* pFace = NULL;

  nRes = fscanf_s(pStream, "%d", &nRegCount);
  for(i = 0; i < nRegCount; i++)
  {
    nRes = fscanf(pStream, "%s", cName);
    CRegion* pReg = new CRegion(cName);

    nRes = fscanf_s(pStream, "%d", &nTriCount);
    for(j = 0; j < nTriCount; j++)
    {
      nRes = fscanf_s(pStream, "%d %d %d", &n0, &n1, &n2);
      n0--;
      n1--;
      n2--;
      pFace = new CFace(m_vNodes.at(n0), m_vNodes.at(n1), m_vNodes.at(n2));
      pReg->vFaces.push_back(pFace);
      nFace = pReg->vFaces.size() - 1;
      CRegFacePair pair(i, nFace);
      m_vNodes.at(n0)->vNbrFaces.push_back(pair);
      m_vNodes.at(n1)->vNbrFaces.push_back(pair);
      m_vNodes.at(n2)->vNbrFaces.push_back(pair);
    }

    nRes = fscanf_s(pStream, "%d", &nQuadCount);
    for(j = 0; j < nQuadCount; j++)
    {
      nRes = fscanf_s(pStream, "%d %d %d %d", &n0, &n1, &n2, &n3);
      n0--;
      n1--;
      n2--;
      n3--;
      pFace = new CFace(m_vNodes.at(n0), m_vNodes.at(n1), m_vNodes.at(n2));
      pReg->vFaces.push_back(pFace);
      nFace = pReg->vFaces.size() - 1;
      CRegFacePair pair1(i, nFace);
      m_vNodes.at(n0)->vNbrFaces.push_back(pair1);
      m_vNodes.at(n1)->vNbrFaces.push_back(pair1);
      m_vNodes.at(n2)->vNbrFaces.push_back(pair1);

      pFace = new CFace(m_vNodes.at(n0), m_vNodes.at(n2), m_vNodes.at(n3));
      pReg->vFaces.push_back(pFace);
      nFace = pReg->vFaces.size() - 1;
      CRegFacePair pair2(i, nFace);
      m_vNodes.at(n0)->vNbrFaces.push_back(pair2);
      m_vNodes.at(n2)->vNbrFaces.push_back(pair2);
      m_vNodes.at(n3)->vNbrFaces.push_back(pair2);
    }

    pReg->bounding_box();
    m_vRegions.push_back(pReg);

// Progress bar and termination support:
    if(i % 100 == 0)
      set_progress(int(0.5 + 100. * i / nRegCount));
    if(m_bTerminate)
      return abort(pStream);
  }

  fclose(pStream);
  return true;
}

bool CTracker::read_gasdyn_data()
{
  set_job_name("Reading gas-dynamic data...");
  set_progress(0);

  std::string cBase = COutputEngine::get_base_name(get_filename());;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (cBase + "var").c_str(), (const char*)("r"));
  if(nErr != 0 || pStream == 0)
    return false;

  int nRes = 0;
  char sHeader[16];
  for(UINT j = 0; j < 24; j++)
    nRes = fscanf(pStream, "%s", sHeader);

  float fPress, fDens, fTemp, fDynVisc, fThermCond, fCp, fVx, fVy, fVz, fEx, fEy, fEz, fRFEx, fRFEy, fRFEz;

  CNode3D* pNode = NULL;
  Vector3D vVel, vFieldDC, vFieldRF;
  size_t nNodeCount = m_vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    nRes = fscanf_s(pStream, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",
      &fPress, &fDens, &fTemp, &fDynVisc, &fThermCond, &fCp, &fVx, &fVy, &fVz, &fEx, &fEy, &fEz, &fRFEx, &fRFEy, &fRFEz);

    if(nRes == EOF)
      return false;

    if(m_bConv2CGS)
      conv_to_cgs(fPress, fDens, fDynVisc, fThermCond, fCp, fVx, fVy, fVz, fEx, fEy, fEz, fRFEx, fRFEy, fRFEz);

    vVel = Vector3D(fVx, fVy, fVz);
    vFieldDC = Vector3D(fEx, fEy, fEz);
    vFieldRF = Vector3D(fRFEx, fRFEy, fRFEz);
// Mesh transformation:
    if(m_Transform.get_enable())
    {
      m_Transform.transform(vVel, false);
      m_Transform.transform(vFieldDC, false);
      m_Transform.transform(vFieldRF, false);
    }

    pNode = m_vNodes.at(i);
    pNode->set_data(fPress, fDens, fTemp, fDynVisc, fThermCond, fCp, vVel, vFieldDC, vFieldRF);

// Progress bar and termination support:
    if(i % 100 == 0)
      set_progress(int(0.5 + 100. * i / nNodeCount));
    if(m_bTerminate)
      return abort(pStream);
  }

  fclose(pStream);
  return true;
}

void CTracker::add_tetra(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3)
{
  CTetra* pTetra = new CTetra(p0, p1, p2, p3);
  pTetra->nInd = m_vElems.size(); // this index is used in integrators as well as in export of the mesh to the OpenFOAM format.
  m_vElems.push_back((CElem3D*)pTetra);

  p0->vNbrElems.push_back((CElem3D*)pTetra);
  p1->vNbrElems.push_back((CElem3D*)pTetra);
  p2->vNbrElems.push_back((CElem3D*)pTetra);
  p3->vNbrElems.push_back((CElem3D*)pTetra);
}

void CTracker::add_pyramid(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4)
{
  CPyramid* pPyr = new CPyramid(p0, p1, p2, p3, p4);
  pPyr->nInd = m_vElems.size();
  m_vElems.push_back((CElem3D*)pPyr);

  p0->vNbrElems.push_back((CElem3D*)pPyr);
  p1->vNbrElems.push_back((CElem3D*)pPyr);
  p2->vNbrElems.push_back((CElem3D*)pPyr);
  p3->vNbrElems.push_back((CElem3D*)pPyr);
  p4->vNbrElems.push_back((CElem3D*)pPyr);
}

void CTracker::add_wedge(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4, CNode3D* p5)
{
  CWedge* pWedge = new CWedge(p0, p1, p2, p3, p4, p5);
  pWedge->nInd = m_vElems.size();
  m_vElems.push_back((CElem3D*)pWedge);

  p0->vNbrElems.push_back((CElem3D*)pWedge);
  p1->vNbrElems.push_back((CElem3D*)pWedge);
  p2->vNbrElems.push_back((CElem3D*)pWedge);
  p3->vNbrElems.push_back((CElem3D*)pWedge);
  p4->vNbrElems.push_back((CElem3D*)pWedge);
  p5->vNbrElems.push_back((CElem3D*)pWedge);
}

void CTracker::add_hexa(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4, CNode3D* p5, CNode3D* p6, CNode3D* p7)
{
  CHexa* pHexa = new CHexa(p0, p1, p2, p3, p4, p5, p6, p7);
  pHexa->nInd = m_vElems.size();
  m_vElems.push_back((CElem3D*)pHexa);

  p0->vNbrElems.push_back((CElem3D*)pHexa);
  p1->vNbrElems.push_back((CElem3D*)pHexa);
  p2->vNbrElems.push_back((CElem3D*)pHexa);
  p3->vNbrElems.push_back((CElem3D*)pHexa);
  p4->vNbrElems.push_back((CElem3D*)pHexa);
  p5->vNbrElems.push_back((CElem3D*)pHexa);
  p6->vNbrElems.push_back((CElem3D*)pHexa);
  p7->vNbrElems.push_back((CElem3D*)pHexa);
}

void CTracker::bounding_box()
{
  size_t nNodesCount = m_vNodes.size();
  if(nNodesCount == 0)
    return;

  m_Box.vMin = m_vNodes.at(0)->pos;
  m_Box.vMax = m_vNodes.at(0)->pos;
  
  for(size_t i = 1; i < nNodesCount; i++)
  {
    CNode3D* pNode = m_vNodes.at(i);

    if(pNode->pos.x < m_Box.vMin.x)
      m_Box.vMin.x = pNode->pos.x;
    if(pNode->pos.x > m_Box.vMax.x)
      m_Box.vMax.x = pNode->pos.x;

    if(pNode->pos.y < m_Box.vMin.y)
      m_Box.vMin.y = pNode->pos.y;
    if(pNode->pos.y > m_Box.vMax.y)
      m_Box.vMax.y = pNode->pos.y;

    if(pNode->pos.z < m_Box.vMin.z)
      m_Box.vMin.z = pNode->pos.z;
    if(pNode->pos.z > m_Box.vMax.z)
      m_Box.vMax.z = pNode->pos.z;
  }
}

double CTracker::get_full_current_at(UINT nIter)
{
  const UINT nMaxCurrentIterCount = 5;
  const UINT nVarIterCount = m_nIterCount > nMaxCurrentIterCount ? m_nIterCount - nMaxCurrentIterCount : 0;
  if(nVarIterCount == 0)
    return m_fFullCurrent;

  if(nIter > nVarIterCount)
    return m_fFullCurrent;

  double fKsi = double(nIter) / nVarIterCount;
  return m_fFullCurrent * fKsi * fKsi * (3 - 2 * fKsi);
}

//-------------------------------------------------------------------------------------------------
// Evaporation:
//-------------------------------------------------------------------------------------------------
void CTracker::create_evapor_model()
{
  if(m_pEvaporModel != NULL)
    delete m_pEvaporModel;

  switch(m_nEvaporModelType)
  {
    case emNone: m_pEvaporModel = new CEvaporationModel(); break;
    case emMaxwell: m_pEvaporModel = new CMaxwellModel(); break;
    case emSteadyDiffusive: m_pEvaporModel = new CSteadyDiffusiveModel(); break;
    case emDiffusive: m_pEvaporModel = new CDiffusiveModel(); break;
  }
}

//-------------------------------------------------------------------------------------------------
// Streaming:
//-------------------------------------------------------------------------------------------------
void CTracker::save(CArchive& ar)
{
  set_data(); // data are copied from the properties list to the model.

  const UINT nVersion = 22;  // 21 - saving fields; 20 - m_vFieldPtbColl; 19 - m_Transform; 16 - m_nIntegrType; 15 - saving tracks; 14 - m_OutputEngine; 11 - Coulomb for non-axial cases; 10 - Calculators; 9 - RF in flatapole; 8 - m_bVelDependent; 7 - m_bByInitRadii and m_nEnsByRadiusCount; 6 - Data Importer; 5 - Coulomb effect parameters; 4 - m_bOnlyPassedQ00 and m_fActEnergy; 3 - the export OpenFOAM object is saved since this version.
  ar << nVersion;

// General parameters:
  CString cFileName(m_sDataFile.c_str());
  ar << cFileName;
  CString cFieldFileName(m_sFieldDataFile.c_str());
  ar << cFieldFileName;   // since version 1.

  ar << m_bUseMultiThread;
  ar << m_nType;
  ar << m_fTimeStep;
  ar << m_fMaxIntegrTime;
  ar << m_fInitD;
  ar << m_nSymPlanes;

  m_Transform.save(ar);

// Electrostatics:
  ar << m_bEnableField;
  ar << m_fAmplDC;
  ar << m_fCharge;

// Ion parameters:
  ar << m_fIonMass;
  ar << m_fIonMobility;
// RF data:
  ar << m_bEnableRF;
  ar << m_fAmplRF;
  ar << m_fOmega;
  ar << m_fAmplRF_Q00;
  ar << m_fOmega_Q00;
  ar << m_fX_Q00;

// Particle source:
  m_Src.save(ar);

// Evaporation:
  ar << m_nEvaporModelType;
  if(m_pEvaporModel != NULL)
    m_pEvaporModel->save(ar);

  EvaporatingParticle::CExportOpenFOAM* pExpObj = CParticleTrackingApp::Get()->GetExporter();
  pExpObj->save(ar);

  ar << m_fActEnergy;

  ar << m_bEnableCoulomb;
  ar << m_fInitBunchRadius;
  ar << m_fFullCurrent;
// Since version 11:
  ar << m_bAxialSymm;
  ar << m_nIterCount;
//  ar << m_nPlanesCount;
  ar << m_fBHDist;
// Since version 12:
  ar << m_fBHCritR;
  ar << m_nMaxRecDepth;

// Data importer since version 6:
  m_Importer.save(ar);

  ar << m_bVelDependent;

// RF in flatapole (Q0 region) since version 9:
  ar << m_fAmplRF_Q0;
  ar << m_fOmega_Q0;
  ar << m_fX_Q0;

// Calculators:
  CCalcCollection* pCalcColl = CParticleTrackingApp::Get()->GetCalcs();
  pCalcColl->save(ar);

  m_OutputEngine.save(ar);  // since version 14.

  ar << m_bSaveTracks;
  if(m_bSaveTracks)
    save_tracks(ar);

  int nIntegrType = (int)m_nIntegrType;
  ar << nIntegrType;

  ar << m_bOldIntegrator;
  ar << m_bUseRadialCoulomb;
  ar << m_fRadialCoulombX;

  m_SpaceChargeDist.save(ar);
  m_vFieldPtbColl.save(ar); // since version 20.

  CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
  pFields->save(ar);  // since version 21.

  CCrossSectColl* pCollCS = CParticleTrackingApp::Get()->GetPlanes();
  pCollCS->save(ar);  // since version 22.
}

static UINT __stdcall read_data_thread_func(LPVOID pData)
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  pObj->read_data();
  
  if(!pObj->get_terminate_flag()) // if the user has not terminated the dialog manually, do it after calculations have ended.
  {
    CDialog* pDlg = (CDialog*)pData;
    pDlg->PostMessage(WM_CLOSE);
  }

  return 0;
}

void CTracker::load(CArchive& ar)
{
  clear();
  clear_tracks();
  set_default();

// Parameters of the output are now properties of m_OutputEngine. 
  double fOutputTimeStep;
  bool bEnableFileOutput, bOnlyPassedQ00, bByInitRadii;
  UINT nEnsByRadiusCount;

  UINT nVersion;
  ar >> nVersion;

// General parameters:
  CString cFileName;
  ar >> cFileName;
  set_filename((const char*)cFileName);

  if(nVersion >= 1)
  {
    CString cFieldFileName;
    ar >> cFieldFileName;
//    set_field_filename((const char*)cFieldFileName);  // this function has been removed, but the variable is reserved.
  }

  ar >> m_bUseMultiThread;
  ar >> m_nType;
  double fTimeStep;
  ar >> fTimeStep;
  set_time_step(fTimeStep);
  ar >> m_fMaxIntegrTime;

  if(nVersion < 14)
    ar >> fOutputTimeStep;

  ar >> m_fInitD;
  ar >> m_nSymPlanes;

  if(nVersion >= 19)
    m_Transform.load(ar);

// Electrostatics:
  ar >> m_bEnableField;
  ar >> m_fAmplDC;
  ar >> m_fCharge;

// Ion parameters:
  double fIonMass;
  ar >> fIonMass;
  set_ion_mass(fIonMass); // to make sure that the run-time variable m_fChargeMassRatio has a correct value.
  ar >> m_fIonMobility;
// RF data:
  ar >> m_bEnableRF;
  ar >> m_fAmplRF;
  ar >> m_fOmega;

  if(nVersion >= 2)
  {
    ar >> m_fAmplRF_Q00;
    ar >> m_fOmega_Q00;
    ar >> m_fX_Q00;
  }

// Particle source:
  m_Src.load(ar);

// Evaporation:
  ar >> m_nEvaporModelType;
  create_evapor_model();
  if(m_pEvaporModel != NULL)
    m_pEvaporModel->load(ar);

  if((nVersion >= 2) && (nVersion < 14))
    ar >> bEnableFileOutput;

  if(nVersion >= 3)
  {
    CExportOpenFOAM* pExpObj = CParticleTrackingApp::Get()->GetExporter();
    pExpObj->set_default();
    pExpObj->load(ar);
  }

  if(nVersion >= 4)
  {
    if(nVersion < 14)
      ar >> bOnlyPassedQ00;

    ar >> m_fActEnergy;
  }

  if(nVersion >= 5)
  {
    ar >> m_bEnableCoulomb;
    ar >> m_fInitBunchRadius;
    ar >> m_fFullCurrent;
  }

  if(nVersion >= 11)
  {
    ar >> m_bAxialSymm;
    ar >> m_nIterCount;

    if(nVersion >= 13)
    {
      if(nVersion < 17)
      {
        UINT nPlanesCount;
        ar >> nPlanesCount;
      }
    }
    else
    {
      double fChargeTimeStep;
      ar >> fChargeTimeStep;
    }

    ar >> m_fBHDist;
    if(nVersion >= 12)
    {
      ar >> m_fBHCritR;
      ar >> m_nMaxRecDepth;
    }
  }

  if(nVersion >= 6)
    m_Importer.load(ar);

  if((nVersion >= 7) && (nVersion < 14))
  {
    ar >> bByInitRadii;
    ar >> nEnsByRadiusCount;
  }

  if(nVersion >= 8)
    ar >> m_bVelDependent;

  if(nVersion >= 9)
  {
    ar >> m_fAmplRF_Q0;
    ar >> m_fOmega_Q0;
    ar >> m_fX_Q0;
  }

// Calculators:
  if(nVersion >= 10)
  {
    CCalcCollection* pCalcColl = CParticleTrackingApp::Get()->GetCalcs();
    pCalcColl->load(ar);
  }

  if(nVersion < 14)
  {
    m_OutputEngine.set_enable_file_output(bEnableFileOutput);
    m_OutputEngine.set_enable_ens_by_radius(bByInitRadii);
    m_OutputEngine.set_restrict_output(bOnlyPassedQ00);
    m_OutputEngine.set_ens_by_radius_count(nEnsByRadiusCount);
    m_OutputEngine.set_output_time_step(fOutputTimeStep);
  }
  else
  {
    m_OutputEngine.load(ar);
  }

  if(nVersion >= 15)
  {
    ar >> m_bSaveTracks;
    if(m_bSaveTracks)
    {
      load_tracks(ar);

      std::string sFileName = COutputEngine::get_file_name(get_filename());   // take the file name from the data file.
      std::string sNewPath = COutputEngine::get_full_path((const char*)ar.m_strFileName); // take the path from the archive.
      m_sDataFile = sNewPath + sFileName;

      CExecutionDialog dlg(&read_data_thread_func, (CObject*)this);
      INT_PTR nRes = dlg.DoModal();

      if(m_bReady)
      {
        CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
        pDrawObj->set_hidden_reg_names();
        pDrawObj->draw();
      }

      m_Src.set_data_from_tracks();
    }
  }

  if(nVersion >= 16)
    ar >> m_nIntegrType;

  if(nVersion >= 18)
  {
    ar >> m_bOldIntegrator;
    ar >> m_bUseRadialCoulomb;
    ar >> m_fRadialCoulombX;
  }

  if(nVersion >= 17)
    m_SpaceChargeDist.load(ar);

  if(nVersion >= 20)
    m_vFieldPtbColl.load(ar); // collection of the external DC field perturbations.

  if(nVersion >= 21)
  {
    CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
    pFields->load(ar);
  }

  if(nVersion >= 22)
  {
    CCrossSectColl* pCollCS = CParticleTrackingApp::Get()->GetPlanes();
    pCollCS->load(ar);  // since version 22.
  }
  
// Derived variables:
  m_fInitMass = get_particle_mass(m_fInitD);

  update_interface();
}

void CTracker::save_track_const(CArchive& ar, const CTrack& track)
{
  const UINT nVersion = 1;  // Two different types of track items - for ions and for droplets.
  ar << nVersion;

  ar << track.get_type();
  ar << track.get_index();
  ar << track.get_phase();
  ar << track.get_current();
}

void CTracker::load_track_const(CArchive& ar, CTrack& track)
{
  UINT nVersion;
  ar >> nVersion;

  int nType;
  UINT nIndex;
  double fPhase, fCharge, fCurr, fGasMolMass, fMass;
  if(nVersion < 1)
  {
    ar >> fPhase;
    ar >> fCharge;
    ar >> fGasMolMass;
    ar >> fCurr;
    ar >> nIndex;

    nType = m_nType;
    if(nType == CTrack::ptIon)
      ar >> fMass;
  }
  else
  {
    ar >> nType;
    ar >> nIndex;
    ar >> fPhase;
    ar >> fCurr; 
  }

  track.set_type(nType);
  track.set_index(nIndex);
  track.set_current(fCurr);
  track.set_phase(fPhase);
}

void CTracker::save_tracks(CArchive& ar)
{
  const UINT nVersion = 2;  // Two different types of track items - for ions and for droplets.
  ar << nVersion;

  size_t nTrackCount = m_Tracks.size();
  ar << nTrackCount;
  if(nTrackCount == 0)
    return;

  for(size_t i = 0; i < nTrackCount; i++)
  {
    const CTrack& track = m_Tracks.at(i);
    size_t nItemCount = track.size();
    ar << nItemCount;
    if(nItemCount == 0)
      continue;

    save_track_const(ar, track);
    for(size_t j = 0; j < nItemCount; j++)
    {
      CBaseTrackItem* pItem = track.at(j);
      pItem->save(ar);
    }
  }
}

void CTracker::load_tracks(CArchive& ar)
{
  clear_tracks();

  UINT nVersion;
  ar >> nVersion;

  size_t nTrackCount;
  ar >> nTrackCount;
  if(nTrackCount == 0)
    return;

  Vector3D vDummy;
  for(size_t i = 0; i < nTrackCount; i++)
  {
    set_status("Loading tracks", 100 * i / nTrackCount);

    CTrack track;
    size_t nItemCount;
    ar >> nItemCount;
    if(nItemCount == 0)
      continue;

    load_track_const(ar, track);

    int nType = track.get_type();
    for(size_t j = 0; j < nItemCount; j++)
    {
      CBaseTrackItem* pItem = CBaseTrackItem::create(nType);
      if(nVersion < 2)
      {
        ar >> pItem->pos.x;
        ar >> pItem->pos.y;
        ar >> pItem->pos.z;
        ar >> pItem->vel.x;
        ar >> pItem->vel.y;
        ar >> pItem->vel.z;

        if(nVersion < 1)
        {
          ar >> vDummy.x;
          ar >> vDummy.y;
          ar >> vDummy.z;
        }

        double fDiffCoeff, fKappa, fTemp, fEqTemp;
        ar >> fDiffCoeff;
        ar >> fKappa;
        ar >> fTemp;
        ar >> fEqTemp;
        ar >> pItem->time;

        if(nType == CTrack::ptIon)
        {
          CIonTrackItem* pIonItem = (CIonTrackItem*)pItem;
          pIonItem->temp = fTemp;
          pIonItem->tempinf = fEqTemp;
        }
      }
      else
      {
        pItem->load(ar);
      }

      track.push_back(pItem);
    }

    m_Tracks.push_back(track);
  }

  set_status("Ready", -1);
}

void CTracker::update_interface()
{
  CMainFrame* pMainWnd = (CMainFrame*)(CParticleTrackingApp::Get()->m_pMainWnd);
  if(pMainWnd == NULL)
    return;

  CPropertiesWnd* pPropWnd = pMainWnd->GetWndProperties();
  pPropWnd->set_update_all();
}

void CTracker::set_data()
{
  CMainFrame* pMainWnd = (CMainFrame*)(CParticleTrackingApp::Get()->m_pMainWnd);
  if(pMainWnd == NULL)
    return;

  CPropertiesWnd* pPropWnd = pMainWnd->GetWndProperties();
  pPropWnd->set_data_to_model();
}

void CTracker::invalidate_calculators()
{
  CCalcCollection* pCalcColl = CParticleTrackingApp::Get()->GetCalcs();
  size_t nCalcCount = pCalcColl->size();
  for(size_t i = 0; i < nCalcCount; i++)
    pCalcColl->at(i)->invalidate();
}

const char* CTracker::get_integr_name(int nType) const
{
  switch(nType)
  {
    case intExplEuler: return _T("Boost Explicit Euler");
    case intModMidpnt: return _T("Boost Modified Midpoint");
    case intRK2: return _T("Runge-Kutta 2");
    case intRK4: return _T("Boost Runge-Kutta 4");
  }

  return _T("None");
}

CRegionsCollection& CTracker::get_regions()
{
  CCrossSectColl* pColl = CParticleTrackingApp::Get()->GetPlanes();
  size_t nPlanesCount = pColl->size();
  if(nPlanesCount == 0)
    return m_vRegions;

  m_vExtRegions.clear();
  size_t nRegCount = m_vRegions.size();
  size_t nExtRegCount = nRegCount + nPlanesCount;
  m_vExtRegions.reserve(nExtRegCount);

  for(size_t i = 0; i < nRegCount; i++)
    m_vExtRegions.push_back(m_vRegions.at(i));

  for(size_t j = 0; j < nPlanesCount; j++)
  {
    CDomainCrossSection* pPlane = pColl->at(j);
    CRegion* pReg = pPlane->get_region();
    if(pReg != NULL)
      m_vExtRegions.push_back(pReg);
  }

  return m_vExtRegions;
}

};  // namespace EvaporatingParticle