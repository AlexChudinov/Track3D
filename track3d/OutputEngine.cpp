
#include "stdafx.h"

#include <direct.h>

#include "OutputEngine.h"
#include "ParticleTracking.h"
#include "EvaporationModel.h"
#include "BeamCrossSection.h"
#include "ColorContour.h"

#include "vector2d.hpp"
#include "mathematics.h"
#include "constant.hpp"
#include "float.h"
#include "math.h"


namespace EvaporatingParticle
{

COutputEngine::COutputEngine()
{
  set_default();
}

COutputEngine::~COutputEngine()
{
  m_AverEngine.clear();
}

void COutputEngine::set_default()
{
  m_bEnableFileOutput = false;

  m_fOutputTimeStep = 1.0e-6; // sec.
  
  m_bOnlyPassedQ00 = false;
  m_bByInitRadii = true;  // output ensemble will be formed by the initial ion radii instead if ensemble index.

  m_nEnsByRadiusCount = 8;

  m_fCrssSctX = 1.0;  // cm, a default cross-section position.
}

std::string COutputEngine::get_full_path(const char* pFullPath)
{
  std::string cPath(pFullPath);
  while(cPath[cPath.size() - 1] != '\\')
    cPath.erase(cPath.end() - 1);

  return cPath;
}

std::string COutputEngine::get_file_name(const char* pFullPath)
{
  std::string cPath(pFullPath);
  size_t nSlashPos = cPath.size() - 1;
  while(cPath[nSlashPos] != '\\')
    nSlashPos--;

  cPath.erase(cPath.begin(), cPath.begin() + nSlashPos + 1);
  return cPath;
}

std::string COutputEngine::get_base_name(const char* pFullPath)
{
  std::string cPath(pFullPath);
  while(cPath[cPath.size() - 1] != '.')
    cPath.erase(cPath.end() - 1);

  return cPath;
}

bool COutputEngine::output_droplet_track(UINT nTrackIndex)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CTrackVector& vTracks = pObj->get_tracks();
  if(nTrackIndex >= vTracks.size())
    return false;

  std::string cPath = get_full_path(pObj->get_filename());

  char buff[8];
  std::string cIndex("");
  if(_itoa_s(nTrackIndex, buff, 8, 10) == 0)
    cIndex = std::string(buff);

  std::string cTrack("track_drop_");
  std::string cExt(".dat");
  std::string cFileName = cPath + cTrack + cIndex + cExt;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return false;

  CEvaporationModel* pEvaporModel = pObj->get_evapor_model();
  double fInitD = pObj->get_init_diameter();
  if(fInitD < Const_Almost_Zero)
    return false;

// Header:
  fprintf(pStream, "%s\n", "time(ms), x(cm), y(cm), d/d0, Twater(K), Tenv(K), Tsat(K), Q/Qr");

  CNode3D node;
  CElem3D* pElem = NULL;
  double fTime_ms, fDropD, fD_part, fDropTemp, fEnvTemp, fSatTemp, fRayleighLim, fRatioQ;

  CTrackItem item;
  const CTrack& track = vTracks.at(nTrackIndex);
  size_t nCount = track.size();
  for(size_t i = 0; i < nCount; i++)
  {
    track.get_track_item(i, item);

    fTime_ms = item.time * 1000;
    fDropD = pObj->get_particle_diameter(item.mass);
    fD_part = fDropD / fInitD;
    fDropTemp = item.temp;

    pObj->interpolate(item.pos, 0, 0, node, pElem);
    fEnvTemp = node.temp;
    pEvaporModel->get_saturation_temp(node.press, pEvaporModel->m_H2O_mf, fSatTemp);

    fRayleighLim = pObj->get_max_charge(fDropTemp, fDropD);
    fRatioQ = fRayleighLim > Const_Almost_Zero ? pObj->get_particle_charge() / fRayleighLim : FLT_MAX;

    fprintf(pStream, "%f, %f, %f, %f, %f, %f, %f, %f\n",
      fTime_ms, item.pos.x, item.pos.y, fD_part, fDropTemp, fEnvTemp, fSatTemp, fRatioQ);
  }

  fclose(pStream);
  return true;
}

UINT COutputEngine::get_max_ensemble_index() const
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CTrackVector& vTracks = pObj->get_tracks();

  UINT nMaxIndex = 0;
  size_t nTrackCount = vTracks.size();
  for(size_t i = 0; i < nTrackCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    if(nMaxIndex < track.get_index())
      nMaxIndex = track.get_index();
  }

  return nMaxIndex;
}

void COutputEngine::get_output_range(double& fArgMin, double& fArgMax)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CTrackVector& vTracks = pObj->get_tracks();

  fArgMin = FLT_MAX;
  fArgMax = -FLT_MAX;

  size_t nTrackCount = vTracks.size();
  for(size_t i = 0; i < nTrackCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    size_t nSize = track.size();
    if(nSize < 2)
      continue;

    if(fArgMin > track.at(0)->pos.x)
      fArgMin = track.at(0)->pos.x;
    if(fArgMax < track.at(nSize - 1)->pos.x)
      fArgMax = track.at(nSize - 1)->pos.x;
  }
}

//---------------------------------------------------------------------------------------
// Output ion ensembles. There will be as many files as many ensembles have started.
//---------------------------------------------------------------------------------------
void COutputEngine::output_ion_tracks()
{
  set_job_name("Writing data to files...");
  set_progress(0);

  m_AverEngine.clear();
  m_AverEngine.set_type(atAverage);

  double fMinX, fMaxX;
  get_output_range(fMinX, fMaxX);
  m_AverEngine.set_range(fMinX, fMaxX, 100);

  if(m_bByInitRadii)
    output_ensemble_by_initial_radius();
  else
    output_ensemble_by_ens_index();

  if(m_bTerminate)
    return;

  output_average_tracks();  // support of averaging over all the tracks to be output.
}

void COutputEngine::output_ensemble_by_ens_index()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CTrackVector& vTracks = pObj->get_tracks();
  double fX_Q00 = pObj->get_Q00_trans();

// Progress bar support:
  int nPercentStart = 0, nPercentEnd;

  CTrackItem item;
  UINT nEnsIndex = 0, nMaxIndex = get_max_ensemble_index();
  while(nEnsIndex <= nMaxIndex)
  {
    CTrackVector vEnsTracks;
    size_t nTrackCount = vTracks.size();
    nPercentEnd = nMaxIndex > 0 ? int(0.5 + 100.0 * nEnsIndex / nMaxIndex) : 100;
    for(size_t i = 0; i < nTrackCount; i++)
    {
      const CTrack& track = vTracks.at(i);
      size_t nSize = track.size();
      if(nSize < 2)
        continue;

      if(m_bOnlyPassedQ00)  // if m_bOnlyPassedQ00 == true only those tracks will be output which have passed to the Q00 region.
      {
        track.get_track_item(nSize - 1, item);
        if(item.pos.x < fX_Q00)
          continue;
      }

      if(track.get_index() == nEnsIndex)
        vEnsTracks.push_back(track);  // collection of ions with one and the same ensemble index.
    }

    if(vEnsTracks.size() == 0)
    {
      nEnsIndex++;
      set_progress(int(0.5 + 100. * nEnsIndex / (1 + nMaxIndex)));
      continue;
    }

    std::string cPath = get_full_path(pObj->get_filename());

    char buff[8];
    std::string cIndex("");
    if(_itoa_s(nEnsIndex, buff, 8, 10) == 0)
      cIndex = std::string(buff);

    std::string cTrack("track_ion_");
    std::string cExt(".dat");
    std::string cFileName = cPath + cTrack + cIndex + cExt;

    FILE* pStream;
    errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
    if(nErr != 0 || pStream == 0)
      continue;

    output_ion_ensemble(vEnsTracks, pStream, nPercentStart, nPercentEnd);

    fclose(pStream);
    nEnsIndex++;

    nPercentStart = nPercentEnd;
    if(m_bTerminate)
      return;
  }
}

//---------------------------------------------------------------------------------------
// This function also uses COutputEngine::output_ion_ensemble() but forms an ensemble in a
// different way: here an ensemble contains ions which starting radii satisfy ra < r0 < rb.
//---------------------------------------------------------------------------------------
void COutputEngine::output_ensemble_by_initial_radius()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CTrackVector& vTracks = pObj->get_tracks();
  bool bEnableCoulomb = pObj->get_enable_coulomb();
  double fX_Q00 = pObj->get_Q00_trans();

  size_t nTrackCount = vTracks.size();
  if(m_nEnsByRadiusCount == 0 || nTrackCount == 0)
    return;

// Progress bar support:
  int nPercentStart = 0, nPercentEnd; 

  CTrackItem item;
  CTrackVector vEnsTracks;
  double fRmax = pObj->get_src()->get_radius();
  double fR, fRa, fRb, fdR = fRmax / m_nEnsByRadiusCount;
  Vector2D vInitPos;
  for(UINT i = 0; i < m_nEnsByRadiusCount; i++)
  {
    fRa = i * fdR;
    fRb = fRa + fdR;
    nPercentEnd = int(0.5 + 100.0 * (i + 1) / m_nEnsByRadiusCount);

    vEnsTracks.clear();
    for(size_t j = 0; j < nTrackCount; j++)
    {
      const CTrack& track = vTracks.at(j);
      size_t nSize = track.size();
      if(nSize < 2)
        continue;

      if(m_bOnlyPassedQ00)  // if m_bOnlyPassedQ00 == true only those tracks will be output which have passed to the Q00 region.
      {
        track.get_track_item(nSize - 1, item);
        if(item.pos.x < fX_Q00)
          continue;
      }

      track.get_track_item(0, item);
      vInitPos = Vector2D(item.pos.y, item.pos.z);
      fR = vInitPos.length();
      if(fR < fRa || fR > fRb)
        continue;

      vEnsTracks.push_back(track);
    }

    if(vEnsTracks.size() == 0)
      continue;

    std::string cPath = get_full_path(pObj->get_filename());

    char buff[8];
    int nFullCurr = bEnableCoulomb ? int(0.5 + pObj->get_full_current() / Const_nA_to_CGSE) : 0;
    std::string cCurr(itoa(nFullCurr, buff, 10));

    std::string cEnd("nA\\");
    std::string cName("I=");
    std::string cDirName = cPath + cName + cCurr + cEnd;
    if(_mkdir(cDirName.c_str()) == ENOENT)
      continue;

    std::string cIndex("");
    if(_itoa_s(i, buff, 8, 10) == 0)
      cIndex = std::string(buff);

    std::string cTrack("track_ens_by_r0_");
    std::string cExt(".dat");
    std::string cFileName = cDirName + cTrack + cIndex + cExt;

    FILE* pStream;
    errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
    if(nErr != 0 || pStream == 0)
      continue;

    fprintf(pStream, "%f < r < %f\n", fRa, fRb);

    output_ion_ensemble(vEnsTracks, pStream, nPercentStart, nPercentEnd);

    fclose(pStream);
    if(m_bTerminate)
      return;

    nPercentStart = nPercentEnd;
  }
}

void COutputEngine::output_ion_ensemble(const CTrackVector& vEnsTracks, FILE* pStream, int nPercentStart, int nPercentEnd)
{
  CNode3D node;
  CElem3D* pElem = NULL;
  double fTime_ms, fIonMob, fGasTemp, fAbsRelVel, fNumberDens, fTi;
  Vector3D vDriftVel;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  double fCrossSection = pObj->get_ion_cross_section();
  double fActEnergy_eV = pObj->get_act_energy();
  double fIonMass = pObj->get_ion_mass();

  const double fIonMassAMU = fIonMass / Const_AMU_CGS;
  const double fMassRatio = Const_Molar_Mass_Air / (Const_Molar_Mass_Air + fIonMassAMU);   // dimensionless.
  const double fReducedMass = fIonMass * fMassRatio;
  const double fTempCoeff = 1. / Const_Boltzmann;
  const double fAirMass = Const_Molar_Mass_Air * Const_AMU_CGS;   // mass of the air molecule, g.
  const double fM_ovr_m = fIonMassAMU / Const_Molar_Mass_Air;
  const double fActEnergy = fActEnergy_eV / Const_Erg_to_EV;       // activation energy in CGS.

  const double fMeanVel300K = sqrt(8 * Const_Boltzmann * 300. / (Const_PI * fAirMass)); // <V> ~ 4.68e+4 cm/s.
  const double fRelVelSigma = fMeanVel300K * fCrossSection;

  double fCollCount = 0, fPrevTime = 0, fFragmConst = 0, fFragmPart, fCollFreqIncr;

  size_t nTimeCount = vEnsTracks.at(0).size();
  size_t nEnsSize = vEnsTracks.size();

// Header:
  fprintf(pStream, "%s\n", "time(ms),    x(mm),    y(mm),    z(mm),    Tion(K), Prob(E>1eV), Tenv(K),   Coll_Freq, Fragm_Percent");

  CTrackItem item;
  double fDiffVelSqr;
  for(size_t i = 1; i < nTimeCount; i++)  // loop over time steps.
  {
    if(i % 3 != 0)
      continue;

    double fKsi = double(i + 1) / nTimeCount;
    set_progress(int((1. - fKsi) * nPercentStart + fKsi * nPercentEnd));
    if(m_bTerminate)
      return;

    double fCurrTime = vEnsTracks.at(0).at(i)->time;
    double fTime_ms = fCurrTime * 1000;
    double fdT = fCurrTime - fPrevTime;

    double fAverTi = 0;
    double fCollFreq = 0;
    double fFragmProb = 0;
    double fAverGasT = 0;
    Vector3D vAverPos(0, 0, 0); // average position in millimeters.

    UINT nIonsCount = 0;  // count of existing ions in the ensemble at the i-th time step.
    bool bVelDependent = pObj->get_vel_depend_flag();
    for(size_t j = 0; j < nEnsSize; j++)  // averaging over the ensemble.
    {
      const CTrack& track = vEnsTracks.at(j);
      if(i >= track.size())
        continue;

      track.get_track_item(i, item);
      if(!pObj->interpolate(item.pos, 0, 0, node, pElem))
        continue;

      fGasTemp = node.temp;
      fDiffVelSqr = (item.vel - node.vel).sqlength();
      fTi = fGasTemp + 0.2 * fTempCoeff * fAirMass * fDiffVelSqr;  // 0.2 = 1/5, the coefficient must be m/5kT.

      fNumberDens = node.dens / fAirMass;
      fAbsRelVel = sqrt(fDiffVelSqr + 3 * Const_Boltzmann * node.temp / fAirMass);

      fAverTi += fTi;
      fCollFreqIncr = bVelDependent ? fNumberDens * fRelVelSigma : fNumberDens * fAbsRelVel * fCrossSection;
      fCollFreq += fCollFreqIncr;
      fFragmProb += CMath::energy_probability_integral(fActEnergy * fTempCoeff / fTi);
      vAverPos += item.pos;
      nIonsCount++;
    }

    if(nIonsCount == 0)
      continue;

    double fCoeff = 1. / (double)nIonsCount;
    fAverTi *= fCoeff;
    fCollFreq *= fCoeff;
    fFragmProb *= fCoeff;
    vAverPos *= fCoeff;

    fCollCount += fCollFreq * fdT;

// dn/dt = -K;
// K = n * fCollFreq * fFragmProb;
// dn/n = -fCollFreq * fFragmProb * dt;
// n = n0 * exp(-Int(fCollFreq * fFragmProb * dt));
// fFragmPart = 1 - n/n0.
//
    fFragmConst += fFragmProb * fCollFreq * fdT;
    fFragmPart = 1. - exp(-fFragmConst);

// Averaging over all tracks support:
    CAverValue vAverVal;
    vAverVal.fX = vAverPos.x;
    vAverVal.fY1 = fFragmPart;
    vAverVal.fY2 = fFragmProb;
    m_AverEngine.add_value(vAverVal);

    vAverPos *= 10; // to make mm from cm.

    fprintf(pStream, "%f, %f, %f, %f, %f, %e, %f, %e, %f\n", 
      fTime_ms, vAverPos.x, vAverPos.y, vAverPos.z, fAverTi, fFragmProb, fGasTemp, fCollFreq, 100 * fFragmPart);

    fPrevTime = fCurrTime;
  }  
}

void COutputEngine::output_average_tracks()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  std::string cPath = get_full_path(pObj->get_filename());

  std::string cTrack("track_ion_average");
  std::string cExt(".dat");
  std::string cFileName = cPath + cTrack + cExt;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  fprintf(pStream, "%s\n", "x(mm),  Aver_Fragm_Percent, Aver_Fragm_Prob");

  UINT nBinsCount = m_AverEngine.get_count();
  for(UINT i = 0; i < nBinsCount; i++)
  {
    const CAverBin& bin = m_AverEngine.get_aver_bin(i);
    if(bin.get_count() == 0)
      continue;

    const CAverValue& vAverVal = bin.get();

    fprintf(pStream, "%f, %f, %f\n", 10 * vAverVal.fX, 100 * vAverVal.fY1, vAverVal.fY2);
  }

  fclose(pStream);
}

void COutputEngine::output_ion_current()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  std::string cPath = get_full_path(pObj->get_filename());
  std::string cName("ion_current_");

  char buff[4];
  int nFullCurr = int(0.5 + pObj->get_full_current() / Const_nA_to_CGSE);
  std::string cCurr(itoa(nFullCurr, buff, 10));

  std::string cExt("nA.dat");
  std::string cFileName = cPath + cName + cCurr + cExt;

  if(pObj->get_axial_symm())
    output_axial_current(cFileName.c_str());
  else
    output_averaged_current(cFileName.c_str());
}

void COutputEngine::output_axial_current(const char* cFileName)
{
  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName, (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CTrackVector& vTracks = pObj->get_tracks();

  fprintf(pStream, "%s\n", "x(mm),  current(nA)");  // header.

  double fMinX, fMaxX;
  get_output_range(fMinX, fMaxX);

  size_t nTrackCount = vTracks.size();
  double fFullCurrent = pObj->get_full_current();
  double fElemCurr = fFullCurrent / nTrackCount;

  const double fStep = 0.01; // cm.
  double fX = fMinX, fCurr;
  while(fX <= fMaxX)
  {
    fCurr = fFullCurrent;
    for(size_t i = 0; i < nTrackCount; i++)
    {
      const CTrack& track = vTracks.at(i);
      size_t nTimeCount = track.size();
      if(nTimeCount == 0 || track.at(nTimeCount - 1)->pos.x < fX)
        fCurr -= fElemCurr;
    }

    fprintf(pStream, "%f, %f\n", 10 * fX, fCurr / Const_nA_to_CGSE);
    fX += fStep;
  }

  fclose(pStream);
}

void COutputEngine::output_averaged_current(const char* cFileName)
{
  UINT nBinsCount = m_AverEngine.get_count();
  if(nBinsCount == 0)
    return;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName, (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

  fprintf(pStream, "%s\n", "x(mm),  current(nA)");  // header.

  double fMinX, fMaxX;
  m_AverEngine.get_range(fMinX, fMaxX);

  double fX = 10 * fMinX;
  double fCurr = pObj->get_full_current() / Const_nA_to_CGSE;
  fprintf(pStream, "%f, %f\n", fX, fCurr);
  
  for(UINT i = 0; i < nBinsCount; i++)
  {
    const CAverBin& bin = m_AverEngine.get_aver_bin(i);
    if(bin.get_count() == 0)
      continue;

    const CAverValue& val = bin.get();
    fX = 10 * val.fX;   // position in mm.
    fCurr = val.fY1 / Const_nA_to_CGSE; // current in nA.
    fprintf(pStream, "%f, %f\n", fX, fCurr);
  }

  fclose(pStream);
}

void COutputEngine::prepare_current_output()
{
  m_AverEngine.clear();
  double fArgMin, fArgMax;
  get_output_range(fArgMin, fArgMax);
  m_AverEngine.set_range(fArgMin, fArgMax, 150);
  m_AverEngine.set_type(atAverage);
}

void COutputEngine::add_current()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CTrackVector& vTracks = pObj->get_tracks();
  size_t nTrackCount = vTracks.size();
  if(nTrackCount == 0)
    return;

  double fMinX, fMaxX;
  get_output_range(fMinX, fMaxX);

  double fFullCurrent = pObj->get_full_current();
  double fElemCurr = fFullCurrent / nTrackCount;

  const double fStep = 0.01; // cm.
  double fX = fMinX, fCurr;
  while(fX <= fMaxX)
  {
    fCurr = fFullCurrent;
    for(size_t i = 0; i < nTrackCount; i++)
    {
      const CTrack& track = vTracks.at(i);
      size_t nTimeCount = track.size();
      if(nTimeCount == 0 || track.at(nTimeCount - 1)->pos.x < fX)
        fCurr -= fElemCurr;
    }

    CAverValue val(fX, fCurr);
    m_AverEngine.add_value(val);

    fX += fStep;
  }
}

void COutputEngine::out_cross_section()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  std::string cPath = get_full_path(pObj->get_filename());
  std::string cName("cross-section_x_=_");

  double fHalf = m_fCrssSctX >= 0 ? 0.5 : -0.5;
  double fCrssSctX = 0.001 * int(fHalf + 10000 * m_fCrssSctX);  // cross-section X in mm with accuracy 2 digits after decimal point.
  char buff[_CVTBUFSIZE];
  if(_gcvt(fCrssSctX, 3, buff) == NULL) // convertation double into array of char.
    return;

  std::string cX(buff);
  std::string cExt("_mm.dat");
  std::string cFileName = cPath + cName + cX + cExt;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  fprintf(pStream, "%s\n\n", " z(mm),  y(mm),  Color_R0,  Color_Ti");   // header.

  CColoredCrossSection cs1;
  cs1.set_cross_sect_pos(m_fCrssSctX);
  cs1.set_var(CColoredCrossSection::varStartRadius);
  cs1.set_levels_count(20);
  cs1.draw();

  CColoredCrossSection cs2;
  cs2.set_cross_sect_pos(m_fCrssSctX);
  cs2.set_var(CColoredCrossSection::varIonTemp);
  cs2.set_levels_count(20);
  cs2.draw();

  double fVal1, fVal2;
  Vector3D vPoint;
  COLORREF nColor;
  RGB_Color clr1, clr2;
  size_t nCount = min(cs1.get_points_count(), cs2.get_points_count());
  for(size_t i = 0; i < nCount; i++)
  {
    cs1.get_colored_point(i, vPoint, fVal1, clr1);
    cs2.get_colored_point(i, vPoint, fVal2, clr2);
    fprintf(pStream, "%f, %f, %d, %d\n", 10 * vPoint.z, 10 * vPoint.y, RGB(clr1.red, clr1.green, clr1.blue), RGB(clr2.red, clr2.green, clr2.blue));
  }

  fclose(pStream);
}

// Streams:
void COutputEngine::save(CArchive& ar)
{
  const UINT nVersion = 1;
  ar << nVersion;

  ar << m_fOutputTimeStep;
  ar << m_bEnableFileOutput;
  ar << m_bOnlyPassedQ00;
  ar << m_bByInitRadii;
  ar << m_nEnsByRadiusCount;
  ar << m_fCrssSctX;
}

void COutputEngine::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_fOutputTimeStep;
  ar >> m_bEnableFileOutput;
  ar >> m_bOnlyPassedQ00;
  ar >> m_bByInitRadii;
  ar >> m_nEnsByRadiusCount;
  if(nVersion >= 1)
    ar >> m_fCrssSctX;
}

};  // namespace EvaporatingParticle