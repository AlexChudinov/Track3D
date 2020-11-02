#include "stdafx.h"

#include "BatchSim.h"
#include "Tracker.hpp"
#include "ParticleTracking.h"
#include "BarnesHut.h"
#include "Calculator.h"
#include "OutputEngine.h"

#include <sstream>
#include <algorithm>
#include <direct.h>

namespace EvaporatingParticle
{

static double fCoeffCS = Const_Angstrem_CGS * Const_Angstrem_CGS;
//-------------------------------------------------------------------------------------------------
// CIonSimParams - an auxiliary class for batch simulations. 
//-------------------------------------------------------------------------------------------------
CIonSimParams::CIonSimParams()
{
  fMass = 508 * Const_AMU_CGS;
  fCharge = Const_Charge_CGS;
  fMob = 1.05 / SI_to_CGS_Voltage;
  fCrossSect = 160 * fCoeffCS;
  fFullCurr = 100 * Const_nA_to_CGSE;
  sLegend = CString(" ");
}

//-------------------------------------------------------------------------------------------------
// CBatchSim - an auxiliary class for batch simulations. 
//-------------------------------------------------------------------------------------------------
CBatchSim::CBatchSim()
  : m_pObj(NULL), m_pPrevClmb(NULL), m_pPrevPhi(NULL), m_pCalc(NULL)
{
  set_default();
}

CBatchSim::~CBatchSim()
{
  clear_prev_clmb();
}

void CBatchSim::set_default()
{
  m_bEnable = false;
  m_bUseSlidingAverage = false;
  m_nStepsCount = 0;

  m_fIncrPerIter = 1.0 * Const_nA_to_CGSE;
  m_sFileName = CString(_T(""));
  set_aver_width(10);
}

void CBatchSim::set_curr_incr_iter(double fCurrIncr)
{
  m_fIncrPerIter = fCurrIncr;
  m_pObj = CParticleTrackingApp::Get()->GetTracker();
  m_pObj->set_iter_count(m_pObj->get_full_iter_count(m_fIncrPerIter));
}

bool CBatchSim::batch_calc_init()
{
  m_pObj = CParticleTrackingApp::Get()->GetTracker();

  backup_sim_params();

  if(!read_sim_params())
    return false;

// DEBUG
  debug_output();
// END DEBUG

  if(!get_calc())
    return false;

  if(m_bUseSlidingAverage)
    init_prev_clmb();

  release_calc();
  return true;
}

void CBatchSim::batch_calc_relax()
{
  if(m_pObj == NULL)
    return;

  set_sim_params_to_tracker(m_BackUpParams);  // restore the original ion and field parameters.
  clear_prev_clmb();
}

bool CBatchSim::read_sim_params()
{
  CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (const char*)m_sFileName, (const char*)("r"));
  if(nErr != 0 || pStream == 0)
    return false;

  int nRes = 0, nFieldsCount;
  nRes = fscanf_s(pStream, "%d, %d", &m_nStepsCount, &nFieldsCount);
  if(m_nStepsCount < 1 || nFieldsCount != pFields->size())
  {
    fclose(pStream);
    return false;
  }

  m_vSimParams.clear();

  char sHeader[16];
  for(UINT j = 0; j < 10; j++)
    nRes = fscanf(pStream, "%s", sHeader);

  double fMass, fCharge, fMob, fCS, fCurr, fScale, fFreq;
  for(UINT i = 0; i < m_nStepsCount; i++)
  {
    nRes = fscanf_s(pStream, "%lf, %lf, %lf, %lf, %lf,", &fMass, &fCharge, &fMob, &fCS, &fCurr);
    CIonSimParams param(fMass * Const_AMU_CGS, fCharge * Const_Charge_CGS, fMob / SI_to_CGS_Voltage, fCS * fCoeffCS, fCurr * Const_nA_to_CGSE);
    param.vFieldPar.resize(nFieldsCount);
    for(UINT k = 0; k < nFieldsCount; k++)
    {
//      nRes = (k == nFieldsCount - 1) ? fscanf_s(pStream, "%lf, %lf", &fScale, &fFreq) : fscanf_s(pStream, "%lf, %lf,", &fScale, &fFreq);
      nRes = fscanf_s(pStream, "%lf, %lf,", &fScale, &fFreq);
      param.vFieldPar[k].fAmpl = fScale;
      param.vFieldPar[k].fFreq = 1000 * fFreq;  // in the input file frequency is in kHz.
    }

    nRes = fscanf(pStream, "%s", sHeader);
    param.sLegend = CString(sHeader);

    m_vSimParams.push_back(param);
    if(nRes == EOF)
      break;
  }

  fclose(pStream);
  return true;
}

void CBatchSim::backup_sim_params()
{
  m_BackUpParams.fCharge = m_pObj->get_particle_charge();
  m_BackUpParams.fMass = m_pObj->get_ion_mass();
  m_BackUpParams.fMob = m_pObj->get_ion_mobility();

  if(!m_pObj->get_user_def_cs())
    m_pObj->calc_cross_section();

  m_BackUpParams.fCrossSect = m_pObj->get_ion_cross_section();
  m_BackUpParams.fFullCurr = m_pObj->get_full_current();

  CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
  size_t nFieldsCount = pFields->size();
  m_BackUpParams.vFieldPar.resize(nFieldsCount);
  for(size_t i = 0; i < nFieldsCount; i++)
  {
    CElectricFieldData* pData = pFields->at(i);
    m_BackUpParams.vFieldPar[i].fAmpl = pData->get_scale();
    m_BackUpParams.vFieldPar[i].fFreq = pData->get_freq();
  }
}

bool CBatchSim::set_sim_params_to_tracker(const CIonSimParams& param)
{
  m_pObj->set_particle_charge(param.fCharge);
  m_pObj->set_ion_mass(param.fMass);
  m_pObj->set_ion_mobility(param.fMob);
  if(!m_pObj->get_user_def_cs() || param.fCrossSect == 0)
    m_pObj->calc_cross_section();
  else
    m_pObj->set_ion_cross_section(param.fCrossSect);

  m_pObj->set_full_current(param.fFullCurr);

  CFieldDataColl* pFields = CParticleTrackingApp::Get()->GetFields();
  size_t nFieldsCount = pFields->size();
  if(nFieldsCount != param.vFieldPar.size())
  {
    AfxMessageBox("Count of fields in the scene and Ampl./Freq. pairs count in the input file do not match.");
    return false;
  }

  for(size_t i = 0; i < nFieldsCount; i++)
  {
    CElectricFieldData* pData = pFields->at(i);
    pData->set_scale(param.vFieldPar[i].fAmpl);
    pData->set_freq(param.vFieldPar[i].fFreq);
  }

  return true;
}

bool CBatchSim::prepare_step(UINT nStep)
{
  if(nStep >= m_vSimParams.size())
    return false;

  m_sLegend = m_vSimParams.at(nStep).sLegend;
  return set_sim_params_to_tracker(m_vSimParams.at(nStep));
}

//-------------------------------------------------------------------------------------------------
// Calculations and intermediate output section.
//-------------------------------------------------------------------------------------------------
bool CBatchSim::intermediate_output()
{
  if(m_pObj->get_terminate_flag())
    return true;

  CString sTransFile, sFragmFile, sIonTempFile, sSpaceChargeFile;
  if(!get_filenames(sTransFile, sFragmFile, sIonTempFile, sSpaceChargeFile))
    return false;

  if(!get_calc())
    return false;

  m_pCalc->set_clc_var_type(CTrackCalculator::clcCurrent);
  m_pCalc->set_filename((const char*)sTransFile);
  m_pCalc->do_sequence_calc();

  m_pCalc->set_clc_var_type(CTrackCalculator::clcFragment);
  m_pCalc->set_filename((const char*)sFragmFile);
  m_pCalc->do_sequence_calc();

  m_pCalc->set_clc_var_type(CTrackCalculator::clcIonTemp);
  m_pCalc->set_filename((const char*)sIonTempFile);
  m_pCalc->do_sequence_calc();

  release_calc();

  m_pObj->save_coulomb_field((const char*)sSpaceChargeFile); 
  return true;
}

bool CBatchSim::intermediate_save()
{
  CFile file;
  if(!file.Open((const char*)m_sTskFile, CFile::modeCreate | CFile::modeWrite))
  {
    AfxMessageBox("Can not open file for intermediate project output.");
    return false;
  }

  CArchive archive(&file, CArchive::store);

  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->save(archive);
  m_pObj->save(archive);
  archive.Flush();
  file.Close();

  return true;
}

bool CBatchSim::get_filenames(CString& sTransFile, CString& sFragmFile, CString& sIonTempFile, CString& sClmbFile)
{
  CString sPath(COutputEngine::get_full_path(m_pObj->get_filename()).c_str());

// Directory name by M/Z:
  double fMassAMU = m_pObj->get_ion_mass() / Const_AMU_CGS;
  double fCharge = m_pObj->get_particle_charge() / Const_Charge_CGS;
  int nMovrZ = int(0.5 + fMassAMU / fCharge);

  char buff[8];
  CString sMassOvrCharge(itoa(nMovrZ, buff, 10));
  CString sSlash("\\");
  CString sDirName = sPath + sMassOvrCharge + sSlash;
  if(_mkdir(sDirName) == ENOENT)
    return false;

  int nCurr = int(0.5 + m_pObj->get_full_current() / Const_nA_to_CGSE);
  CString sCurr(itoa(nCurr, buff, 10));

  CString sExt(".csv");
  sTransFile = sDirName + sMassOvrCharge + CString("_Trans_I=") + sCurr + CString("nA_") + m_sLegend + sExt;
  sFragmFile = sDirName + sMassOvrCharge + CString("_Fragm_I=") + sCurr + CString("nA_") + m_sLegend + sExt;
  sIonTempFile = sDirName + sMassOvrCharge + CString("_Ion_Temp_I=") + sCurr + CString("nA_") + m_sLegend + sExt;
  sClmbFile = sDirName + CString("Space_Charge_I=") + sCurr + CString("nA_") + m_sLegend + sExt;

  m_sTskFile = sPath + sMassOvrCharge + CString("_I=") + sCurr + CString("nA_") + m_sLegend + CString(".tsk");

  return true;
}

bool CBatchSim::get_calc()
{
  CCalcCollection* pCalcsColl = CParticleTrackingApp::Get()->GetCalcs();
  size_t nCalcCount = pCalcsColl->size();
  for(size_t i = 0; i < nCalcCount; i++)
  {
    CCalculator* pCalc = pCalcsColl->at(i);
    if(pCalc->type() == CCalculator::ctTrackCalc)
    {
      HWND hJobName, hProgress, hDlgWnd;
      m_pObj->get_handlers(hJobName, hProgress, hDlgWnd);
      pCalc->set_handlers(hJobName, hProgress, hDlgWnd);
      m_pCalc = (CTrackCalculator*)pCalc;
      return true;
    }
  }

  AfxMessageBox("Track calculator was not found. Ctreate the calculator and set proper parameters for sequential calculations.");
  return false;
}

void CBatchSim::release_calc()
{
  if(m_pCalc != NULL)
  {
    m_pCalc->set_handlers(NULL, NULL, NULL);
    m_pCalc = NULL;
  }
}

void CBatchSim::debug_output()
{
  CString sPath(COutputEngine::get_full_path(m_pObj->get_filename()).c_str());
  CString sFileName = sPath + CString("Batch_Data_Read.txt");

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (const char*)sFileName, (const char*)("w"));
  if(nErr != 0 || pStream == 0)
  {
    AfxMessageBox("Cannot open 'Batch_Data_Read.txt' file for writing.");
    return;
  }

  int nStepsCount = m_vSimParams.size();
  int nFieldCount = nStepsCount > 0 ? m_vSimParams.at(0).vFieldPar.size() : 0;
  fprintf(pStream, "%d,  %d\n", nStepsCount, nFieldCount);
  for(UINT i = 0; i < nStepsCount; i++)
  {
    const CIonSimParams& params = m_vSimParams.at(i);

    fprintf(pStream, "%10.2lf,%7.1lf,%8.2lf,%8.1lf,%8.1lf,",
      params.fMass / Const_AMU_CGS, params.fCharge / Const_Charge_CGS, params.fMob * SI_to_CGS_Voltage, params.fCrossSect / fCoeffCS, params.fFullCurr / Const_nA_to_CGSE);

    for(UINT j = 0; j < nFieldCount; j++)
      fprintf(pStream, "%10.2lf,%10.2lf,", params.vFieldPar[j].fAmpl, 0.001 * params.vFieldPar[j].fFreq);

    fputs((const char*)(params.sLegend), pStream);
    fputs("\n", pStream);
  }

  fclose(pStream);
}

//-------------------------------------------------------------------------------------------------
// Auxiliary arrays allocation / deallocation.
//-------------------------------------------------------------------------------------------------
void CBatchSim::init_prev_clmb()
{
  if(m_nAverWidth == 0)
    return;

  clear_prev_clmb();

  size_t nNodesCount = m_pObj->get_nodes().size();
  m_pPrevClmb = new Vector3F*[m_nAverWidth];
  m_pPrevPhi = new float*[m_nAverWidth];
  for(UINT i = 0; i < m_nAverWidth; i++)
  {
    m_pPrevClmb[i] = new Vector3F[nNodesCount];
    m_pPrevPhi[i] = new float[nNodesCount];
  }
}

void CBatchSim::clear_prev_clmb()
{
  if(m_nAverWidth == 0 || m_pPrevClmb == NULL)
    return;

  size_t nNodesCount = m_pObj->get_nodes().size();
  for(UINT i = 0; i < m_nAverWidth; i++)
  {
    delete[] m_pPrevClmb[i];
    delete[] m_pPrevPhi[i];
  }

  delete[] m_pPrevClmb;
  delete[] m_pPrevPhi;
  m_pPrevClmb = NULL;
  m_pPrevPhi = NULL;
}

//-------------------------------------------------------------------------------------------------
// Streaming:
//-------------------------------------------------------------------------------------------------
void CBatchSim::save(CArchive& ar)
{
  const UINT nVersion = 1;  // arbitrary simulations list instead of a fixed currents set.
  ar << nVersion;

  ar << m_bEnable;
  ar << m_bUseSlidingAverage;

  ar << m_fIncrPerIter;
  ar << m_nAverWidth;

// Since version 1:
  ar << m_sFileName;
}

void CBatchSim::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_bEnable;
  ar >> m_bUseSlidingAverage;

  ar >> m_fIncrPerIter;
  ar >> m_nAverWidth;

  if(nVersion == 0)
  {
    size_t nStepsCount;
    ar >> nStepsCount;
    double fCurr;
    for(size_t i = 0; i < nStepsCount; i++)
      ar >> fCurr;
  }

  if(nVersion >= 1)
    ar >> m_sFileName;
}

void CBatchSim::accum_clmb_field_in_node(CNode3D* pNode, CBarnesHut* pBHObj, CElectricFieldData* pData, UINT nIter)
{
  bool bFieldExists = pData != NULL;
  bool bMirrEnabled = bFieldExists && pData->get_enable_field();
  size_t nInd = pNode->nInd;

  Vector3D vMirrField = bMirrEnabled ? pData->get_field(nInd) : vNull;
  Vector3D vClmbField = pBHObj->coulomb_force(pNode->pos) + vMirrField;

// Visualization of the Coulomb potential.
  float fMirrPhi = bMirrEnabled ? pData->get_phi(nInd) : 0.0f;
  float fClmbPhi = bFieldExists ? CGS_to_SI_Voltage * (pBHObj->coulomb_phi(pNode->pos) + fMirrPhi) : 0.0f;
  float fClmbPhiStored = fClmbPhi;

  pNode->clmb = vClmbField;
  if(nIter <= m_nAverWidth) // In the course of the first m_nAverWidth iterations m_pPrevClmb[i] are being defined.
  {
    UINT nAverWidth = nIter - 1;
    double fNormCoeff = 1. / (double)nIter;
    for(UINT i = 0; i < nAverWidth; i++)
      pNode->clmb += m_pPrevClmb[i][nInd];

    pNode->clmb *= fNormCoeff;

    for(int j = nAverWidth; j >= 0; j--)
      m_pPrevClmb[j][nInd] = j == 0 ? vClmbField : m_pPrevClmb[j - 1][nInd];

    if(bFieldExists)
    {
// Visualization of the Coulomb potential.
      for(UINT i = 0; i < nAverWidth; i++)
        fClmbPhi += m_pPrevPhi[i][nInd];

      fClmbPhi *= fNormCoeff;
      pData->set_clmb_phi(nInd, fClmbPhi);

      for(int j = nAverWidth; j >= 0; j--)
        m_pPrevPhi[j][nInd] = j == 0 ? fClmbPhiStored : m_pPrevPhi[j - 1][nInd];
    }
  }
  else  // nIter > m_nAverWidth and all m_pPrevClmb[i] are already defined.
  {
    for(UINT i = 0; i < m_nAverWidth; i++)
      pNode->clmb += m_pPrevClmb[i][nInd];

    pNode->clmb *= m_fNormCoeff;

// Re-assigning the previous values of field for correct averaging at the next iteration.
    for(int j = m_nAverWidth - 1; j >= 0; j--)
      m_pPrevClmb[j][nInd] = j == 0 ? vClmbField : m_pPrevClmb[j - 1][nInd];

    if(bFieldExists)
    {
// Visualization of the Coulomb potential.
      for(UINT i = 0; i < m_nAverWidth; i++)
        fClmbPhi += m_pPrevPhi[i][nInd];

      fClmbPhi *= m_fNormCoeff;
      pData->set_clmb_phi(nInd, fClmbPhi);

      for(int j = m_nAverWidth - 1; j >= 0; j--)
        m_pPrevPhi[j][nInd] = j == 0 ? fClmbPhiStored : m_pPrevPhi[j - 1][nInd];
    }
  }
}


};  // namespace EvaporatingParticle
