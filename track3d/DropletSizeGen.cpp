#include "stdafx.h"

#include "DropletSizeGen.h"
#include "ParticleTracking.h"
#include "OutputEngine.h"     // for debug_output() only.

namespace EvaporatingParticle
{

CDropletSizeGen::CDropletSizeGen()
{
  set_default();
}

void CDropletSizeGen::set_default()
{
  m_nRandSeed = 14101963;
  m_nDistrType = dstGauss;
// Gaussian distribution parameters:
  m_fMeanSize = 1.0e-4;   // cm, the default droplet diameter is 1 mcm.
  m_fStdDev = 2.0e-5;
// Uniform distribution parameters:
  m_fMinSize = 0.8e-4;
  m_fMaxSize = 1.2e-4;

  size_t nCount = 100;
  set_droplets_count(nCount);
}

void CDropletSizeGen::generate()
{
  m_RndGen.seed(m_nRandSeed);
  switch(m_nDistrType)
  {
    case dstUniform:
    {
      std::uniform_real_distribution<double> uniDistr(m_fMinSize, m_fMaxSize);
      size_t nCount = m_vSizeArr.size();
      for(size_t i = 0; i < nCount; i++)
        m_vSizeArr[i] = uniDistr(m_RndGen);

      break;
    }
    case dstGauss:
    {
      std::normal_distribution<double> normDistr(m_fMeanSize, m_fStdDev);
      size_t nCount = m_vSizeArr.size();
      for(size_t i = 0; i < nCount; i++)
        m_vSizeArr[i] = normDistr(m_RndGen);

      break;
    }
  }

  debug_output();
}

double CDropletSizeGen::get_droplet_size(size_t i) const
{
  if(m_nDistrType == dstConst)
    return m_fMeanSize;

  return i < m_vSizeArr.size()? m_vSizeArr.at(i) : 0;
}

const char* CDropletSizeGen::distr_name(int nType)
{
  switch(nType)
  {
    case dstConst: return _T("Constant Value");
    case dstUniform: return _T("Uniform Distribution");
    case dstGauss: return _T("Gaussian Distribution");
  }

  return _T("");
}

void CDropletSizeGen::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  ar << m_nRandSeed;

  size_t nCount = get_droplets_count();
  ar << nCount;

  ar << m_fMeanSize;
  ar << m_fStdDev;
  ar << m_fMinSize;
  ar << m_fMaxSize;

  ar << m_nDistrType;
}

void CDropletSizeGen::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_nRandSeed;

  size_t nCount;
  ar >> nCount;
  set_droplets_count(nCount);

  ar >> m_fMeanSize;
  ar >> m_fStdDev;
  ar >> m_fMinSize;
  ar >> m_fMaxSize;

  ar >> m_nDistrType;
}

void CDropletSizeGen::debug_output()
{
  if(m_nDistrType == dstConst)
    return;

  double fD, fMinD = FLT_MAX, fMaxD = -FLT_MAX;
  size_t nCount = m_vSizeArr.size();
  for(size_t i = 0; i < nCount; i++)
  {
    fD = m_vSizeArr.at(i);
    if(fD < fMinD)
      fMinD = fD;
    if(fD > fMaxD)
      fMaxD = fD;
  }

  UINT nBinCount = 20;
  double fStep = (fMaxD - fMinD) / nBinCount;
  if(fStep < Const_Almost_Zero)
    return;

  double* pHist = new double[nBinCount];
  for(UINT j = 0; j < nBinCount; j++)
    pHist[j] = 0.0;

  for(size_t i = 0; i < nCount; i++)
  {
    fD = m_vSizeArr.at(i);
    UINT j = UINT((fD - fMinD) / fStep);
    if(j == nBinCount)
      j = nBinCount - 1;

    pHist[j] += 1.0;
  }

  double fCount, fMax = 0;
  for(UINT j = 0; j < nBinCount; j++)
  {
    fCount = pHist[j];
    if(fCount > fMax)
      fMax = fCount;
  }

  if(fMax < Const_Almost_Zero)
    return;

  for(UINT j = 0; j < nBinCount; j++)
    pHist[j] /= fMax;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CString sPath(COutputEngine::get_full_path(pObj->get_filename()).c_str());
  CString sFileName = sPath + CString("Droplets_Diam_Distr.csv");

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (const char*)sFileName, (const char*)("w"));
  if(nErr != 0 || pStream == 0)
  {
    AfxMessageBox("Cannot open 'Droplets_Diam_Distr.csv' file for writing.");
    return;
  }

  fputs("D(mcm), Count\n\n", pStream);
  fD = fMinD + 0.5 * fStep;
  for(UINT j = 0; j < nBinCount; j++)
  {
    fprintf(pStream, "%7.4lf, %7.4lf\n", 1e+4 * fD, pHist[j]);
    fD += fStep;
  }

  fclose(pStream);
}

};  // namespace EvaporatingParticle
