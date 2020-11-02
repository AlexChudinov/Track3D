
#include "stdafx.h"
#include "ColorImage.h"
#include "ParticleTracking.h"
#include <algorithm>


namespace EvaporatingParticle
{

void RGBF_Color::clamp(float& rgb)
{
  if(rgb < 0.0f)
    rgb = 0.0f;
  if(rgb > 1.0f)
    rgb = 1.0f;
}

CColorImage::CColorImage()
  : m_bEnable(false)
{
}

CColorImage::~CColorImage()
{
  m_vRegNames.clear();
  m_vValues.clear();
  m_vColors.clear();
}

void CColorImage::set_default()
{
  m_bReady = false;

  m_nLevelsCount = 11;
  m_nColorMapType = cmRainbow;
  m_nScaleType = scLinear;

  m_bUserDefRange = false;
  m_fMinVal = 0;
  m_fMaxVal = 0;
}

CRegionsCollection CColorImage::get_regions() const
{
  CRegionsCollection vDrawnReg;
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CRegionsCollection& vRegColl = pObj->get_regions();
  size_t nRegCount = vRegColl.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    CRegion* pReg = vRegColl.at(i);
    if(std::find(m_vRegNames.begin(), m_vRegNames.end(), pReg->sName) == m_vRegNames.end())
      continue;

    vDrawnReg.push_back(pReg);
  }

  return vDrawnReg;
}

void CColorImage::get_values_array()
{
  m_vValues.clear();
  m_vValues.reserve(m_nLevelsCount);
  UINT nIntervCount = m_nLevelsCount + 1;
  if(m_nScaleType == scLinear)  // linear scale.
  {
    double fStep = (m_fMaxVal - m_fMinVal) / nIntervCount;
    if(fStep < Const_Almost_Zero)
      return;

    for(UINT i = 0; i < m_nLevelsCount; i++)
      m_vValues.push_back(m_fMinVal + (i + 1) * fStep);
  }
  else if(m_fMinVal > Const_Almost_Zero)  // logarithmic scale.
  {
    double fStep = log(m_fMaxVal / m_fMinVal) / nIntervCount;
    if(fStep < Const_Almost_Zero)
      return;

    for(UINT i = 0; i < m_nLevelsCount; i++)
      m_vValues.push_back(m_fMinVal * exp((i + 1) * fStep));
  }
}

void CColorImage::get_colors_array()
{
  double fXi, fKsi;
  unsigned char r, g, b;
  UINT nIntervCount = m_nLevelsCount + 1;

  m_vColors.clear();
  m_vColors.reserve(nIntervCount);

  for(UINT i = 0; i < nIntervCount; i++)
  {
    fKsi = (double)i / m_nLevelsCount, fXi;

    switch(m_nColorMapType)
    {
      case cmRainbow: // 4 intervals.
      {
        if(fKsi < 0.25)
        {
          fXi = 4 * fKsi;
          r = 0;
          g = clamp(255 * fXi);
          b = 255;
        }
        else if(fKsi < 0.5)
        {
          fXi = 4 * (fKsi - 0.25);
          r = 0;
          g = 255;
          b = clamp(255 * (1. - fXi));
        }
        else if(fKsi < 0.75)
        {
          fXi = 4 * (fKsi - 0.5);
          r = clamp(255 * fXi);
          g = 255;
          b = 0;
        }
        else
        {
          fXi = 4 * (fKsi - 0.75);
          r = 255;
          g = clamp(255 * (1. - fXi));
          b = 0;
        }
        break;
      }
      case cmRainbowExt:  // 5 intervals.
      {
        if(fKsi < 0.2)
        {
          fXi = 5 * fKsi;
          r = 0;
          g = clamp(255 * fXi);
          b = 255;
        }
        else if(fKsi < 0.4)
        {
          fXi = 5 * (fKsi - 0.2);
          r = 0;
          g = 255;
          b = clamp(255 * (1. - fXi));
        }
        else if(fKsi < 0.6)
        {
          fXi = 5 * (fKsi - 0.4);
          r = clamp(255 * fXi);
          g = 255;
          b = 0;
        }
        else if(fKsi < 0.8)
        {
          fXi = 5 * (fKsi - 0.6);
          r = 255;
          g = clamp(255 * (1. - fXi));
          b = 0;
        }
        else
        {
          fXi = 5 * (fKsi - 0.8);
          r = 255;
          g = 0;
          b = clamp(255 * fXi);
        }
        break;
      }
    }

    m_vColors.push_back(RGB_Color(r, g, b));
  }
}

RGB_Color CColorImage::get_color(double fVal) const
{
  if(m_vValues.size() == 0 || m_vColors.size() == 0)
    return RGB_Color();

  for(size_t i = 0; i < m_nLevelsCount; i++)
  {
    if(fVal < m_vValues.at(i))
      return m_vColors.at(i);
  }

  return m_vColors.at(m_nLevelsCount);
}

void CColorImage::save(CArchive& ar)
{
  const UINT nVersion = 1;
  ar << nVersion;

  ar << m_bEnable;

  ar << m_nScaleType;
  ar << m_nLevelsCount;
  ar << m_nColorMapType;

  CString cStr;
  size_t nCount = m_vRegNames.size();
  ar << nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    cStr = CString(m_vRegNames.at(i).c_str());
    ar << cStr;
  }

  ar << m_bUserDefRange;
  ar << m_fMinVal;
  ar << m_fMaxVal;
}

void CColorImage::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_bEnable;

  ar >> m_nScaleType;
  ar >> m_nLevelsCount;
  ar >> m_nColorMapType;

  CString cStr;
  size_t nCount;
  ar >> nCount;
  m_vRegNames.reserve(nCount);

  for(size_t i = 0; i < nCount; i++)
  {
    ar >> cStr;
    std::string sRegName((const char*)cStr);
    m_vRegNames.push_back(sRegName);
  }

  if(nVersion >= 1)
  {
    ar >> m_bUserDefRange;
    ar >> m_fMinVal;
    ar >> m_fMaxVal;
  }

  m_bReady = false;
}

void CColorImage::set_enable_user_range(bool bEnable)
{
  if(m_bUserDefRange != bEnable)
  {
    m_bUserDefRange = bEnable;
    if(m_bUserDefRange)
      restore_user_range();

    m_bReady = false;
  }
}

};  // namespace EvaporatingParticle