
#include "stdafx.h"

#include "matrix3d.hpp"
#include "constant.hpp"
#include <math.h>

namespace EvaporatingParticle
{

//-----------------------------------------------------------------------------
//                       Matrix3D
//-----------------------------------------------------------------------------
Matrix3D Matrix3D::rot(const Vector3D& vDir, double fAngle, bool bRad)
{
  Vector3D v = vDir.normalized();

  double fArg = bRad ? fAngle : Const_PI * fAngle / 180.;
  double fSinArg = sin(fArg);
  double fCosArg = cos(fArg);
  double fOneMinusCosArg = 1.- fCosArg;
  Vector3D vDirSinArg(v.x * fSinArg, v.y * fSinArg, v.z * fSinArg);

  double m00 = fCosArg + fOneMinusCosArg * v.x * v.x;
  double m01 = fOneMinusCosArg * v.x * v.y;
  double m02 = fOneMinusCosArg * v.x * v.z;

  double m10 = m01;
  double m11 = fCosArg + fOneMinusCosArg * v.y * v.y;
  double m12 = fOneMinusCosArg * v.y * v.z;

  double m20 = m02;
  double m21 = m12;
  double m22 = fCosArg + fOneMinusCosArg * v.z * v.z;

  m01 -= vDirSinArg.z;
  m02 += vDirSinArg.y;

  m10 += vDirSinArg.z;
  m12 -= vDirSinArg.x;

  m20 -= vDirSinArg.y;
  m21 += vDirSinArg.x;

  return Matrix3D(m00, m01, m02,
                  m10, m11, m12,
                  m20, m21, m22);
}

Matrix3D Matrix3D::euler(double fAlpha, double fBeta, double fGamma, bool bRad)
{
  double fArgAlp = bRad ? fAlpha : Const_PI * fAlpha / 180.;
  double fSinAlp = sin(fArgAlp);
  double fCosAlp = cos(fArgAlp);

  double fArgBet = bRad ? fBeta : Const_PI * fBeta / 180.;
  double fSinBet = sin(fArgBet);
  double fCosBet = cos(fArgBet);

  double fArgGam = bRad ? fGamma : Const_PI * fGamma / 180.;
  double fSinGam = sin(fArgGam);
  double fCosGam = cos(fArgGam);

  double m00 =  fCosAlp * fCosGam - fSinAlp * fCosBet * fSinGam;
  double m01 = -fCosAlp * fSinGam - fSinAlp * fCosBet * fCosGam;
  double m02 =  fSinAlp * fSinBet;

  double m10 =  fSinAlp * fCosGam + fCosAlp * fCosBet * fSinGam;
  double m11 = -fSinAlp * fSinGam + fCosAlp * fCosBet * fCosGam;
  double m12 = -fCosAlp * fSinBet;

  double m20 =  fSinBet * fSinGam;
  double m21 =  fSinBet * fCosGam;
  double m22 =  fCosBet;

  return Matrix3D(m00, m01, m02,
                  m10, m11, m12,
                  m20, m21, m22);
}

void Matrix3D::set_elem(UINT nRow, UINT nCol, double fVal)
{
  if(nRow > 2 || nCol > 2)
    return;

  switch(nRow)
  {
    case 0: nCol == 0 ? m00 = fVal : (nCol == 1 ? m01 = fVal : m02 = fVal); break;
    case 1: nCol == 0 ? m10 = fVal : (nCol == 1 ? m11 = fVal : m12 = fVal); break;
    case 2: nCol == 0 ? m20 = fVal : (nCol == 1 ? m21 = fVal : m22 = fVal); break;
  }
}

void Matrix3D::set_col(const Vector3D& v, UINT nCol)
{
  switch(nCol)
  {
    case 0: m00 = v.x; m10 = v.y; m20 = v.z; break;
    case 1: m01 = v.x; m11 = v.y; m21 = v.z; break;
    case 2: m02 = v.x; m12 = v.y; m22 = v.z; break;
  }
}

void Matrix3D::transpose()
{
  double a01 = m01;
  double a02 = m02;
  double a12 = m12;

  m01 = m10;
  m02 = m20;
  m12 = m21;

  m10 = a01;
  m20 = a02;
  m21 = a12;
}

//-----------------------------------------------------------------------------
//                       CTransform
//-----------------------------------------------------------------------------
CTransform::CTransform()
  : m_bReady(false)
{
  set_default();
}

void CTransform::set_default()
{
  m_bEnable = false;
  m_nAxis = CTransform::axisX;
  m_fRotAngle = 0;
}

void CTransform::transform(Vector3D& v, bool bRadiusVector)
{
  if(!m_bReady)
    build_matrix();

  v = m_TransMtx * v;

/*
  if(bRadiusVector)
// TO DO: apply shift.
*/
}

void CTransform::build_matrix()
{
  switch(m_nAxis)
  {
    case axisX: m_vRotAxis = Vector3D(1, 0, 0); break;
    case axisY: m_vRotAxis = Vector3D(0, 1, 0); break;
    case axisZ: m_vRotAxis = Vector3D(0, 0, 1); break;
  }

  m_TransMtx = Matrix3D::rot(m_vRotAxis, m_fRotAngle);
  m_bReady = true;
}

void CTransform::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  ar << m_bEnable;
  ar << m_nAxis;
  ar << m_fRotAngle;
}

void CTransform::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_bEnable;
  ar >> m_nAxis;
  ar >> m_fRotAngle;
}

}; // namespace EvaporatingParticle