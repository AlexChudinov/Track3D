#ifndef _Matrix3D_
#define _Matrix3D_

#include "vector3d.hpp"

namespace EvaporatingParticle
{

//-----------------------------------------------------------------------------
//                       Matrix3D
//-----------------------------------------------------------------------------
class Matrix3D
{
  double  m00, m01, m02,
          m10, m11, m12,
          m20, m21, m22;

public:
  Matrix3D()
    : m00(1), m01(0), m02(0),
      m10(0), m11(1), m12(0),
      m20(0), m21(0), m22(1)
  {
  }

  Matrix3D(double a00, double a01, double a02,
           double a10, double a11, double a12,
           double a20, double a21, double a22)
    : m00(a00), m01(a01), m02(a02),
      m10(a10), m11(a11), m12(a12),
      m20(a20), m21(a21), m22(a22)
  {
  }

  Matrix3D(const Vector3D& v0, const Vector3D& v1, const Vector3D& v2)
    : m00(v0.x), m01(v1.x), m02(v2.x),
      m10(v0.y), m11(v1.y), m12(v2.y),
      m20(v0.z), m21(v1.z), m22(v2.z)
  {
  }

  Matrix3D(const Matrix3D& m)
    : m00(m.m00), m01(m.m01), m02(m.m02),
      m10(m.m10), m11(m.m11), m12(m.m12),
      m20(m.m20), m21(m.m21), m22(m.m22)
  {
  }

  void            set_elem(UINT nRow, UINT nCol, double fVal);

  Vector3D        get_col(UINT nCol) const;
  void            set_col(const Vector3D& v, UINT nCol);

  double          det() const;

  void            transpose();
  Matrix3D        transposed() const;

  static Matrix3D rot(const Vector3D& vDir, double fAngle, bool bRad = true);
  static Matrix3D euler(double fAlpha, double fBeta, double fGamma, bool bRad = true);

  Vector3D operator * (const Vector3D& v) const;
  Matrix3D operator * (const Matrix3D& m) const;
};

//-----------------------------------------------------------------------------
//                       CTransform
//-----------------------------------------------------------------------------
class CTransform
{
public:
  CTransform();

  enum
  {
    axisX     = 0,
    axisY     = 1,
    axisZ     = 2,
    axesCount = 3
  };

  void            set_default();

  bool            get_enable() const;
  DWORD_PTR       get_enable_ptr() const;
  void            set_enable(bool bEnable);

  int             get_rot_axis() const;
  DWORD_PTR       get_rot_axis_ptr() const;
  void            set_rot_axis(int nAxis);

  double          get_rot_angle() const;  // radian.
  DWORD_PTR       get_rot_angle_ptr() const;
  void            set_rot_angle(double fAngle);

  void            transform(Vector3D& v, bool bRadiusVector = true);

  void            save(CArchive& archive);
  void            load(CArchive& archive);

  const char*     axis_name(int nAxis) const;

protected:
  void            build_matrix();

  bool            m_bEnable;

  int             m_nAxis;
  double          m_fRotAngle;

private:
// Run-time variables:
  Vector3D        m_vRotAxis;
  Matrix3D        m_TransMtx;
  bool            m_bReady;
};

//-----------------------------------------------------------------------------
//                Inline functions
//-----------------------------------------------------------------------------
inline double Matrix3D::det() const
{
  return m00 * (m11 * m22 - m12 * m21) + m01 * (m12 * m20 - m10 * m22) + m02 * (m10 * m21 - m11 * m20);
}

inline Matrix3D Matrix3D::transposed() const
{
  Matrix3D m(*this);
  m.transpose();
  return m;
}

inline Vector3D Matrix3D::operator * (const Vector3D& v) const
{ 
  return Vector3D(m00 * v.x + m01 * v.y + m02 * v.z, m10 * v.x + m11 * v.y + m12 * v.z, m20 * v.x + m21 * v.y + m22 * v.z);
}

inline Matrix3D Matrix3D::operator * (const Matrix3D& m) const
{
  return Matrix3D(m00 * m.m00 + m01 * m.m10 + m02 * m.m20, m00 * m.m01 + m01 * m.m11 + m02 * m.m21, m00 * m.m02 + m01 * m.m12 + m02 * m.m22,
                  m10 * m.m00 + m11 * m.m10 + m12 * m.m20, m10 * m.m01 + m11 * m.m11 + m12 * m.m21, m10 * m.m02 + m11 * m.m12 + m12 * m.m22,
                  m20 * m.m00 + m21 * m.m10 + m22 * m.m20, m20 * m.m01 + m21 * m.m11 + m22 * m.m21, m20 * m.m02 + m21 * m.m12 + m22 * m.m22);
}

inline Vector3D Matrix3D::get_col(UINT nCol) const
{
  if(nCol > 2)
    return Vector3D();

  return nCol == 0 ? Vector3D(m00, m10, m20) : (nCol == 1 ? Vector3D(m01, m11, m21) : Vector3D(m02, m12, m22));
}


inline bool CTransform::get_enable() const
{
  return m_bEnable;
}

inline DWORD_PTR CTransform::get_enable_ptr() const
{
  return (DWORD_PTR)&m_bEnable;
}

inline void CTransform::set_enable(bool bEnable)
{
  if(m_bEnable != bEnable)
  {
    m_bEnable = bEnable;
    m_bReady = false;
  }
}

inline int CTransform::get_rot_axis() const
{
  return m_nAxis;
}

inline DWORD_PTR CTransform::get_rot_axis_ptr() const
{
  return (DWORD_PTR)&m_nAxis;
}

inline void CTransform::set_rot_axis(int nAxis)
{
  if(m_nAxis != nAxis)
  {
    m_nAxis = nAxis;
    m_bReady = false;
  }
}

inline double CTransform::get_rot_angle() const
{
  return m_fRotAngle;
}

inline DWORD_PTR CTransform::get_rot_angle_ptr() const
{
  return (DWORD_PTR)&m_fRotAngle;
}

inline void CTransform::set_rot_angle(double fAngle)
{
  if(m_fRotAngle != fAngle)
  {
    m_fRotAngle = fAngle;
    m_bReady = false;
  }
}

inline const char* CTransform::axis_name(int nAxis) const
{
  switch(nAxis)
  {
    case axisX: return _T("Global X");
    case axisY: return _T("Global Y");
    case axisZ: return _T("Global Z");
  }

  return _T("None");
}

};  // namespace EvaporatingParticle

#endif // #ifndef _Matrix3D_
