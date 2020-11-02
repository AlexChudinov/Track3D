#pragma once
#ifndef _PRIMITIVES_
#define _PRIMITIVES_

#include "AnsysMesh.h"
#include "vector2d.hpp"
#include "mathematics.h"


namespace EvaporatingParticle
{

//-------------------------------------------------------------------------------------------------
// An abstract class - a base class for all bounding primitives.
//-------------------------------------------------------------------------------------------------
class CPrimitive
{
public:
  enum
  {
    primBox     = 0,
    primCyl     = 1,
    primEllipse = 2,
    primCount   = 3
  };

  virtual Vector2D    model_coord(const Vector3D& vPos) const = 0;
  virtual bool        contains(const Vector3D& vPos) const = 0;
  virtual int         get_type() const = 0;

  virtual Vector3D    get_loc_normal(const Vector3D& vPos) const = 0;

  void                set_height(double fH);
  bool                is_ready() const;

  CNodesCollection*   get_reg_nodes();

protected:
  void                set_default();

  bool                get_regions(const CStringVector& vRegNames);
  bool                collect_reg_nodes();

  double              m_fHeight;

  CRegionsCollection  m_vRegions;
  CNodesCollection    m_vRegNodes;

  bool                m_bReady;
};

//-------------------------------------------------------------------------------------------------
// CPrimitiveBox - an arbitrarily oriented orthogonal box. 
//-------------------------------------------------------------------------------------------------
class CPrimitiveBox : public CPrimitive
{
public:
  CPrimitiveBox(const CStringVector& vRegNames, double fHeight);
  CPrimitiveBox(const Vector3D& vOrig, const Vector3D& vX, const Vector3D& vY, const Vector3D& vZ);

  virtual Vector2D    model_coord(const Vector3D& vPos) const;
  virtual bool        contains(const Vector3D& vPos) const;
  virtual int         get_type() const { return primBox; }

  virtual Vector3D    get_loc_normal(const Vector3D& vPos) const;

protected:
  void                calc_box_params();
  void                create_box(const Vector3D& vOrig, const Vector3D& vX, const Vector3D& vY, const Vector3D& vZ);
  void                calc_loc_triad();   // local Z is assumed to coincide with the world Z.

private:
  CPlane              m_Planes[6];

  Vector3D            m_vOrig,  // the box origin.
                      m_vLocX,  // m_vLocY * m_vLocZ.
                      m_vLocY,  // normalized average surface normal taken with opposite sign (the surface is assumed to be flat).
                      m_vLocZ;  // always (0, 0, 1).
};

class CMatrix_NxM;
//-------------------------------------------------------------------------------------------------
// CPrimCylSector - a sectorial volume between two cylinders.
//-------------------------------------------------------------------------------------------------
class CPrimCylSector : public CPrimitive
{
public:
  CPrimCylSector();
  CPrimCylSector(const CStringVector& vRegNames, double fX0, double fY0, double fR, double fHeight);

  virtual Vector2D    model_coord(const Vector3D& vPos) const;
  virtual bool        contains(const Vector3D& vPos) const;
  virtual int         get_type() const { return primCyl; }

  virtual Vector3D    get_loc_normal(const Vector3D& vPos) const;

  double              get_cyl_radius() const;
  Vector3D            get_orig_dir() const;
  Vector3D            get_end_dir() const;

protected:
/*
  void                calc_cyl_params();

  virtual void        calc_sums();
  virtual bool        calc_center_and_radius();

  void                grad_F(double x, double y, double r, double& dFdx, double& dFdy, double& dFdr) const;

  bool                get_initial_guess(double& a, double& b, double& r) const;

  double              initial_estimate_alpha(double x, double y, double r) const;

  bool                criterion(double dFdx, double dFdy, double dFdr) const;
*/

  void                assign(double a, double b, double r);
/*
// Test output:
  virtual FILE*       open_out_file();
  void                test_output(FILE* pStream, UINT nStep, double a, double b, double r, double grad);
*/
  double              m_fBigR,
                      m_fSmallR;

  CPlane              m_vSectBounds[2];   // two planes bounding the cylinder sector.

  CRay                m_vCylAxis;
  float               m_fEndZ;

// Normalized directions from the cylinder axis to the point on the cylinder surface with min(Z,[vP,vP0]) and max(Z,[vP,vP0]).
  Vector3D            m_vOrigDir,   // ... m_vOrigDir.z and m_vEndDir.z are always zero!
                      m_vEndDir;

// Run-time variables
/*
  double              m_fSX3,
                      m_fSY3,
                      m_fSXY2,
                      m_fSYX2,
                      m_fSX,
                      m_fSY,
                      m_fSXY,
                      m_fSX2,
                      m_fSY2;
*/
};

//-------------------------------------------------------------------------------------------------
// CEllipseData - contains coordinates of the initial point and semi-axes of an ellipse .
//-------------------------------------------------------------------------------------------------
struct CEllipseData
{
  CEllipseData()
    : x0(0), y0(0), a(1), b(1)
  {
  }

  CEllipseData(double _x0, double _y0, double _a, double _b)
    : x0(_x0), y0(_y0), a(_a), b(_b)
  {
  }

  double x0, y0, a, b;

  bool operator == (const CEllipseData& data) { return (x0 == data.x0) && (y0 = data.y0) && (a == data.a) && (b == data.b); }

  CEllipseData& operator += (const CEllipseData& data) { x0 += data.x0; y0 += data.y0; a += data.a; b += data.b; return *this; }
  CEllipseData& operator -= (const CEllipseData& data) { x0 -= data.x0; y0 -= data.y0; a -= data.a; b -= data.b; return *this; }

  CEllipseData& operator *= (double alpha) { x0 *= alpha; y0 *= alpha; a *= alpha; b *= alpha; return *this; }
  CEllipseData& operator /= (double alpha) { x0 /= alpha; y0 /= alpha; a /= alpha; b /= alpha; return *this; }

  CEllipseData operator + (const CEllipseData& data) const { return CEllipseData(x0 + data.x0, y0 + data.y0, a + data.a, b + data.b); }
  CEllipseData operator - (const CEllipseData& data) const { return CEllipseData(x0 - data.x0, y0 - data.y0, a - data.a, b - data.b); }

  CEllipseData operator * (double ksi) const { return CEllipseData(x0 * ksi, y0 * ksi, a * ksi, b * ksi); }
  CEllipseData operator / (double ksi) const { return CEllipseData(x0 / ksi, y0 / ksi, a / ksi, b / ksi); }

// Direct multiplication:
  CEllipseData operator * (const CEllipseData& data) const { return CEllipseData(x0 * data.x0, y0 * data.y0, a * data.a, b * data.b); }
};

typedef std::vector<Vector3F> CNodeLocVector;
//-------------------------------------------------------------------------------------------------
// CEllipticalCylSector - a sectorial volume between two ellipses.
//-------------------------------------------------------------------------------------------------
class CEllipticalCylSector : public CPrimCylSector
{
public:
  CEllipticalCylSector(const CStringVector& vRegNames, const CEllipseData& data, double fHeight, double fDelta = 1e-4);
  
  virtual int         get_type() const { return primEllipse; }

  virtual bool        contains(const Vector3D& vPos) const;

  virtual Vector2D    model_coord(const Vector3D& vPos) const;
  virtual Vector3D    get_loc_normal(const Vector3D& vPos) const;

  double              get_h(const Vector3D& vPos) const;  // returns approximate value of the height of vPos over the base surface (ellipse).

protected:
  virtual void        calc_sums();
  virtual bool        calc_center_and_radius();

//  void                collect_opt_nodes(double fdZ = 0.05); // collect nodes in a narrow stripe -fdZ < z < +fdZ for ellipse parameters calculation.

// Square of the discrepancy.
//  double              get_F(const CEllipseData& data) const;
//  double              get_F(double x0, double y0, double a, double b) const;

//  bool                grad2_F(const CEllipseData& data, CMatrix_NxM& mtx) const;
//  void                grad_F(const CEllipseData& data, CEllipseData& grad) const;

// Both these functions return the number of steps done. If zero is returned, then something wrong occured.
//  UINT                one_iter_Newton(const CEllipseData& start, CEllipseData& end, double& fMinF) const;
//  UINT                one_iter_steepest_descent(const CEllipseData& start, CEllipseData& end, double& fMinF) const;

//  UINT                min_along_line(const CEllipseData& beg, CEllipseData& end, double& fMinF) const;

//  bool                get_initial_guess(CEllipseData& data) const;

//  bool                criterion(double fEps2) const;

  void                assign(const CEllipseData& data);
//  bool                get_solution(const CMatrix_NxM& mtx, CEllipseData& data) const;

  bool                prepare_ellipse();  // initilize the ellipse length spline.
  bool                get_x_ranges(double& fXmin, double& fXmax) const;

  double              ellipse(double x, double y, double a, double b) const;  // computes [(x-x0)/a]^2 + [(y-y0)/b]^2.
  double              ellipse_arc_length(double x1, double x2) const; // numerically computes the length of the elliptical arc, for spline initialization.

  static double       ellipse_elem_arc_length(double alpha, double* par);

  void                clamp(double& x) const;

// Test output:
//  void                test_output(FILE* pStream, UINT nStep, const CEllipseData& data, double grad);

 // Run-time variables
//  CNodeLocVector      m_vOptNodes;  // nodes used for ellipse parameters determination only.
  CubSpline1D         m_LenSpl;     // length along the ellipse L(x) beginning from the minimal x and ending with the maximal x.

  double              m_fA,   // x and y semi-axes, respectively.
                      m_fB;

  double              m_fD;   // shift for all the variables during derivatives calculation.
};


//-------------------------------------------------------------------------------------------------
// Inlines:
//-------------------------------------------------------------------------------------------------
inline void CPrimitive::set_height(double fH)
{
  m_fHeight = fH;
}

inline bool CPrimitive::is_ready() const
{
  return m_bReady;
}

inline CNodesCollection* CPrimitive::get_reg_nodes()
{
  return &m_vRegNodes;
}

inline double CPrimCylSector::get_cyl_radius() const
{
  return m_fSmallR;
}

inline Vector3D CPrimCylSector::get_orig_dir() const
{
  return m_vOrigDir;
}

inline Vector3D CPrimCylSector::get_end_dir() const
{
  return m_vEndDir;
}

inline void CEllipticalCylSector::clamp(double& x) const
{
  if(x < m_LenSpl.get_x(0))
    x = m_LenSpl.get_x(0);

  int N = m_LenSpl.get_nodes_count() - 1;
  if(x > m_LenSpl.get_x(N))
    x = m_LenSpl.get_x(N);
}

};  // namespace EvaporatingParticle

#endif
