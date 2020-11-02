#ifndef _Mathematics_
#define _Mathematics_

#pragma once

#include "vector3D.hpp"

namespace EvaporatingParticle
{

class CMath
{
public:
  static double root(double a, double b, double(*f)(double x, double* p), double eps, double* par);

// Returns exp(-x) * sqrt(x).
  static double energy_probability_func(double x);

// Computes the integral from E0 to +infinity of the dWe = (2 / sqrt(pi * kT^3)) * exp(-E / kT) * sqrt(E) * dE;
// the probability for the atom to have the kinetic energy in the interval from E to E+dE. Landau, v.5, page 104.
// Input (dimensionless): Threshold energy E0/kT.
  static double energy_probability_integral(double fE0_ovr_kT);

// Transformation of radius-vector pos to the cylinder c.s. and back. If bCylAxisX == false then the cylinder and the cartesian Z-axes coincide.
  static bool   cart_to_cyl(const Vector3D& vPos, double& z, double& r, double& azim, bool bCylAxisX = true);
  static void   cyl_to_cart(double z, double r, double azim, Vector3D& vPos, bool bCylAxisX = true);

  static double erf(double x);

// Numeric integration of f(x, par) dx from a to b using the trapeze method. N is the number of points (not intervals) on the [a, b] line segment.
  static double integr(double a, double b, double(*f)(double x, double* p), UINT N, double* par);
};


//--------------------------------------------------------------------------------------------------
// CMatrix_NxM - an auxiliary class for Gauss elimination solver
//--------------------------------------------------------------------------------------------------
class CMatrix_NxM
{
public:
  CMatrix_NxM(UINT nRowCount, UINT nColCount);
  ~CMatrix_NxM();

  double      get(UINT nRow, UINT nCol) const;
  bool        set(UINT nRow, UINT nCol, double fVal);

  UINT        get_row_count() const;
  UINT        get_col_count() const;

  bool        permute_rows(UINT nRow1, UINT nRow2);

  bool        get_abs_max_row(UINT nDiag, UINT& nRowMax) const;

  bool        is_good() const;

protected:
  void        allocate_arrays();
  void        delete_arrays();

private:
  double**    m_pData;
  UINT        m_nRowCount,
              m_nColCount;
};

//--------------------------------------------------------------------------------------------------
// CGaussEliminationSolver
//--------------------------------------------------------------------------------------------------
class CGaussEliminationSolver
{
public:
  CGaussEliminationSolver(UINT nRowCount, UINT nColCount)
    : m_Mtx(nRowCount, nColCount)
  {
  }

  enum  // result flags.
  {
    resOK           = 0,
    resMtxDegen     = 1,
    resMtxBadStruct = 2
  };

  int           solve();  // returns one of the result flags.

  CMatrix_NxM&  get_matrix();

  void          output_mtx(const char* pFileName);

protected:
// Initial, pre-solve check of the matrix. Normally, the number of columns is greater than the number of rows.
// Additional columns are interpreted as different right parts, for which solutions are to be found.
  bool          bad_matrix_structure();

  bool          run_forward();
  bool          run_backward();

private:
  CMatrix_NxM   m_Mtx;

// Run-time variable:
  int           m_nRes;
};

//--------------------------------------------------------------------------------------------------
// CMatrix_NxM inlines
//--------------------------------------------------------------------------------------------------
inline double CMatrix_NxM::get(UINT nRow, UINT nCol) const
{
  return ((m_pData != NULL) && (nRow < m_nRowCount) && (nCol < m_nColCount)) ? m_pData[nRow][nCol] : -FLT_MAX;
}

inline UINT CMatrix_NxM::get_row_count() const
{
  return m_nRowCount;
}

inline UINT CMatrix_NxM::get_col_count() const
{
  return m_nColCount;
}

inline bool CMatrix_NxM::is_good() const
{
  return (m_pData != NULL) && (m_nRowCount > 0) && (m_nColCount > 0);
}

//--------------------------------------------------------------------------------------------------
// CGaussEliminationSolver inlines
//--------------------------------------------------------------------------------------------------
inline bool CGaussEliminationSolver::bad_matrix_structure()
{
  return m_Mtx.get_col_count() <= m_Mtx.get_row_count();
}

inline CMatrix_NxM& CGaussEliminationSolver::get_matrix()
{
  return m_Mtx;
}

//-----------------------------------------------------------------------------
//  CSplineInfo
//-----------------------------------------------------------------------------
struct CSplineInfo  // an auxiliary structure for 2D and 3D spline coefficients building.
{
  CSplineInfo() { nStage = nBoundType = i = j = k = 0; }
  int nStage,     // a step of spline building; for Spline2D this can be stepVx, stepVy or stepVxy.
      nBoundType, // different boundary types are supported.
      i,  //
      j,  // Corresponding fixed indices, defining which column, row, etc is being processed.
      k;  //
};

//-----------------------------------------------------------------------------
//  CubicSplineNode
//-----------------------------------------------------------------------------
struct CubicSplineNode
{
  CubicSplineNode() { V = 0; }
  double V;   // the given value at this node.
};

typedef std::vector<double> CPointVector;
typedef std::vector<CubicSplineNode*> CSplineNodeColl;
//------------------------------------------------------------------------------
// CubicSpline - an abstract class for cubic splines family with different boundary
// conditions. The computed coefficients s(i) are, in fact, 1/6 of the second
// derivative of the spline at i-th point. See J. Forsight, M. Malcolm, C. Moler,
// "Computer Methods of Mathematical Computations", 1980, p.86 of Russian edition.
//------------------------------------------------------------------------------
class CubicSpline
{
public:
  CubicSpline();
  virtual ~CubicSpline();

  enum  // types of boundary conditions.
  {
    btFirstDeriv  = 0,
    btSecondDeriv = 1,
    btPeriodic    = 2,  // for periodic functions.
    btFree        = 3   // the same as btSecondDeriv, but with zero second derivatives at both ends.
  };

  double            get_node_val(int nNodeIndex) const;
  void              set_node_val(int nNodeIndex, double fValue);

  double            get_x(int i) const;
  void              set_x(int i, double x);

  void              get_range(double& fMinVal, double& fMaxVal) const;
  
protected:
// First derivatives at the ends of the X-interval.
  double            m_fVx0,
                    m_fVxN;
// Second derivatives at the ends of the X-interval.
  double            m_fVxx0,
                    m_fVxxN;

  double            m_fMinX,
                    m_fMaxX,
                    m_fRangeX;

  CPointVector      m_vX;
  CSplineNodeColl   m_vSplData;

  bool              m_bReady;

  virtual bool      check() = 0;
  virtual void      prepare() = 0;
  virtual void      set_default() = 0;

  virtual int       index(int i, int j, int k) = 0;
  
  virtual void      allocate_arrays() = 0;
  virtual void      delete_arrays();

// Computation of spline coefficients for a periodic spline.
  void              periodic_coeff(double* pAux, int nIntervCount, const CSplineInfo& info);

// Computation of spline coefficients for a non-periodic spline.
  void              coeff(double* pAux, int nIntervCount, const CSplineInfo& info);

// Re-calculate sigma-coefficients to the final spline coefficients.
  virtual void      assign(double* dx, double* dv, double* s, int nIntervCount, const CSplineInfo& info) = 0;

  int               interval(double x, const CPointVector& v, int& nInterv, bool bEven) const;
  int               even_interval(double x, const CPointVector& v) const;

  int               divide_by_two(double x, const CPointVector& v) const;

// Careful: no bounds check is performed here!
  bool              inside(double x, const CPointVector& v, int nInterv) const;

  void              set_equidistant(double fMin, double fMax, CPointVector& v);
};

//------------------------------------------------------------------------------
// Inlines of CubicSpline
//------------------------------------------------------------------------------
inline double CubicSpline::get_node_val(int nIndex) const
{
  return (nIndex >= 0) && (nIndex < m_vSplData.size()) ? m_vSplData[nIndex]->V : FLT_MAX;
}

inline double CubicSpline::get_x(int i) const
{
  return (i >= 0) && (i < m_vX.size()) ? m_vX[i] : FLT_MAX;
}

//------------------------------------------------------------------------------
//  CubSpline1DNode
//------------------------------------------------------------------------------
struct CubSpline1DNode : public CubicSplineNode
{
  CubSpline1DNode() { b = c = d = 0; }
  double b,   // S(x) = V + b*(x-x(i)) + c*(x-x(i))^2 + d*(x-x(i))^3 and
         c,   // is computed by the nest multiplication method, see the code.
         d;
};

//------------------------------------------------------------------------------
//  CubSpline1D - a simple cubic spline with different boundary conditions.
//------------------------------------------------------------------------------
class CubSpline1D : public CubicSpline
{
public:
  CubSpline1D(int nBoundType = CubicSpline::btFree);
  virtual ~CubSpline1D();

//------------------------------------------------------------------------------
// A usual sequence of operations needed to initialize the spline:
//  1. set_boundary_type(int nBoundType);
//  2. set_nodes_count(int nNodesCount);
//  3. set_range(double fXmin, double fXmax, bool bEven);
//  4. for(int i = 0; i < nNodesCount; i++)
//       set_node_val(i, double fVal);
//  5. prepare();
//------------------------------------------------------------------------------
// Note: the user must call this function himself before using the spline!
  virtual void      prepare();

// CubSpline1D specific:
  double            get(double x, int nOrder, int& nInterv) const;

  double            integral(double x1, double x2) const;

  int               get_nodes_count() const;
  void              set_nodes_count(int nCount);

  bool              get_even() const;
  void              set_even(bool bEven = true);

  void              set_range(double fMinX, double fMaxX, bool bEven = true);

  void              set_boundary_type(int nType);

// Boundary conditions for btFirstDeriv and btSecondDeriv types. How the input values are
// interpreted depends on m_nBoundTypeX value: they can be either first or second derivatives.
  void              set_bound_derivatives(double fV0, double fVN);

// Attention! The following two functions will have no effect if m_bEven == true!
  
// Returns false and does nothing if this node cannot be inserted. For example,
// if (fX, fValue) is too close to an existing node. Attention! The existing
// array of x is supposed to be ordered!
  bool              insert_node(double fX, double fValue);

// Returns false and does nothing if fX <= the last existing x. Normally, this 
// function adds the (fX, fValue) node to the end of the spline.
  bool              push_back_node(double fX, double fValue);
  
protected:
  bool              m_bEven;
  int               m_nBoundType;

  int               m_nNodesCount;
  
  virtual bool      check();
  virtual void      set_default();
  virtual void      allocate_arrays();

// Re-calculate sigma-coefficients to the final spline coefficients.
  virtual void      assign(double* dx, double* dv, double* s, int nIntervCount, const CSplineInfo& info);

  virtual int       index(int i, int j, int k);
};

//------------------------------------------------------------------------------
// Inlines of CubSpline1D
//------------------------------------------------------------------------------
inline void CubSpline1D::set_boundary_type(int nType)
{
  if(m_nBoundType != nType)
  {
    m_nBoundType = nType;
    m_bReady = false;
  }
}

inline int CubSpline1D::get_nodes_count() const
{
  return m_vX.size();
}

inline bool CubSpline1D::get_even() const
{
   return m_bEven;
}

inline void CubSpline1D::set_even(bool bEven)
{
  if(m_bEven != bEven)
  {
    m_bEven = bEven;
    m_bReady = false;
  }
}

};  // namespace EvaporatingParticle

#endif // _Mathematics_