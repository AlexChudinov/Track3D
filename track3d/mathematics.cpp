
#include "stdafx.h"

#include "math.h"
#include "mathematics.h"
#include "constant.hpp"
#include "vector3d.hpp"


namespace EvaporatingParticle
{

#define sign(a) (a > 0.0 ? 1 : -1)
double CMath::root(double a, double b, double(*f)(double x, double* p), double eps, double* par)
{
  double xa = a;
  double xb = b;
  double fa = f(xa, par);
  double fb = f(xb, par);
  if(sign(fa) == sign(fb))
    return xa;

  UINT nIterCount = 0;
  double xc, fc, dx, yeps;
  while(nIterCount < 10)
  {
// First two iterations do by the bisection method, the next ones - by secant method.
    xc = nIterCount < 2 ? 0.5 * (xa + xb) : (xa * fb - xb * fa) / (fb - fa);
    fc = f(xc, par);
    if(sign(fc) == sign(fb))
    {
      xb = xc;
      fb = fc;
    }
    else
    {
      xa = xc;
      fa = fc;
    }

    dx = xb - xa;
    if(fabs(dx) < eps)
      return 0.5 * (xa + xb);

    yeps = eps * fabs((fb - fa) / dx);
    if(fabs(fc) < yeps)
      return xc;

    nIterCount++;
  }

  return 0.5 * (xa + xb);
}

double CMath::energy_probability_func(double x)
{
  return exp(-x) * sqrt(x);
}

static const double scfTwoOvrSqrt2 = 2. / sqrt(Const_PI);
//---------------------------------------------------------------------------------------
//                     x0
// Tabulated values of I exp(-x)*sqrt(x)*dx for equidistant x0 from 0 to 15 with a step dx = 0.1.
//                     0
//---------------------------------------------------------------------------------------
static const double scFunc[151] = {
  0.0000000e+000, 1.9662873e-002, 5.2755656e-002, 9.1579042e-002, 1.3319962e-001, 1.7592824e-001, 2.1868601e-001, 2.6075443e-001, 3.0164903e-001,
  3.4104597e-001, 3.7873553e-001, 4.1459074e-001, 4.4854540e-001, 4.8057815e-001, 5.1070070e-001, 5.3894900e-001, 5.6537647e-001, 5.9004890e-001,
  6.1304042e-001, 6.3443050e-001, 6.5430155e-001, 6.7273713e-001, 6.8982055e-001, 7.0563383e-001, 7.2025693e-001, 7.3376717e-001, 7.4623883e-001,
  7.5774289e-001, 7.6834684e-001, 7.7811465e-001, 7.8710671e-001, 7.9537986e-001, 8.0298748e-001, 8.0997959e-001, 8.1640293e-001, 8.2230115e-001,
  8.2771489e-001, 8.3268198e-001, 8.3723756e-001, 8.4141423e-001, 8.4524221e-001, 8.4874950e-001, 8.5196197e-001, 8.5490356e-001, 8.5759635e-001,
  8.6006075e-001, 8.6231555e-001, 8.6437807e-001, 8.6626429e-001, 8.6798889e-001, 8.6956537e-001, 8.7100618e-001, 8.7232272e-001, 8.7352548e-001,
  8.7462411e-001, 8.7562744e-001, 8.7654358e-001, 8.7737997e-001, 8.7814344e-001, 8.7884024e-001, 8.7947609e-001, 8.8005625e-001, 8.8058553e-001,
  8.8106831e-001, 8.8150863e-001, 8.8191018e-001, 8.8227632e-001, 8.8261014e-001, 8.8291445e-001, 8.8319184e-001, 8.8344465e-001, 8.8367505e-001,
  8.8388499e-001, 8.8407628e-001, 8.8425056e-001, 8.8440932e-001, 8.8455393e-001, 8.8468565e-001, 8.8480561e-001, 8.8491485e-001, 8.8501432e-001,
  8.8510489e-001, 8.8518735e-001, 8.8526242e-001, 8.8533075e-001, 8.8539296e-001, 8.8544957e-001, 8.8550110e-001, 8.8554799e-001, 8.8559066e-001,
  8.8562948e-001, 8.8566481e-001, 8.8569696e-001, 8.8572620e-001, 8.8575280e-001, 8.8577700e-001, 8.8579901e-001, 8.8581903e-001, 8.8583724e-001,
  8.8585380e-001, 8.8586886e-001, 8.8588256e-001, 8.8589501e-001, 8.8590634e-001, 8.8591663e-001, 8.8592600e-001, 8.8593451e-001, 8.8594224e-001,
  8.8594928e-001, 8.8595567e-001, 8.8596148e-001, 8.8596677e-001, 8.8597157e-001, 8.8597593e-001, 8.8597990e-001, 8.8598351e-001, 8.8598678e-001,
  8.8598976e-001, 8.8599246e-001, 8.8599492e-001, 8.8599716e-001, 8.8599919e-001, 8.8600103e-001, 8.8600271e-001, 8.8600423e-001, 8.8600561e-001,
  8.8600687e-001, 8.8600801e-001, 8.8600904e-001, 8.8600999e-001, 8.8601084e-001, 8.8601162e-001, 8.8601232e-001, 8.8601297e-001, 8.8601355e-001,
  8.8601408e-001, 8.8601456e-001, 8.8601499e-001, 8.8601539e-001, 8.8601575e-001, 8.8601608e-001, 8.8601637e-001, 8.8601664e-001, 8.8601689e-001,
  8.8601711e-001, 8.8601731e-001, 8.8601749e-001, 8.8601766e-001, 8.8601781e-001, 8.8601795e-001, 8.8601807e-001 };

//---------------------------------------------------------------------------------------
// Using formula (860.04) from H.B.Dwight, "Tables of integrals and other mathematical data" 
// the infinite integral can be reduced to an integral with finite limits from 0 to x0 = E0/kT.
//                                inf                                      x0
// Thus, instead of computation of I exp(-x)*sqrt(x)*dx we have to compute I exp(-x)*sqrt(x)*dx.
//                                 x0                                      0
//---------------------------------------------------------------------------------------
double CMath::energy_probability_integral(double fE0_ovr_kT)
{
  if(fE0_ovr_kT >= 15)
    return 0;

  UINT i0 = 10 * fE0_ovr_kT;
  UINT i1 = i0 + 1;

  double x0 = 0.1 * i0;
  double ksi = fE0_ovr_kT - x0;
  double r = (1. - ksi) * scFunc[i0] + ksi * scFunc[i1]; 

  return 1. - scfTwoOvrSqrt2 * r;
}

bool CMath::cart_to_cyl(const Vector3D& vPos, double& z, double& r, double& azim, bool bCylAxisX)
{
  Vector3D p = bCylAxisX ? Vector3D(vPos.y, vPos.z, vPos.x) : vPos;
  z = p.z;
  r = sqrt(p.x * p.x + p.y * p.y);
  if(r < Const_Almost_Zero)
    return false;

  if(p.y > Const_Almost_Zero)
    azim = acos(p.x / r);
  else if(p.y < -Const_Almost_Zero)
    azim = Const_2PI - acos(p.x / r);
  else if(p.x < 0)
    azim = Const_PI;
  else
    azim = 0;

  return true;
}

void CMath::cyl_to_cart(double z, double r, double azim, Vector3D& vPos, bool bCylAxisX)
{
  Vector3D p(r * cos(azim), r * sin(azim), z);
  vPos = bCylAxisX ? Vector3D(p.z, p.x, p.y) : p;
}

double CMath::erf(double ax)
{
  const double a1 = 0.254829592;
  const double a2 = -0.284496736;
  const double a3 = 1.421413741;
  const double a4 = -1.453152027;
  const double a5 = 1.061405429;
  const double p = 0.3275911;

  double x = fabs(ax);
  double T = 1. / (1. + p * x);
  double res = 1. - (((((a5 * T + a4) * T) + a3) * T + a2) * T + a1) * T * exp(-x * x);
  if(ax < 0)
    res = -res;

  return res;
}

double CMath::integr(double a, double b, double(*f)(double x, double* p), UINT N, double* par)
{
  if(N <= 1)
    return 0;

  double fX = a;
  double fdX = (b - a) / (N - 1);

  double fRes = 0.5 * f(fX, par);
  for(UINT i = 1; i < N - 1; i++)
  {
    fX += fdX;
    fRes += f(fX, par);
  }

  fX = b;
  fRes += 0.5 * f(fX, par);

  fRes *= fdX;
  return fRes;
}

//--------------------------------------------------------------------------------------------------
// CMatrix_NxM - an auxiliary class for Gauss elimination solver
//--------------------------------------------------------------------------------------------------
CMatrix_NxM::CMatrix_NxM(UINT nRowCount, UINT nColCount)
  : m_pData(NULL), m_nRowCount(nRowCount), m_nColCount(nColCount)
{
  allocate_arrays();
}

CMatrix_NxM::~CMatrix_NxM()
{
  delete_arrays();
}

bool CMatrix_NxM::set(UINT nRow, UINT nCol, double fVal)
{
  if((m_pData == NULL) || (nRow >= m_nRowCount) || (nCol >= m_nColCount))
    return false;

  m_pData[nRow][nCol] = fVal;
  return true;
}

bool CMatrix_NxM::permute_rows(UINT nRow1, UINT nRow2)
{
  if((m_pData == NULL) || (nRow1 >= m_nRowCount) || (nRow2 >= m_nRowCount))
    return false;

  if(nRow1 == nRow2)
    return true;

  double* pTmp = m_pData[nRow1];
  m_pData[nRow1] = m_pData[nRow2];
  m_pData[nRow2] = pTmp;
  return true;
}

bool CMatrix_NxM::get_abs_max_row(UINT nDiag, UINT& nRowMax) const
{
  if((m_pData == NULL) || (nDiag >= m_nRowCount) || (nDiag >= m_nColCount))
    return false;

  nRowMax = nDiag;
  double fAbsVal, fAbsMax = fabs(m_pData[nRowMax][nDiag]);
  for(UINT i = nDiag + 1; i < m_nRowCount; i++)
  {
    fAbsVal = fabs(m_pData[i][nDiag]);
    if(fAbsVal > fAbsMax)
    {
      fAbsMax = fAbsVal;
      nRowMax = i;
    }
  }

  return true;
}

void CMatrix_NxM::allocate_arrays()
{
  delete_arrays();
  if(m_nRowCount < 1 || m_nColCount < 1)
    return;

  m_pData = new double*[m_nRowCount];
  for(UINT i = 0; i < m_nRowCount; i++)
    m_pData[i] = new double[m_nColCount];

  for(UINT i = 0; i < m_nRowCount; i++)
    for(UINT j = 0; j < m_nColCount; j++)
      m_pData[i][j] = (i == j) ? 1.0 : 0.0;
}

void CMatrix_NxM::delete_arrays()
{
  if(m_pData == NULL)
    return;

  for(UINT i = 0; i < m_nRowCount; i++)
    delete[] m_pData[i];

  delete[] m_pData;
}

//--------------------------------------------------------------------------------------------------
// CGaussEliminationSolver
//--------------------------------------------------------------------------------------------------
int CGaussEliminationSolver::solve()
{
  if(run_forward())
    run_backward();

  return m_nRes;
}

// Reducing the matrix to the upper-triangle form.
bool CGaussEliminationSolver::run_forward()
{
  UINT N = m_Mtx.get_row_count();
  UINT M = m_Mtx.get_col_count();

  UINT nAbsMaxRow;
  double fMainElem, fCoeff, fElem;
  for(UINT nRow = 0; nRow < N - 1; nRow++)
  {
    m_Mtx.get_abs_max_row(nRow, nAbsMaxRow);
    m_Mtx.permute_rows(nRow, nAbsMaxRow);

    fMainElem = m_Mtx.get(nRow, nRow);
    if(fabs(fMainElem) < Const_Almost_Zero)
    {
      m_nRes = resMtxDegen;
      return false;
    }

    for(UINT i = nRow + 1; i < N; i++)
    {
      fCoeff = -m_Mtx.get(i, nRow) / fMainElem;
      for(UINT j = nRow; j < M; j++)
      {
        fElem = m_Mtx.get(nRow, j) * fCoeff + m_Mtx.get(i, j);
        m_Mtx.set(i, j, fElem);
      }
    }
  }

// Final check:
  for(UINT k = 0; k < N; k++)
  {
    if(fabs(m_Mtx.get(k, k)) < Const_Almost_Zero)
    {
      m_nRes = resMtxDegen;
      return false;
    }
  }

  m_nRes = resOK;
  return true;
}

// Solving the upper-triangle system. After this function the right-most columns of m_Mtx will contain 
// solution vectors for corresponding right parts. Kahaner, Moler, Nash "Numerical Methods and Software", 
// Moscow, Mir, 1998, p.62.
bool CGaussEliminationSolver::run_backward()
{
  UINT N = m_Mtx.get_row_count();
  UINT M = m_Mtx.get_col_count();

  double fSum, fX;
  for(int k = N; k < M; k++) // loop over all right part columns of this system.
  {
    for(int i = N - 1; i >= 0; i--)
    {
      fSum = 0;
      for(int j = i + 1; j < N; j++)
        fSum += m_Mtx.get(i, j) * m_Mtx.get(j, k);

      fX = (m_Mtx.get(i, k) - fSum) / m_Mtx.get(i, i);
      m_Mtx.set(i, k, fX);
    }
  }

  return true;
}

// DEBUG
void CGaussEliminationSolver::output_mtx(const char* pFileName)
{
  std::string cFileName(pFileName);

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return;

  UINT N = m_Mtx.get_row_count();
  UINT M = m_Mtx.get_col_count();

  for(UINT i = 0; i < N; i++)
  {
    for(UINT j = 0; j < M - 1; j++)
      fprintf_s(pStream, "%12.3e", m_Mtx.get(i, j));

    fprintf_s(pStream, "%12.3e\n\n", m_Mtx.get(i, M - 1));
  }

  fclose(pStream);
}
// END DEBUG

//------------------------------------------------------------------------------
//  CubicSpline - an abstract class for cubic splines family with different
//                boundary conditions.
//------------------------------------------------------------------------------
CubicSpline::CubicSpline()
  : m_bReady(false)
{
}

CubicSpline::~CubicSpline()
{
  delete_arrays();
}

void CubicSpline::periodic_coeff(double* pAux, int N, const CSplineInfo& info)
{
  double* dx = pAux;
  double* dv = pAux + N;
  double* alp = dv + N;
  double* bet = alp + N + 1;
  double* gam = bet + N + 1;
// Run forward.
  double A, B, C, F, D;
  alp[0] = 0;
  bet[0] = 0;
  gam[0] = 1;
  for(int i = 1; i <= N; i++)
  {
    A = -dx[i-1];
    B = (i < N) ? -dx[i] : -dx[0];
    C = -2 * (A + B);
    F = (i < N) ? dv[i] / dx[i] - dv[i-1] / dx[i-1] : dv[0] / dx[0] - dv[N-1] / dx[N-1];

    D = C - A * alp[i-1];   // assert: fabs(D) > scConst_Almost_Zero!
    alp[i] = B / D;
    bet[i] = (F + A * bet[i-1]) / D;
    gam[i] = A * gam[i-1] / D;
  }

  double* u = gam + N + 1;
  double* v = u + N + 1;  // u and v are solutions of two auxiliary tasks (Samarski & Nikolayev).
// Run backward.  
  u[N] = 0;
  v[N] = 1;
  for(int i = N - 1; i >= 1; i--)
  {
    u[i] = alp[i] * u[i+1] + bet[i];
    v[i] = alp[i] * v[i+1] + gam[i];
  }

  double* s = v + N + 1;  // the resulting coefficients.
  
  s[0] = (bet[N] + alp[N] * u[1]) / (1.- gam[N] - alp[N] * v[1]);
  s[N] = s[0];
  for(int i = 1; i < N; i++)
    s[i] = u[i] + s[0] * v[i];

  assign(dx, dv, s, N, info);
}

void CubicSpline::coeff(double* pAux, int N, const CSplineInfo& info)
{
  double* dx = pAux;
  double* dv = pAux + N;
  double* alp = dv + N;
  double* bet = alp + N + 1;

  double A, B, C, F, D;
// Run forward:
  switch(info.nBoundType)
  {
    case btFirstDeriv: B = -dx[0]; C = 2 * dx[0]; F = dv[0] / dx[0] - m_fVx0; break;
    case btSecondDeriv: B = 0; C = 1; F = m_fVxx0 / 6.; break;
    case btFree: B = 0; C = 1; F = 0; break; // zero second derivative.
  }

  alp[0] = B / C;
  bet[0] = F / C;
  for(int i = 1; i < N; i++)
  {
    A = -dx[i-1];
    B = -dx[i];
    C = -2 * (A + B);
    F = dv[i] / dx[i] - dv[i-1] / dx[i-1];

    D = C - A * alp[i-1];   // assert: fabs(D) > scConst_Almost_Zero!
    alp[i] = B / D;
    bet[i] = (F + A * bet[i-1]) / D; 
  }

  switch(info.nBoundType)
  {
    case btFirstDeriv: A = -dx[N-1]; B = 0; C = 2 * dx[N-1]; F = m_fVxN - dv[N-1] / dx[N-1]; break;
    case btSecondDeriv: A = 0; B = 0; C = 1; F = m_fVxxN / 6.; break;
    case btFree: A = 0; B = 0; C = 1; F = 0; break; // zero second derivative.
  }
  
// Run backward:
  double* s = bet + N + 1;  // the resulting coefficients.
  s[N] = (F + A * bet[N-1]) / (C - A * alp[N-1]);
  for(int i = N - 1; i >= 0; i--)
    s[i] = alp[i] * s[i+1] + bet[i];

  assign(dx, dv, s, N, info);
}

void CubicSpline::set_node_val(int nNodeIndex, double fValue)
{
  if(m_vSplData.size() == 0)
    allocate_arrays();

  if(nNodeIndex < 0 || nNodeIndex >= m_vSplData.size())
    return;

  m_vSplData[nNodeIndex]->V = fValue;
  m_bReady = false;
}

/*virtual*/void CubicSpline::delete_arrays()
{
  std::vector<CubicSplineNode*>::iterator it;
  for(it = m_vSplData.begin(); it != m_vSplData.end(); it++)
    delete *it;

  m_vSplData.erase(m_vSplData.begin(), m_vSplData.end());
}

int CubicSpline::even_interval(double x, const CPointVector& v) const
{
  int nIntervCount = v.size() - 1;
  if(nIntervCount <= 0)
    return -1;

  double fMinX = v[0];
  double fMaxX = v[nIntervCount];
  if(x < fMinX || x > fMaxX)
    return -1;

  int nEvenInterv = int(nIntervCount * (x - fMinX) / (fMaxX - fMinX));
  if(nEvenInterv == nIntervCount)
    nEvenInterv--;
      
  return nEvenInterv;
}

int CubicSpline::interval(double x, const CPointVector& v, int& nInterv, bool bEven) const
{
  if(bEven)
    return even_interval(x, v);
  
  int nIntervCount = v.size() - 1;
  if(nIntervCount <= 0)
    return -1;

  double fMinX = v[0];
  double fMaxX = v[nIntervCount];
  if(x < fMinX || x > fMaxX)
    return -1;

  if((nInterv >= 0) && (nInterv < nIntervCount))
  {
// First of all, try the same interval.
    if(inside(x, v, nInterv))
      return nInterv;

// Try neighbouring intervals.    
    if(nInterv >= 1)
    {
      nInterv--;
      if(inside(x, v, nInterv))
        return nInterv;
    }

    if(nInterv < nIntervCount - 1)
    {
      nInterv++;
      if(inside(x, v, nInterv))
        return nInterv;
    }
  }

  nInterv = divide_by_two(x, v);
  return nInterv;
}

bool CubicSpline::inside(double x, const CPointVector& v, int nInterv) const
{
  return (x >= v[nInterv]) && (x < v[nInterv + 1]);
}

int CubicSpline::divide_by_two(double x, const CPointVector& v) const
{
  int nLeft = 0;
  int nRight = v.size() - 1;
  while(nRight - nLeft > 1)
  {
    int i = (nLeft + nRight) / 2;
    if(x > v[i])
      nLeft = i;
    else if(x < v[i])
      nRight = i;
    else
      return i < v.size() - 1 ? i : i - 1; // x is equal to v[i] exactly.
  }

  return nLeft;
}

void CubicSpline::set_equidistant(double fMin, double fMax, CPointVector& v)
{
  double fKsi;
  int Ni = v.size();
  for(int i = 0; i < Ni; i++)
  {
    fKsi = double(i) / (Ni - 1);
    v[i] = fMin * (1.- fKsi) + fMax * fKsi;
  }
}

void CubicSpline::set_x(int nNodeX, double fValX)
{
  if(nNodeX < 0 || nNodeX >= m_vX.size())
    return;

  m_vX[nNodeX] = fValX;

  if(nNodeX == 0)
    m_fMinX = fValX;

  if(nNodeX == m_vX.size() - 1)
    m_fMaxX = fValX;
    
  m_bReady = false;
}

void CubicSpline::get_range(double& fMinVal, double& fMaxVal) const
{
  fMinVal = FLT_MAX, fMaxVal = -FLT_MAX;
  for(int i = m_vSplData.size() - 1; i >= 0; i--)
  {
    const CubicSplineNode* pNode = m_vSplData.at(i);

    if(fMinVal > pNode->V)
      fMinVal = pNode->V;

    if(fMaxVal < pNode->V)
      fMaxVal = pNode->V;
  }
}

//------------------------------------------------------------------------------
//  CubSpline1D - a simple cubic spline with different boundary conditions.
//------------------------------------------------------------------------------
CubSpline1D::CubSpline1D(int nBoundType)
  : CubicSpline(),
    m_nBoundType(nBoundType),
    m_bEven(false)
{
  set_default();
}

CubSpline1D::~CubSpline1D()
{
}

double CubSpline1D::get(double x, int nOrder, int& nInterv) const
{
  if(!m_bReady)
    return FLT_MAX;
  
  nInterv = interval(x, m_vX, nInterv, m_bEven);
  if(nInterv < 0)
    return FLT_MAX;
    
  const CubSpline1DNode* pNode = (CubSpline1DNode*)m_vSplData.at(nInterv);
  double dx = x - m_vX[nInterv];
  
  switch(nOrder)
  {
    case 0: return pNode->V + dx * (pNode->b + dx * (pNode->c + dx * pNode->d));
    case 1: return pNode->b + dx * (2 * pNode->c + 3 * dx * pNode->d);
    case 2: return 2 * pNode->c + 6 * dx * pNode->d;
    case 3: return 6 * pNode->d;
  }

  return 0;
}

double CubSpline1D::integral(double x1, double x2) const
{
  int nLast = m_vX.size() - 1;
  if(nLast < 1)
    return 0.;

  double fLeft = m_vX[0];
  if((x1 <= fLeft) && (x2 <= fLeft))
    return 0.;

  double fRight = m_vX[nLast];
  if((x1 >= fRight) && (x2 >= fRight))
    return 0.;
    
  if(!m_bReady)
    return 0.;

  int i1 = 0, i2 = 0;
  i1 = interval(x1, m_vX, i1, m_bEven);
  if(i1 < 0)
    i1 = 0;

  i2 = interval(x2, m_vX, i2, m_bEven);
  if(i2 < 0)
    i2 = nLast - 1;

  static const double c2 = 0.5, c3 = 1./3., c4 = 0.25;

  double dx;
  const CubSpline1DNode* pNode = (CubSpline1DNode*)m_vSplData.at(i1);
  if(i1 == i2)
  {
    dx = x2 - x1;
    return dx * (pNode->V + dx * (c2 * pNode->b + dx * (c3 * pNode->c + dx * c4 * pNode->d)));
  }

// The first interval:
  dx = m_vX[i1 + 1] - x1;
  double S = dx * (pNode->V + dx * (c2 * pNode->b + dx * (c3 * pNode->c + dx * c4 * pNode->d)));

// Intermediate intervals (if any):
  for(int i = i1 + 1; i < i2; i++)
  {
    dx = m_vX[i + 1] - m_vX[i];
    pNode = (CubSpline1DNode*)m_vSplData.at(i);
    S += dx * (pNode->V + dx * (c2 * pNode->b + dx * (c3 * pNode->c + dx * c4 * pNode->d)));
  }

// The last interval:
  dx = x2 - m_vX[i2];
  pNode = (CubSpline1DNode*)m_vSplData.at(i2);
  S += dx * (pNode->V + dx * (c2 * pNode->b + dx * (c3 * pNode->c + dx * c4 * pNode->d)));

  return S;
}

void CubSpline1D::set_nodes_count(int nCount)
{
  if((nCount != m_vX.size()) && (nCount > 1))
  {
    m_nNodesCount = nCount;
    m_vX.resize(m_nNodesCount);
    allocate_arrays();
    m_bReady = false;
  }
}

void CubSpline1D::set_bound_derivatives(double fV0, double fVN)
{
  if(m_nBoundType == btFirstDeriv)
  {
    if(m_fVx0 != fV0 || m_fVxN != fVN)
    {
      m_fVx0 = fV0;
      m_fVxN = fVN;
      m_bReady = false;
    }
  }
  else if(m_nBoundType == btSecondDeriv)
  {
    if(m_fVxx0 != fV0 || m_fVxxN != fVN)
    {
      m_fVxx0 = fV0;
      m_fVxxN = fVN;
      m_bReady = false;
    }
  }
}

// Returns false and does nothing if this node cannot be inserted. For example,
// if (fX, fValue) is too close to an existing node. Attention! The existing
// array of x is supposed to be ordered!
bool CubSpline1D::insert_node(double fX, double fValue)
{
  if(m_bEven)
    return false; // the user must call set_node_val(int, double) if m_bEven is true.
    
  int i;
  for(i = m_vX.size() - 1; i >= 0; i--)
  {
    if(fX > m_vX.at(i))
      break;
  }

  if(m_vX.size() > 0)
  {
    if((i > 0) && (fX - m_vX.at(i) < Const_Almost_Zero))
      return false; // the new node is too close to its left neighbour.

    if((i != m_vX.size() - 1) && (m_vX.at(i + 1) - fX < Const_Almost_Zero))
      return false; // the new node is too close to its right neighbour.
  }

  CubSpline1DNode* pNode = new CubSpline1DNode();
  pNode->V = fValue;
  
  m_vX.insert(m_vX.begin() + i + 1, fX);
  m_vSplData.insert(m_vSplData.begin() + i + 1, (CubicSplineNode*)pNode);

  if(m_vX.size() != m_vSplData.size())
  {
    AfxMessageBox("CubSpline1D: Sizes of X and Y arrays differ!");
    return false;
  }
  
  return true;
}

// Returns false and does nothing if fX <= the last existing x. Normally, this 
// function adds the (fX, fValue) node to the end of the spline nodes collection.
bool CubSpline1D::push_back_node(double fX, double fValue)
{
  if(m_bEven)
    return false; // the user must call set_node_val(int, double) if m_bEven is true.

  if(m_vX.size() > 0)
  {
    double fLastX = m_vX.at(m_vX.size() - 1);
    if(fX < fLastX || fX - fLastX < Const_Almost_Zero)
      return false;
  }

  CubSpline1DNode* pNode = new CubSpline1DNode();
  pNode->V = fValue;

  m_vX.push_back(fX);
  m_vSplData.push_back(pNode);

  if(m_vX.size() != m_vSplData.size())
  {
    AfxMessageBox("CubSpline1D: Sizes of X and Y arrays differ!");
    return false;
  }
  
  return true;
}

/*virtual*/bool CubSpline1D::check()
{
  int N = m_vX.size();
  for(int i = 1; i < N; i++)  // the X-array must be an increasing series.
  {
    if(m_vX[i - 1] >= m_vX[i])
      return false;
  }

  m_fMinX = m_vX[0];
  m_fMaxX = m_vX[N-1];
  m_fRangeX = m_vX[N-1] - m_vX[0];
  return true;
}

/*virtual*/void CubSpline1D::prepare()
{  
  if(m_bEven)
    set_equidistant(m_fMinX, m_fMaxX, m_vX);

  if(!check())
  {
    AfxMessageBox("CubSpline1D: Bad X-array!");
    return;
  }

  CSplineInfo info;
  info.nBoundType = m_nBoundType;
  
  int N = m_vX.size() - 1;  // number of intervals.
  int nAuxDim = info.nBoundType == btPeriodic ? 8 * N + 6 : 5 * N + 3;
  double* pAux = new double[nAuxDim]; // 2*N + 6*(N+1)  or  2*N + 3*(N+1)
  
  double* dx = pAux;
  double* dv = pAux + N;
  for(int i = 0; i < N; i++)
  {
    dx[i] = m_vX[i+1] - m_vX[i];
    dv[i] = m_vSplData[i+1]->V - m_vSplData[i]->V;  // this spline is one-dimensional.
  }

  if(info.nBoundType == btPeriodic)
    periodic_coeff(pAux, N, info);
  else
    coeff(pAux, N, info);

  delete[] pAux;
  m_bReady = true;
}

/*virtual*/void CubSpline1D::allocate_arrays()
{
  delete_arrays();

  int Ni = m_vX.size();
  m_vSplData.resize(Ni, 0);
  for(int i = 0; i < Ni; i++)
  {
    m_vSplData[i] = new CubSpline1DNode();
  }
}

/*virtual*/void CubSpline1D::assign(double* dx, double* dv, double* s, int N, const CSplineInfo& /*info*/)
{
  for(int i = 0; i < N; i++)
  {
    CubSpline1DNode* pNode = (CubSpline1DNode*)m_vSplData.at(i);
    pNode->b = dv[i] / dx[i] - dx[i] * (s[i+1] + 2 * s[i]);
    pNode->c = 3 * s[i];
    pNode->d = (s[i+1] - s[i]) / dx[i];
  }
}

/*virtual*/int CubSpline1D::index(int i, int /*j*/, int /*k*/)
{
  return i;
}

/*virtual*/void CubSpline1D::set_default()
{
  m_fVx0 = 0;
  m_fVxN = 0;
  m_fVxx0 = 0;
  m_fVxxN = 0;
}

void CubSpline1D::set_range(double fMinX, double fMaxX, bool bEven)
{
  if(m_fMinX != fMinX || m_fMaxX != fMaxX || m_bEven != bEven)
  {
    m_fMinX = fMinX;
    m_fMaxX = fMaxX;
    m_bEven = bEven;
    m_bReady = false;
  }

  m_fRangeX = m_fMaxX - m_fMinX;
  if(m_bEven && (m_vX.size() > 1))
    set_equidistant(m_fMinX, m_fMaxX, m_vX);
}

};  // namespace EvaporatingParticle