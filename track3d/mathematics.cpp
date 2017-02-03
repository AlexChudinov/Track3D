
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

};  // namespace EvaporatingParticle