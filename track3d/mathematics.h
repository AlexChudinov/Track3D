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
};

};  // namespace EvaporatingParticle

#endif // _Mathematics_