
#ifndef _Constant_
#define _Constant_

const double
  Const_Almost_Zero       = 1.e-14,
  Const_PI                = 3.14159265358979323846,
  Const_2PI               = 6.28318530717958647692,
  Const_1_PI              = 0.318309886183790671538,
  Const_Half_PI           = 1.57079632679489661923,
  Const_E                 = 2.71828182845904523536,
  Const_Sqrt2             = 1.41421356237309504880,
  Const_One_Third         = 0.33333333333333333333,
  Const_One_Sixth         = 0.16666666666666666667,
  Const_RadianToDegree    = 180. / Const_PI,
  Const_DegreeToRadian    = Const_PI / 180.,
  Const_Ln10              = 2.302585092994045684018;

const double
  Const_One_Atm_CGS       = 1.01325e+6,   // Dyn/cm2
  Const_One_Atm_SI        = 1.01325e+5,   // Pa
  Const_Boltzmann         = 1.3806504e-16,// Erg/K
  Const_Boltzmann_SI      = 1.3806504e-23,// J/K
  Const_Avogadro          = 6.022142e+23, // 1/mol
  Const_Molar_Mass_Air    = 28.98,        // g/mol
  Const_Molar_Mass_H2O    = 18.0,
  Const_Mol_Mass_Ratio    = 0.62112,
  Const_Min_H2O_mf        = 1.0e-5 * Const_Mol_Mass_Ratio,
  Const_H2O_per_Gram      = Const_Avogadro / Const_Molar_Mass_H2O,  // 1/g
  Const_H2O_Heat_Capacity = 4.184e+7,     // erg/(g K)
  Const_H2O_Dens          = 1.0,          // g/cm3
  Const_Air_Dens          = 0.001204,     // g/cm3 at P = 1 atm and T = 293 K, room conditions.
  Const_T0                = 273.16;		    // 0 C in K added by AC 23/06/2016
 
const double
  SI_to_CGS_Len         = 100.0,
  SI_to_CGS_Weight      = 1000.0,
  SI_to_CGS_Vel         = 100.0,
  SI_to_CGS_Press       = 10.0,       // 1 Pa = 10 Dyn/cm2
  SI_to_CGS_Dens        = 0.001,      // 1 kg/m3 = 0.001 g/cm3
  SI_to_CGS_NumDens     = 1.e-6,      // 1 m-3 = 1.e-6 cm-3
  SI_to_CGS_DynVisc     = 10.0,       // 1 Pa s = 10 CGS units
  SI_to_CGS_ThermCond   = 1.e+5,      // 1 J/(s m K) = 1e+5 Erg/(s cm K)
  SI_to_CGS_Voltage     = 1./299.8,   // 1 V in CGSE
  SI_to_CGS_ElecField   = 0.01 * SI_to_CGS_Voltage,
  SI_to_CGS_Cp          = 1.0e+4,     // 1 J/(kg K) = 1e+4 Erg/(g K)
  SI_to_CGS_Energy      = 1.0e+7;     // 1 J = 1e+7 Erg

const double
  CGS_to_SI_Len         = 1. / SI_to_CGS_Len,
  CGS_to_SI_Weight      = 1. / SI_to_CGS_Weight,
  CGS_to_SI_Vel         = 1. / SI_to_CGS_Vel,
  CGS_to_SI_Press       = 1. / SI_to_CGS_Press,
  CGS_to_SI_Dens        = 1. / SI_to_CGS_Dens,
  CGS_to_SI_NumDens     = 1. / SI_to_CGS_NumDens,
  CGS_to_SI_DynVisc     = 1. / SI_to_CGS_DynVisc,
  CGS_to_SI_ThermCond   = 1. / SI_to_CGS_ThermCond,
  CGS_to_SI_Voltage     = 1. / SI_to_CGS_Voltage,
  CGS_to_SI_ElecField   = 1. / SI_to_CGS_ElecField,
  CGS_to_SI_Energy      = 1. / SI_to_CGS_Energy;

const double
  Const_Cal_to_Erg      = 4.184e+7; // 1 Cal = 4.184 J, 1 J = 1e+7 Erg.

const double
  Const_AMU_SI          = 1.66042e-27,	// atomic mass unit, kg
  Const_AMU_CGS         = 1.66042e-24,	// atomic mass unit, g
  Const_H2O_CGS         = 2.988756e-23, // g.
  Const_Angstrem_CGS    = 1.0e-8;       // 1 angstrem, cm.

const double 
  Const_Charge_CGS      = 4.803e-10,    // electron charge, CGSE units
  Const_Erg_to_EV				= 6.242e+11;		// 1 erg in eV

const double
  Const_C_CGS           = 2.99792458e+10, // speed of light, cm/s.
  Const_nA_to_CGSE      = 2.99792458;     // current from nA to CGSE.

const COLORREF clMaroon = 128;
const COLORREF clGreen  = 32768;
const COLORREF clOlive  = 32896;
const COLORREF clNavy   = 8388608;
const COLORREF clPurple = 8388736;
const COLORREF clTeal   = 8421376;
const COLORREF clGray   = 8421504;
const COLORREF clSilver = 12632256;
const COLORREF clRed    = 255;
const COLORREF clLime   = 65280;
const COLORREF clYellow = 65535;
const COLORREF clBlue   = 16711680;
const COLORREF clFuchsia = 16711935;
const COLORREF clAqua   = 16776960;
const COLORREF clLtGray = 12632256;
const COLORREF clDkGray = 8421504;
const COLORREF clWhite  = 16777215;
const COLORREF clSkyBlue = 15780518;
const COLORREF clCream  = 15793151;
const COLORREF clMedGray = 10789024;

#endif
