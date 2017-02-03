
#pragma once

#include "constant.hpp"
#include "math.h"

class CArchive;

namespace EvaporatingParticle
{

struct CNode3D;
//-------------------------------------------------------------------------------------------------
// CEvaporationModel - the base class of all evaporation models available in this code.
//-------------------------------------------------------------------------------------------------
class CEvaporationModel
{
public:
  CEvaporationModel();

  enum
  {
    mtmNone         = 0,
    mtmRanzMarshall = 1
  };

// Evaporation rate and cooling rate depend on the droplet movement velocity through the environment, i.e. Re.
  virtual double          get_evaporation_rate(const CNode3D& cNode, double fDropTemp, double fDropDiam, double fRe);
  virtual double          get_cooling_rate(const CNode3D& cNode, double fDropTemp, double fDropDiam, double fRe);

// User interface:
  double                  get_env_humidity() const;
  DWORD_PTR               get_env_humidity_ptr() const;
  void                    set_env_humidity(double fHumid);

  bool                    get_enable_surf_tens() const;
  DWORD_PTR               get_enable_surf_tens_ptr() const;
  void                    set_enable_surf_tens(bool bEnable);

  int                     get_mass_trans_model() const;
  DWORD_PTR               get_mass_trans_model_ptr() const;
  void                    set_mass_trans_model(int nMod);

// Moving through environmental gas support:
  double                  get_Pr(const CNode3D& cNode) const;   // Prandtl number.
  double                  get_Sc(const CNode3D& cNode) const;   // Schmidt number.

  double                  get_Nu(double fRe, double fPr) const; // Nusselt number is a function of both Re and Pr.
  double                  get_Sh(double fRe, double fSc) const; // Sherwood number depends on Re and Sc.

// Streams support:
  virtual void            save(CArchive& ar);
  virtual void            load(CArchive& ar);

protected:
  void                    set_default();

// In the 3 functions below fPress is full pressure, Dyn/cm2; fTemp is absolute temperature, K.
  double                  get_humidity(double fH2O_mf, double fPress, double fTemp) const;
  double                  get_H2O_mf(double fHumidity, double fPress, double fTemp) const;

  double                  get_partial_H2O_press(double fHumidity, double fPress, double fTemp) const;

// Returns (1 + dP/P), where dP is increase of saturated pressure due to curvature of the droplet's surface.
  double                  get_sat_press_correction(double fDropTemp, double fDropDiam) const;

//-------------------------------------------------------------------------------------------------
// Empirical data fits (CGS):
//-------------------------------------------------------------------------------------------------
public:
// Regression curve fit to the empirical data from Bolz and Tuve (1976), H2O vapour diffusion
// coefficient in air, cm2/s.
  static double           get_diffusion_coeff(double fT);

// Regresion curve fit to air heat conductivity, erg/(cm K).
  static double           get_air_heat_cond(double fTemp);

// Cubical regression fit to the experimental data. Input: partial pressure of water vapour, Dyn/cm2.
  static bool             get_saturation_temp(double fPressH2O, double& fSatTemp);
// Cubical regression fit to the experimental data. Input: full pressure, Dyn/cm2.
  static bool             get_saturation_temp(double fPress, double fMassFractH2O, double& fSatTemp);

// Cubical regression fit to the experimental data.
// Input: fTemp - absolute T, K. Output: pressure of saturated water vapour, Dyn/cm2.
  static double           get_saturation_press(double fTemp);

// Cubical regression fit for the specific latent H2O vaporization heat. Input: fPress - full pressure, Dyn/cm2.
// Output - latent specific vaporization heat, Erg/g.
  static double           get_specific_latent_heat(double fPress);
// Cubical regression fit for the specific latent H2O vaporization heat. Input: fTemp - absolute temperature, K.
// Output - latent vaporization heat per one H2O molecule, Erg.
  static double           get_latent_heat_per_molecule(double fTemp);

// Cubical regression fit for the surface tension. 273 < fTemp < 647.1 K, where 647.1 K is the critical point T.
  static double           get_surface_tension(double fTemp);

//-------------------------------------------------------------------------------------------------
// User data (CGS):
//-------------------------------------------------------------------------------------------------
protected:
  bool                    m_bEnableSurfTens;

  double                  m_fEnvHumidity,
                          m_H2O_mf;

  int                     m_nMassTransModel;

  friend class COutputEngine;
  friend class CTracker;
};

//-------------------------------------------------------------------------------------------------
// CMaxwellModel - the simplest of all evaporation models, isothermal.
//-------------------------------------------------------------------------------------------------
class CMaxwellModel : public CEvaporationModel
{
public:
  CMaxwellModel();

  virtual double          get_evaporation_rate(const CNode3D& cNode, double fDropTemp, double fDropDiam, double fRe);
};

//-------------------------------------------------------------------------------------------------
// CSteadyDiffusiveModel - the H2O vapour at the surface (r = R) is saturated, the temperature T(R) is
// found from Q = W*L, where Q = 4*pi*lambda*R*(T(inf) - T(R)), W = 4*pi*D*R*(n(R) - n(inf)), L = L(T)
// is the latent vaporization heat per one H2O molecule.
//-------------------------------------------------------------------------------------------------
class CSteadyDiffusiveModel : public CEvaporationModel
{
public:
  CSteadyDiffusiveModel();

  virtual double          get_evaporation_rate(const CNode3D& cNode, double fDropTemp, double fDropDiam, double fRe);

protected:
  static double           heat_balance_equation(double fT, double* pPar);
};

//-------------------------------------------------------------------------------------------------
// CDiffusiveModel - the H2O vapour at the surface (r = R) is saturated, the temperature T(R) is
// found from the heat balance equation at every time step: M*Cv*dT/dt = Q - W*L.
//-------------------------------------------------------------------------------------------------
class CDiffusiveModel : public CEvaporationModel
{
public:
  CDiffusiveModel();

  virtual double          get_evaporation_rate(const CNode3D& cNode, double fDropTemp, double fDropDiam, double fRe);
  virtual double          get_cooling_rate(const CNode3D& cNode, double fDropTemp, double fDropDiam, double fRe);
};

//-------------------------------------------------------------------------------------------------
// CEvaporationModel inlines:
//-------------------------------------------------------------------------------------------------
inline bool CEvaporationModel::get_enable_surf_tens() const
{
  return m_bEnableSurfTens;
}

inline DWORD_PTR CEvaporationModel::get_enable_surf_tens_ptr() const
{
  return (DWORD_PTR)&m_bEnableSurfTens;
}

inline void CEvaporationModel::set_enable_surf_tens(bool bEnable)
{
  m_bEnableSurfTens = bEnable;
}

inline double CEvaporationModel::get_env_humidity() const
{
  return m_fEnvHumidity;
}

inline DWORD_PTR CEvaporationModel::get_env_humidity_ptr() const
{
  return (DWORD_PTR)&m_fEnvHumidity;
}

inline void CEvaporationModel::set_env_humidity(double fHumid)
{
  m_fEnvHumidity = fHumid;
  m_H2O_mf = get_H2O_mf(fHumid, Const_One_Atm_CGS, 300.0);
}

inline int CEvaporationModel::get_mass_trans_model() const
{
  return m_nMassTransModel;
}

inline DWORD_PTR CEvaporationModel::get_mass_trans_model_ptr() const
{
  return (DWORD_PTR)&m_nMassTransModel;
}

inline void CEvaporationModel::set_mass_trans_model(int nMod)
{
  m_nMassTransModel = nMod;
}

// Regression fit to the empirical data of H2O vapor diffusion in Air from Bolz and Tuve (1976), cm2/s.
inline double CEvaporationModel::get_diffusion_coeff(double fT)
{
  return 2.775e-2 + (4.479e-4 + 1.656e-6 * fT) * fT;
}

// Regression fit, erg/(cm s K).
inline double CEvaporationModel::get_air_heat_cond(double fT)
{
  return fT * (9.693 + fT * (-4.294e-3 + fT * 1.777e-6));
}

// Regression fit. Input: fPress - full pressure, Dyn/cm2. Output - latent specific vaporization heat, Erg/g.
inline double CEvaporationModel::get_specific_latent_heat(double fPress)
{
  double x = log10(fPress / Const_One_Atm_CGS);
  return Const_Cal_to_Erg * (539.551 + x * (-42.132 + x * (-12.805 - 2.6936 * x)));
}

inline double CEvaporationModel::get_latent_heat_per_molecule(double fT)
{
  static const double scfCoeff = Const_Cal_to_Erg / Const_H2O_per_Gram;
  return scfCoeff * (832.4 + fT * (-1.41 + fT * (3.0e-3 - fT * 3.41e-6)));
}

inline double CEvaporationModel::get_surface_tension(double fT)
{
  static const double scfCoeff = 1000.; // 1 N/m = 1000 Dyn/cm.
  if(fT < 273 || fT > 644.5)
    return 0.;  // real critical temperature is 647.1, but the fit intersects the x-axis at T ~ 644.7.

  return scfCoeff * (0.0694 + fT * (2.214e-4 + fT * (-8.888e-7 + fT * 5.87e-10)));
}

inline double CEvaporationModel::get_sat_press_correction(double fDropTemp, double fDropDiam) const
{
  if(fDropTemp < 273)
    return 1.;

  static const double scfVolOneMolec = Const_H2O_CGS / (Const_Boltzmann * Const_H2O_Dens);

  return 1. + 4 * scfVolOneMolec * get_surface_tension(fDropTemp) / (fDropTemp * fDropDiam);
}

};  // namespace EvaporatingParticle