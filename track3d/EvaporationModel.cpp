
#include "stdafx.h"

#include "EvaporationModel.h"
#include "mathematics.h"
#include "Tracker.hpp"

namespace EvaporatingParticle
{
//-------------------------------------------------------------------------------------------------
// CEvaporationModel - the base class of all evaporation models available in this code.
//-------------------------------------------------------------------------------------------------
CEvaporationModel::CEvaporationModel()
{
  set_default();
}

void CEvaporationModel::set_default()
{
  set_env_humidity(0.5);
  set_enable_surf_tens(true);
  set_mass_trans_model(mtmRanzMarshall);
}

double CEvaporationModel::get_evaporation_rate(const CNode3D& node, double fDropTemp, double fDropDiam, double fRe)
{
  return 0.;  // for non-zero rate see the descendants.
}

double CEvaporationModel::get_cooling_rate(const CNode3D& cNode, double fDropTemp, double fDropDiam, double fRe)
{
  return 0.;
}

double CEvaporationModel::get_Pr(const CNode3D& cNode) const   // Prandtl number.
{
  return cNode.cond > Const_Almost_Zero ? cNode.visc * cNode.cp / cNode.cond : 0;
}

double CEvaporationModel::get_Sc(const CNode3D& cNode) const   // Schmidt number.
{
  return cNode.visc / (cNode.dens * get_diffusion_coeff(cNode.temp));
}

double CEvaporationModel::get_Nu(double fRe, double fPr) const // Nusselt number is a function of both Re and Pr.
{
  if(fRe < Const_Almost_Zero || fPr < Const_Almost_Zero)
    return 2;

  switch(m_nMassTransModel)
  {
    case mtmNone:
      return 2;
    case mtmRanzMarshall: 
      return 2 + exp(-0.5108256 + 0.5 * log(fRe) + 0.333 * log(fPr));
  }

  return 2;
}

double CEvaporationModel::get_Sh(double fRe, double fSc) const // Sherwood number depends on Re and Sc.
{
  if(fRe < Const_Almost_Zero || fSc < Const_Almost_Zero)
    return 2;

  switch(m_nMassTransModel)
  {
    case mtmNone:
      return 2;
    case mtmRanzMarshall: 
      return 2 + exp(-0.5108256 + 0.5 * log(fRe) + 0.333 * log(fSc));
  }

  return 2;
}

// Regression fit of the experimental data, K.
bool CEvaporationModel::get_saturation_temp(double fPressH2O, double& fSatTemp)
{
  double r = fPressH2O / Const_One_Atm_CGS;
  if(r < 1.e-5)
    return false; // saturation temperature is not defined if partial pressure of H2O tends to zero.

  double x = log10(r);
  fSatTemp = 373.0 + x * (65.48426 + x * (12.09 + x * 1.26626));
  return true;
}

bool CEvaporationModel::get_saturation_temp(double fPress, double fMassFractH2O, double& fSatTemp)
{
  if(fMassFractH2O < Const_Min_H2O_mf)
    return false; // saturation temperature is not defined if partial pressure of H2O tends to zero.

  double fPressH2O = fPress / (1. + Const_Mol_Mass_Ratio * (1. - fMassFractH2O) / fMassFractH2O);
  return get_saturation_temp(fPressH2O, fSatTemp);
}

// Cubical regression fit to the experimental data.
// Input: fTemp - absolute T, K. Output: pressure of saturated water vapour, Dyn/cm2.
double CEvaporationModel::get_saturation_press(double fTemp)
{
  double fLg10Ratio = -23.59 + fTemp * (0.14065 + fTemp * (-2.861e-4 + fTemp * 2.1085e-7));
  return Const_One_Atm_CGS * exp(Const_Ln10 * fLg10Ratio);
}

double CEvaporationModel::get_H2O_mf(double fHumidity, double fPress, double fTemp) const
{
  if(fHumidity < Const_Almost_Zero)
    return 0.;

  double fRatio = -1. + fPress / (fHumidity * get_saturation_press(fTemp));
  return 1. / (1. + fRatio / Const_Mol_Mass_Ratio);
}

double CEvaporationModel::get_humidity(double fH2O_mf, double fPress, double fTemp) const
{
  if(fH2O_mf < Const_Almost_Zero)
    return 0.;

  double fR1 = fPress / get_saturation_press(fTemp);
  double fR2 = Const_Mol_Mass_Ratio * (1. - fH2O_mf) / fH2O_mf;
  return fR1 / (1. + fR2);
}

double CEvaporationModel::get_partial_H2O_press(double fHumidity, double fPress, double fTemp) const
{
  if(fHumidity < Const_Almost_Zero)
    return 0.;

  double fH2O_mf = get_H2O_mf(fHumidity, fPress, fTemp);
  double fRatio = Const_Mol_Mass_Ratio * (1. - fH2O_mf) / fH2O_mf;

  return fPress / (1. + fRatio);
}

// Streams support:
void CEvaporationModel::save(CArchive& ar)
{
  const UINT nVersion = 1;
  ar << nVersion;

  ar << m_bEnableSurfTens;
  ar << m_fEnvHumidity;
  ar << m_nMassTransModel;
}

void CEvaporationModel::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_bEnableSurfTens;
  ar >> m_fEnvHumidity;

  m_H2O_mf = get_H2O_mf(m_fEnvHumidity, Const_One_Atm_CGS, 300.0);

  if(nVersion >= 1)
    ar >> m_nMassTransModel;
}

//-------------------------------------------------------------------------------------------------
// CMaxwellModel - the simplest of all evaporation models, isothermal.
//-------------------------------------------------------------------------------------------------
CMaxwellModel::CMaxwellModel()
  : CEvaporationModel()
{
}

double CMaxwellModel::get_evaporation_rate(const CNode3D& cNode, double fDropTemp, double fDropDiam, double fRe)
{
  double fDiffCoeff = get_diffusion_coeff(cNode.temp);

  double fEnvPressH2O = get_partial_H2O_press(m_fEnvHumidity, cNode.press, cNode.temp);
  double fEnvNumDens = fEnvPressH2O / (Const_Boltzmann * cNode.temp);

// In the Maxwell theory the temperature of the droplet is always equal to the temperature of the environment
// and the H2O vapours near the droplet boundary are always saturated:
  double fSatPressH2O = get_saturation_press(cNode.temp);
  if(m_bEnableSurfTens)
    fSatPressH2O *= get_sat_press_correction(fDropTemp, fDropDiam);

  double fSatNumDens = fSatPressH2O / (Const_Boltzmann * cNode.temp);

  double fSh = m_nMassTransModel == mtmNone ? 2 : get_Sh(fRe, get_Sc(cNode));

  return Const_PI * fDropDiam * fSh * fDiffCoeff * (fSatNumDens - fEnvNumDens) * Const_H2O_CGS;
}

//-------------------------------------------------------------------------------------------------
// CSteadyDiffusiveModel - the H2O vapour at the surface (r = R) is saturated, the temperature T(R) is
// found from Q = W*L, where Q = 4*pi*lambda*R*(T(inf) - T(R)), W = 4*pi*D*R*(n(R) - n(inf)), L = L(T)
// is the latent vaporization heat per one H2O molecule.
//-------------------------------------------------------------------------------------------------
CSteadyDiffusiveModel::CSteadyDiffusiveModel()
  : CEvaporationModel()
{
}

double CSteadyDiffusiveModel::get_evaporation_rate(const CNode3D& cNode, double fDropTemp, double fDropDiam, double fRe)
{
  double fEnvPressH2O = get_partial_H2O_press(m_fEnvHumidity, cNode.press, cNode.temp);
  double fEnvNumDens = fEnvPressH2O / (Const_Boltzmann * cNode.temp);

  double fSh = m_nMassTransModel == mtmNone ? 2 : get_Sh(fRe, get_Sc(cNode));
  double fNu = m_nMassTransModel == mtmNone ? 2 : get_Nu(fRe, get_Pr(cNode));

  double pParam[4] = { fEnvNumDens, cNode.temp, fSh, fNu };
// In the steady-state diffusion theory the temperature of the droplet is obtained from the heat balance
// equation Q = W*L so that the input parameter fDropTemp is ignored.
  double fT = CMath::root(10., 550., &heat_balance_equation, 0.1, pParam);

// H2O vapours near the droplet boundary are saturated at the temperature obtained above:
  double fSatPressH2O = get_saturation_press(fT);
  if(m_bEnableSurfTens)
    fSatPressH2O *= get_sat_press_correction(fDropTemp, fDropDiam);

  double fSatNumDens = fSatPressH2O / (Const_Boltzmann * fT);
  double fDiffCoeff = get_diffusion_coeff(fT);

  return Const_PI * fDropDiam * fSh * fDiffCoeff * (fSatNumDens - fEnvNumDens) * Const_H2O_CGS;
}

// pPar must contain: pPar[0] = n(inf), 1/cm3; pPar[1] = T(inf), K, pPar[2] = Sherwood, pPar[3] = Nusselt.
double CSteadyDiffusiveModel::heat_balance_equation(double fT, double* pPar)
{
  double fLatentHeat = get_latent_heat_per_molecule(fT);
  double fDiffCoeff = get_diffusion_coeff(fT);

  double fSatPressH2O = get_saturation_press(fT);
  double fSatNumDens = fSatPressH2O / (Const_Boltzmann * fT);

  double fHeatCond = get_air_heat_cond(fT);

  double fEnvNumDens = pPar[0];
  double fEnvTemp = pPar[1];
  double fSh = pPar[2];
  double fNu = pPar[3];

  return fLatentHeat * fDiffCoeff * fSh * (fSatNumDens - fEnvNumDens) - fHeatCond * fNu * (fEnvTemp - fT);
}

//-------------------------------------------------------------------------------------------------
// CDiffusiveModel - the H2O vapour at the surface (r = R) is saturated, the temperature T(R) is
// found from the heat balance equation at every time step: M*Cv*dT/dt = Q - W*L.
//-------------------------------------------------------------------------------------------------
CDiffusiveModel::CDiffusiveModel()
  : CEvaporationModel()
{
}

double CDiffusiveModel::get_evaporation_rate(const CNode3D& cNode, double fDropTemp, double fDropDiam, double fRe)
{
  double fEnvPressH2O = get_partial_H2O_press(m_fEnvHumidity, cNode.press, cNode.temp);
  double fEnvNumDens = fEnvPressH2O / (Const_Boltzmann * cNode.temp);

  double fDiffCoeff = get_diffusion_coeff(fDropTemp);

// H2O vapour near the droplet boundary is saturated at the temperature fDropTemp:
  double fSatPressH2O = get_saturation_press(fDropTemp);
  if(m_bEnableSurfTens)
    fSatPressH2O *= get_sat_press_correction(fDropTemp, fDropDiam);

  double fSatNumDens = fSatPressH2O / (Const_Boltzmann * fDropTemp);

  double fSh = m_nMassTransModel == mtmNone ? 2 : get_Sh(fRe, get_Sc(cNode));

  return Const_PI * fDropDiam * fSh * fDiffCoeff * (fSatNumDens - fEnvNumDens) * Const_H2O_CGS;
}

double CDiffusiveModel::get_cooling_rate(const CNode3D& cNode, double fDropTemp, double fDropDiam, double fRe)
{
  static const double scfHeatCapacityH2O = Const_H2O_Dens * Const_H2O_Heat_Capacity;

  double fEnvPressH2O = get_partial_H2O_press(m_fEnvHumidity, cNode.press, cNode.temp);
  double fEnvNumDens = fEnvPressH2O / (Const_Boltzmann * cNode.temp);

  double fDiffCoeff = get_diffusion_coeff(fDropTemp);
  double fLatentHeat = get_latent_heat_per_molecule(fDropTemp);
  double fAirHeatCond = get_air_heat_cond(fDropTemp);

// H2O vapour near the droplet boundary is saturated at the temperature fDropTemp:
  double fSatPressH2O = get_saturation_press(fDropTemp);
  double fSatNumDens = fSatPressH2O / (Const_Boltzmann * fDropTemp);

  double fdT = cNode.temp - fDropTemp;
  double fdN = fSatNumDens - fEnvNumDens;

  double fSh = m_nMassTransModel == mtmNone ? 2 : get_Sh(fRe, get_Sc(cNode));
  double fNu = m_nMassTransModel == mtmNone ? 2 : get_Nu(fRe, get_Pr(cNode));

  return -6 * (fNu * fAirHeatCond * fdT - fLatentHeat * fSh * fDiffCoeff * fdN) / (fDropDiam * fDropDiam * scfHeatCapacityH2O);
}

};  // namespace EvaporatingParticle