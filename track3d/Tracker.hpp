#pragma once
#ifndef _TRACKER_
#define _TRACKER_

#include "TrackItem.h"
#include "OutputEngine.h"
#include "ParticleSource.h"
#include "ImportOpenFOAM.h"
#include "BeamCrossSection.h"
#include "../RandomProcess/RandomProcess.h"    // random diffusion support.
#include "Perturbation.h"
#include "matrix3d.hpp"
#include "math.h"

class CExecutionDialog;
class DiffusionVelocityJump;

typedef RandomProcess::RandomProcessType CRandomDiffType;

namespace EvaporatingParticle
{

const double cfDCVoltage = 2500.; // ANSYS data supplies DC field computed with 2500 V potential at the emitter.

class CSource;
class CEvaporationModel;
class CalcThreadVector;
class CBarnesHut;
//-------------------------------------------------------------------------------------------------
// CTracker - main class containing data and methods for ions and evaporating particles tracking.
//-------------------------------------------------------------------------------------------------
class CTracker : public CObject
{
public:
  CTracker();
  ~CTracker();

  enum  // Evaporation model type
  {
    emNone            = 0,
    emMaxwell         = 1,
    emSteadyDiffusive = 2,
    emDiffusive       = 3
  };

  enum  // Symmetry planes
  {
    spNone  = 0,
    spXY    = 1,
    spXZ    = 2,
    spYZ    = 4
  };

	enum // Integrator type
	{
    intExplEuler = 0,
    intModMidpnt = 1,
		intRK2       = 2,
		intRK4       = 3,
		intCount     = 4
	};

  static void             GetTimeDeriv(void* pData, const double* pItemState, double* pTimeDeriv, const double* pTime);

  int                     get_integr_type() const;
  DWORD_PTR               get_integr_type_ptr() const;
  void                    set_integr_type(int nType);

  const char*             get_integr_name(int nType) const;

  const char*             get_filename() const;
  DWORD_PTR               get_filename_ptr() const;
  bool                    set_filename(const char* pName);

  bool                    get_convert_to_cgs() const;
  void                    set_convert_to_cgs(bool bEnable);

  int                     get_particle_type() const;
  DWORD_PTR               get_particle_type_ptr() const;
  void                    set_particle_type(int nType);

  double                  get_molar_mass() const;
  DWORD_PTR               get_molar_mass_ptr() const;
  void                    set_molar_mass(double fMass);   // g / mol

// Particle tracking:
  double                  get_time_step() const;
  DWORD_PTR               get_time_step_ptr() const;
  void                    set_time_step(double fTimeStep);

  double                  get_max_track_time() const;
  DWORD_PTR               get_max_track_time_ptr() const;
  void                    set_max_track_time(double fTime);

  bool                    get_enable_field() const;
  DWORD_PTR               get_enable_field_ptr() const;
  void                    set_enable_field(bool bEnable);

  double                  get_particle_charge() const;
  DWORD_PTR               get_particle_charge_ptr() const;
  void                    set_particle_charge(double fCharge);

  double                  get_dc_amplitude() const;
  DWORD_PTR               get_dc_amplitude_ptr() const;
  void                    set_dc_amplitude(double fAmpl);

  int                     get_sym_plane() const;
  DWORD_PTR               get_sym_plane_ptr() const;
  void                    set_sym_plane(int nPlane);

  int                     get_symmetry_type() const;

// Particle source parameters:
  CSource*                get_src();

// Evaporation:
  int                     get_evapor_model_type() const;
  DWORD_PTR               get_evapor_model_type_ptr() const;
  void                    set_evapor_model_type(int nModel);

  double                  get_init_diameter() const;
  DWORD_PTR               get_init_diameter_ptr() const;
  void                    set_init_diameter(double fDiam);

// All evaporation user interface must be called directly using this pointer.
  CEvaporationModel*      get_evapor_model() const;

// Ion type of particles:
  double                  get_ion_mass() const;
  DWORD_PTR               get_ion_mass_ptr() const;
  void                    set_ion_mass(double fMass);

  double                  get_ion_mobility() const;     // user-defined mobility at STP.
  DWORD_PTR               get_ion_mobility_ptr() const;
  void                    set_ion_mobility(double fMob);

  double                  get_act_energy() const;       // user-defined activation energy, eV, for output only.
  DWORD_PTR               get_act_energy_ptr() const;
  void                    set_act_energy(double fEact);

  double                  get_ion_cross_section() const;     // user-defined collision cross-section, cm^2, for output only.
  DWORD_PTR               get_ion_cross_section_ptr() const;
  void                    set_ion_cross_section(double fCrossSect);

  bool                    get_vel_depend_flag() const;
  DWORD_PTR               get_vel_depend_flag_ptr() const;
  void                    set_vel_depend_flag(bool bFlag);

  bool                    get_enable_rf() const;
  DWORD_PTR               get_enable_rf_ptr() const;
  void                    set_enable_rf(bool bEnable);

  double                  get_rf_amplitude() const;
  DWORD_PTR               get_rf_amplitude_ptr() const;
  void                    set_rf_amplitude(double fAmpl);

  double                  get_rf_frequency() const;
  DWORD_PTR               get_rf_frequency_ptr() const;
  void                    set_rf_frequency(double fFreq); // Hz

  double                  get_rf_Q00_ampl() const;
  DWORD_PTR               get_rf_Q00_ampl_ptr() const;
  void                    set_rf_Q00_ampl(double fAmpl);

  double                  get_rf_Q00_freq() const;
  DWORD_PTR               get_rf_Q00_freq_ptr() const;
  void                    set_rf_Q00_freq(double fFreq); // Hz

  double                  get_Q00_trans() const;
  DWORD_PTR               get_Q00_trans_ptr() const;
  void                    set_Q00_trans(double fX);      // cm

  double                  get_rf_flatapole_ampl() const;
  DWORD_PTR               get_rf_flatapole_ampl_ptr() const;
  void                    set_rf_flatapole_ampl(double fAmpl);

  double                  get_rf_flatapole_freq() const;
  DWORD_PTR               get_rf_flatapole_freq_ptr() const;
  void                    set_rf_flatapole_freq(double fFreq); // Hz

  double                  get_flatapole_trans() const;
  DWORD_PTR               get_flatapole_trans_ptr() const;
  void                    set_flatapole_trans(double fX);      // cm

// Coulomb effects:
  bool                    get_enable_coulomb() const;
  DWORD_PTR               get_enable_coulomb_ptr() const;
  void                    set_enable_coulomb(bool bEnable);

  double                  get_full_current() const;
  DWORD_PTR               get_full_current_ptr() const;
  void                    set_full_current(double fCurrent);

  double                  get_bunch_r0() const;
  DWORD_PTR               get_bunch_r0_ptr() const;
  void                    set_bunch_r0(double fR0);

  bool                    get_axial_symm() const;
  DWORD_PTR               get_axial_symm_ptr() const;
  void                    set_axial_symm(bool bAxialSymm);

// Coulomb repulsion in the absense of axial symmetry.
  UINT                    get_iter_count() const;
  DWORD_PTR               get_iter_count_ptr() const;
  void                    set_iter_count(UINT nCount);

  double                  get_BH_dist_par() const;
  DWORD_PTR               get_BH_dist_par_ptr() const;
  void                    set_BH_dist_par(double fDist);

  double                  get_crit_radius() const;
  DWORD_PTR               get_crit_radius_ptr() const;
  void                    set_crit_radius(double fR);

  UINT                    get_max_rec_depth() const;
  DWORD_PTR               get_max_rec_depth_ptr() const;
  void                    set_max_rec_depth(UINT nDepth);

  bool                    get_enable_quad_terms() const;
  DWORD_PTR               get_enable_quad_terms_ptr() const;
  void                    set_enable_quad_terms(bool bEnable);

  bool                    get_use_radial_coulomb() const;
  DWORD_PTR               get_use_radial_coulomb_ptr() const;
  void                    set_use_radial_coulomb(bool bEnable);

  double                  get_radial_coulomb_trans() const;
  DWORD_PTR               get_radial_coulomb_trans_ptr() const;
  void                    set_radial_coulomb_trans(double fX);

  bool                    get_use_pre_calc_coulomb() const;
  DWORD_PTR               get_use_pre_calc_coulomb_ptr() const;

  const char*             get_pre_calc_clmb_file() const;
  DWORD_PTR               get_pre_calc_clmb_file_ptr() const;
  void                    set_pre_calc_clmb_file(const char* pName);

// Random diffusion (for ion type of particles only).
  bool                    get_enable_diffusion() const;
  DWORD_PTR               get_enable_diffusion_ptr() const;

  CRandomDiffType         get_rand_diff_type() const;
  DWORD_PTR               get_rand_diff_type_ptr() const;
  void                    set_rand_diff_type(CRandomDiffType nType);

  long                    get_random_seed() const;
  DWORD_PTR               get_random_seed_ptr() const;
  void                    set_random_seed(long nSeed);

// Different integrators:
  bool                    get_use_old_integrator() const;
  DWORD_PTR               get_use_old_integrator_ptr() const;
  void                    set_use_old_integrator(bool bEnable);

  bool                    get_enable_ansys_field() const;
  DWORD_PTR               get_enable_ansys_field_ptr() const;

// Turn ON/OFF saving the tracks on disk:
  bool                    get_enable_save_tracks() const;
  DWORD_PTR               get_enable_save_tracks_ptr() const;
  void                    set_enable_save_tracks(bool bEnable);

//-------------------------------------------------------------------------------------------------
// Multi-threading support.
//-------------------------------------------------------------------------------------------------
  bool                    get_use_multi_thread() const;
  DWORD_PTR               get_use_multi_thread_ptr() const;
  void                    set_use_multi_thread(bool bEnable);

  void                    single_thread_calculate();
  void                    multi_thread_calculate();

// m_bTerminate becomes true in two cases: 1) after calling terminate() from the CExecutionDialog; 2) after posting WM_CLOSE 
// for the CExecutionDialog by main_thread_func when the calculation have ended successfully. Therefore, this flag can not be
// an indicator of successful run. Use m_bResult flag for this.
  bool                    get_result_flag() const;  // if m_bResult is true, the 3D scene will be drawn even after termination.

  void                    set_tracking_progress();  // works in multi-threading environment for progress of particle tracking only.

  void                    update_interface(); // force the properties list to update after loading the data.

//-------------------------------------------------------------------------------------------------
// Drawing support:
//-------------------------------------------------------------------------------------------------
  CNodesCollection&       get_nodes();
  CElementsCollection&    get_elems();
  CRegionsCollection&     get_regions();
  CTrackVector&           get_tracks();
  CBox&                   get_box();

  CSpaceChargeDistrib&    get_space_charge_dist();

  Vector3D                get_center() const;

  bool                    is_ready();

// OpenFOAM import support:
  CImportOpenFOAM&        get_importer();

// Output to files:
  COutputEngine&          get_output_engine();

// Mesh transformation:
  CTransform&             get_transform();

// DC Field perturbations:
  CFieldPtbCollection&    get_field_ptb();

//-------------------------------------------------------------------------------------------------
// Streaming:
//-------------------------------------------------------------------------------------------------
  void                    save(CArchive& archive);
  void                    load(CArchive& archive);

  void                    save_tracks(CArchive& archive);
  void                    load_tracks(CArchive& archive);

  void                    save_track_const(CArchive& ar, const CTrack& track);
  void                    load_track_const(CArchive& ar, CTrack& track);

  bool                    read_data();  // main function for reading ANSYS data.

  void                    clear_scene();

protected:
  void                    set_default();
  void                    set_data();         // force data from the properties list to be set before saving.

//-------------------------------------------------------------------------------------------------
// Mesh specific interface:
//-------------------------------------------------------------------------------------------------
  bool                    read_geometry();
  bool                    read_2D_regions();

public:
  bool                    read_gasdyn_data();

  bool                    save_coulomb_field(const char* pFile);
  bool                    read_coulomb_field();

protected:
  void                    add_tetra(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3);
  void                    add_pyramid(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4);
  void                    add_wedge(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4, CNode3D* p5);
  void                    add_hexa(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4, CNode3D* p5, CNode3D* p6, CNode3D* p7);

  void                    bounding_box();
  
  void                    clear();

public:
  bool                    interpolate(const Vector3D& vPos, double fTime, double fPhase, CNode3D& node, CElem3D*& pElem) const;

// Reflect the particle's position, velocity and acceleration against the symmetry plane(s)if necessary.
// The function returns "true" if reflection has been done and the back reflection is needed. In this case set
// bForceReflect to "true" in the next call to ensure that the back reflection will occur.
  bool                    sym_corr(Vector3D& vPos, Vector3D& vVel, Vector3D& vAccel, bool bForceReflect = false) const;

protected:
// Returns pointer to the element containing the input point or NULL if no element has been found.
  CElem3D*                find_elem(CElem3D* pPrevElem, const Vector3D& vPos) const;

// Try the nearest neighbors of the element, where the previous location was found, including the element itself.
// Return pointer to the element containing the input point or NULL if no element has been found.
  CElem3D*                try_neighbors(CElem3D* pElem, const Vector3D& vPos) const;

//-------------------------------------------------------------------------------------------------
// Particle tracking interface (everywhere the CGS system is suggested):
//-------------------------------------------------------------------------------------------------
  void                    initial_conditions();
  void                    clear_tracks(bool bFinally = true);

  void                    do_track();
  void                    do_iterations();  // start computing ion motion with Coulomb repulsion in a general case (no axial symmetery).

// Different calls for particles of droplet and ion types.
  bool                    do_time_step(CBaseTrackItem* pItem);
  bool                    do_ion_time_step(CBaseTrackItem* pItem, double fPhase, double fCurr);

  double                  get_Re(const Vector3D& vVel, double fDens, double fDynVisc, double fD) const;
  double                  get_Cd(double fRe) const;

// Droplet acceleration. Input: node - environment parameters, vVel - particle velocity, cm/s,
// fMass - particle mass, g, fD - particle diameter, cm, fTime - time since start of integration, s. 
// Output: acceleration of the particle, cm/s2, fRe - the Reynolds number.
  Vector3D                get_accel(const CNode3D&  node, 
                                    const Vector3D& vVel,
                                    double          fMass,
                                    double          fD,
                                    double          fTime,  // useful if only m_bEnableRF is ON.
                                    double&         fRe) const;

// Ion acceleration. The mass of the ion is constant and is the member of this class.
public:
  Vector3D                get_ion_accel(const CNode3D&  node,
                                        const Vector3D& vVel,
                                        double          fTime,
                                        double          fTimeStep,
                                        double          fPhase,
                                        double          fCurr,
                                        double&         fExpCoeff,
                                        double&         fMob) const;  // 16-05-2017 mobility for diffusion velocity jumps support.

  double                  get_dTi(const CNode3D& node, const Vector3D& vVel, double fIonTemp, double fExpCoeff, double& fTinf) const;

// Ion mobility recalculated from m_fIonMobility at STP using the Chapman-Enskog approximation.
  double                  get_ion_mob(double fPress, double fTemp) const;

protected:
  CBaseTrackItem*         create_track_item(size_t          nElemId,
                                            const Vector3D& vPos,
                                            const Vector3D& vVel,
                                            double          fMass,
                                            double          fTemp,
                                            double          fIonMob,  // random diffusion velocity jumps support.
                                            double          fTime) const;
// Coulomb effects:
  void                    init_currents();  // before integration runs over all track items and initializes item.curr members.

  void                    get_output_freq(UINT& nOutFreq) const;

  bool                    create_BH_object(CalcThreadVector& vThreads, UINT nIter);
  void                    get_BH_cube(Vector3D& c, double& edge, double& fMinX, double& fMaxX) const;

  double                  get_full_current_at(UINT nIter);  // the full current is gradually increasing with iteration number.
  bool                    capture_save_image(UINT nIter);   // capture and save screenshot with tracks after nIter-th iteration.

// These two functions are droplet-type specific:
  double                  get_particle_mass(double fD) const;
  double                  get_particle_diameter(double fMass) const;

// Track limitations:
  bool                    track_is_over(double fTime, double fMass, double fTemp) const;
  bool                    limit_of_Rayleigh(double fT, double fD) const;  // input: particle's temperature (K) and diameter (cm).
  double                  get_max_charge(double fT, double fD) const;

//-------------------------------------------------------------------------------------------------
// Multi-threading support:
//-------------------------------------------------------------------------------------------------
  static UINT __stdcall   main_thread_func(LPVOID pCalcThread);

  bool                    prepare();
  void                    relax();

//-------------------------------------------------------------------------------------------------
// Evaporation:
//-------------------------------------------------------------------------------------------------
  void                    create_evapor_model();

private:
  int                     m_nIntegrType;

// Multi-threading support.
  bool                    m_bUseMultiThread;

  CalcThreadVector*       m_pCalcThreads;   // run-time, for progress bar support only.

// Mesh specific data:
  std::string             m_sDataFile;

  CBox                    m_Box;          // bounding box.

  CNodesCollection        m_vNodes;
  CElementsCollection     m_vElems;

  CRegionsCollection      m_vRegions;     // Regions containing triangular faces, for drawing only.
  CRegionsCollection      m_vExtRegions;  // m_vRegions + cross-section regions, run-time.

  bool                    m_bConv2CGS;    // a flag showing whether ANSYS data, which are normally in SI must be 
                                          // converted to the CGS system; this is always "true" so far.
// Particle data:
  int                     m_nType;        // type of particles, either droplets or ions.

  CTrackVector            m_Tracks;       // output data.

// Saving the tracks on a hard disk:
  bool                    m_bSaveTracks;

public:
  double                  m_fTimeStep,
                          m_fHalfTimeStep,
                          m_fMaxIntegrTime;
private:
  bool                    m_bOldIntegrator, // the first predictor-corrector integrator improved by correct field averaging.
                          m_bAnsysFields;   // if true, the ANSYS calculated electric fields are used.

  double                  m_fInitMass,
                          m_fPartDens,
                          m_fInitD,
                          m_fMinD,      // minimal particle diameter, below which tracking stops automatically, cm ...
                          m_fMinMass,   // ... and corresponding minimal mass, g.

                          m_fMolarMass; // environment gas molar mass.

  int                     m_nSymPlanes;

// Electrostatics:
  bool                    m_bEnableField;
  double                  m_fAmplDC,    // DC voltage at the emitter, V. This is a multiplier equal to U/2500.
                          m_fCharge;    // sum particle's charge, CGSE.

  CSource                 m_Src;

// Evaporation data:
  CEvaporationModel*      m_pEvaporModel;
  int                     m_nEvaporModelType;

// Particles of Ion type:
  double                  m_fIonMass,         // g.
                          m_fChargeMassRatio, // CGSE.
                          m_fIonMobility,     // CGSE, this is a user-defined mobility at STP.
                          m_fCrossSection,    // cm^2, this is a user-defined collision cross-section with the environmental gas.
                          m_fActEnergy;       // activation energy of the ion, eV.

  bool                    m_bVelDependent;    // if true, the collision cross-section is inversely proportional to the relative velocity.

  bool                    m_bEnableRF;
// RF field stuff. The amplitude is in Volts because ANSYS calculated RF field for RF potential amplitude
// equal to 1 V. To scale the field to desirable voltage the m_fAmplRF coefficient must be in Volts.
  double                  m_fAmplRF,      // Amplitude and circular frequency of RF field in the funnel (S-Lens).
                          m_fOmega,

                          m_fAmplRF_Q00,  // Amplitude and circular frequency of RF field in the Q00 region (quadrupole, octopole).
                          m_fOmega_Q00,
                          m_fX_Q00,       // X-coordinate of the transition from the funnel (S-Lens) to Q00.

                          m_fAmplRF_Q0,   // Amplitude and circular frequency of RF field in the Q0 region (flatapole).
                          m_fOmega_Q0,
                          m_fX_Q0;        // X-coordinate of the transition from the Q00 region to Q0.

  void                    get_ampl_freq(const Vector3D& vPos, double& fAmpl, double& fFreq) const;

public:
// Obsolete, for backward compatibility only.
  Vector3D                get_rf_field(const CNode3D& node, double fTime, double fPhase) const;

protected:
  bool                    m_bEnableCoulomb,
                          m_bUsePreCalcCoulomb, // if true, look for the pre-calculated field in a file.
                          m_bAxialSymm;         // if "false" the iterational approach is applied.

  std::string             m_sClmbDataFile;    // file in which pre-calculated Coulomb field is stored.

  double                  m_fFullCurrent,     // full current in the bunch.
                          m_fInitBunchRadius; // initial radius of the bunch for axial symmetry, cm.

// Iterational calculations of the Coulomb repulsion, see also the CBarnesHut class:
  UINT                    m_nIterCount;       // count of iterations.

  double                  m_fBHDist,          // distance parameter for the Barnes-Hut object, dimensionless.
                          m_fBHCritR;         // critical radius for the Barnes-Hut object, cm; if r < m_fBHCritR, f = c / m_fBHCritR^3.

  UINT                    m_nMaxRecDepth;

  bool                    m_bEnableQuadTerms; // if true than the quadrupole terms are added to the field of a distant charge.

  CBarnesHut*             m_pBarnesHut;

  CSpaceChargeDistrib     m_SpaceChargeDist;
// If true, the Gabovich radial formula for the Coulomb field will be used for x > m_fRadialCoulombX.
  bool                    m_bUseRadialCoulomb;
  double                  m_fRadialCoulombX;

// DC Field perturbations:
  CFieldPtbCollection     m_vFieldPtbColl;

// Random diffusion (for ion type of particles only).
  bool                    m_bEnableDiffusion;
  CRandomDiffType         m_nRndDiffType;
  UINT                    m_nRandomSeed;

  RandomProcess*          create_random_jump(UINT nSeed) const;

// Convertation to CGS:
  void                    conv_to_cgs(float& fPress, float& fDens, float& fDynVisc, float& fThermCond, float& fCp,
                            float& fVx, float& fVy, float& fVz, float& fEx, float& fEy, float& fEz, float& fRFEx, float& fRFEy, float& fRFEz);
// Calculators:
  void                    invalidate_calculators(); // called from read_data() and makes the calculators update their internal data.

// Output stuff:
  COutputEngine           m_OutputEngine;

// Import OpenFOAM data:
  CImportOpenFOAM         m_Importer;

// Mesh transformation (applied in read_geometry()).
  CTransform              m_Transform;

// Run-time:
  bool                    m_bReady,     // this flag is set to "false" in set_filename(), to "true" in read_data().
                          m_bResult;

  bool                    abort(FILE* pStream = NULL);

  friend class CExportOpenFOAM;
  friend class CImportOpenFOAM;
  friend class COutputEngine;
  friend class CTrackDraw;
  friend class CSource;
};

//-------------------------------------------------------------------------------------------------
// Inline implementation:
//-------------------------------------------------------------------------------------------------
// Integrator type:
inline int CTracker::get_integr_type() const
{
  return m_nIntegrType;
}

inline DWORD_PTR CTracker::get_integr_type_ptr() const
{
  return (DWORD_PTR)&m_nIntegrType;
}

inline void CTracker::set_integr_type(int nType)
{
  m_nIntegrType = nType;
}

// Multi-threading support.
inline bool CTracker::get_use_multi_thread() const
{
  return m_bUseMultiThread;
}

inline DWORD_PTR CTracker::get_use_multi_thread_ptr() const
{
  return (DWORD_PTR)&m_bUseMultiThread;
}

inline void CTracker::set_use_multi_thread(bool bEnable)
{
  m_bUseMultiThread = bEnable;
}

// Gas-dynamic data input stuff.
inline const char* CTracker::get_filename() const
{
  return m_sDataFile.c_str();
}

inline DWORD_PTR CTracker::get_filename_ptr() const
{
  return (DWORD_PTR)&m_sDataFile;
}

inline bool CTracker::set_filename(const char* pName)
{
  if(strcmp(pName, m_sDataFile.c_str()) == 0)
    return false; // filename will not change, no need to set m_bReady flag to false. 

  m_sDataFile = pName;
  m_bReady = false; // force the data files to be re-read.
  return true;
}

inline bool CTracker::get_convert_to_cgs() const
{
  return m_bConv2CGS;
}

inline void CTracker::set_convert_to_cgs(bool bEnable)
{
  m_bConv2CGS = bEnable;
}

inline int CTracker::get_particle_type() const
{
  return m_nType;
}

inline DWORD_PTR CTracker::get_particle_type_ptr() const
{
  return (DWORD_PTR)&m_nType;
}

inline void CTracker::set_particle_type(int nType)
{
  m_nType = nType;
}

inline double CTracker::get_molar_mass() const
{
  return m_fMolarMass;
}

inline DWORD_PTR CTracker::get_molar_mass_ptr() const
{
  return (DWORD_PTR)&m_fMolarMass;
}

inline void CTracker::set_molar_mass(double fMass)
{
  m_fMolarMass = fMass;
}

inline double CTracker::get_time_step() const
{
  return m_fTimeStep;
}

inline DWORD_PTR CTracker::get_time_step_ptr() const
{
  return (DWORD_PTR)&m_fTimeStep;
}

inline void CTracker::set_time_step(double fTimeStep)
{
  m_fTimeStep = fTimeStep;
  m_fHalfTimeStep = 0.5 * m_fTimeStep;
}

inline double CTracker::get_max_track_time() const
{
  return m_fMaxIntegrTime;
}

inline DWORD_PTR CTracker::get_max_track_time_ptr() const
{
  return (DWORD_PTR)&m_fMaxIntegrTime;
}

inline void CTracker::set_max_track_time(double fTime)
{
  m_fMaxIntegrTime = fTime;
}

// Electrostatics:
inline bool CTracker::get_enable_field() const
{
  return m_bEnableField;
}

inline DWORD_PTR CTracker::get_enable_field_ptr() const
{
  return (DWORD_PTR)&m_bEnableField;
}

inline void CTracker::set_enable_field(bool bEnable)
{
  m_bEnableField = bEnable;
}

inline double CTracker::get_particle_charge() const
{
  return m_fCharge;
}

inline DWORD_PTR CTracker::get_particle_charge_ptr() const
{
  return (DWORD_PTR)&m_fCharge;
}

inline void CTracker::set_particle_charge(double fCharge)
{
  m_fCharge = fCharge;
  m_fChargeMassRatio = m_fCharge / m_fIonMass;
}

inline double CTracker::get_dc_amplitude() const
{
  return m_fAmplDC * cfDCVoltage;
}

inline DWORD_PTR CTracker::get_dc_amplitude_ptr() const
{
  return (DWORD_PTR)&m_fAmplDC;
}

inline void CTracker::set_dc_amplitude(double fAmpl)
{
  m_fAmplDC = fAmpl / cfDCVoltage;
}

inline int CTracker::get_sym_plane() const
{
  return m_nSymPlanes;
}

inline DWORD_PTR CTracker::get_sym_plane_ptr() const
{
  return (DWORD_PTR)&m_nSymPlanes;
}

inline void CTracker::set_sym_plane(int nPlane)
{
  m_nSymPlanes = nPlane;
}

// Particle source parameters:
inline CSource* CTracker::get_src()
{
  return &m_Src;
}

// Particles evaporation:
inline int CTracker::get_evapor_model_type() const
{
  return m_nEvaporModelType;
}

inline DWORD_PTR CTracker::get_evapor_model_type_ptr() const
{
  return (DWORD_PTR)&m_nEvaporModelType;
}

inline void CTracker::set_evapor_model_type(int nModelType)
{
  m_nEvaporModelType = nModelType;
  create_evapor_model();
}

inline CEvaporationModel* CTracker::get_evapor_model() const
{
  return m_pEvaporModel;
}

inline double CTracker::get_init_diameter() const
{
  return m_fInitD;
}

inline DWORD_PTR CTracker::get_init_diameter_ptr() const
{
  return (DWORD_PTR)&m_fInitD;
}

inline void CTracker::set_init_diameter(double fD)
{
  m_fInitD = fD;
  m_fInitMass = get_particle_mass(fD);
}

// Ion type of particles:
inline double CTracker::get_ion_mass() const
{
  return m_fIonMass;
}

inline DWORD_PTR CTracker::get_ion_mass_ptr() const
{
  return (DWORD_PTR)&m_fIonMass;
}

inline void CTracker::set_ion_mass(double fMass)
{
  m_fIonMass = fMass;
  m_fChargeMassRatio = m_fCharge / m_fIonMass;
}

inline double CTracker::get_ion_mobility() const
{
  return m_fIonMobility;
}

inline DWORD_PTR CTracker::get_ion_mobility_ptr() const
{
  return (DWORD_PTR)&m_fIonMobility;
}

inline void CTracker::set_ion_mobility(double fMob)
{
  m_fIonMobility = fMob;
}

inline double CTracker::get_ion_mob(double fPress, double fTemp) const
{
  return m_fIonMobility * sqrt(fTemp / 273.) * (Const_One_Atm_CGS / fPress);
}

inline double CTracker::get_ion_cross_section() const     // user-defined collision cross-section, cm^2.
{
  return m_fCrossSection;
}

inline DWORD_PTR CTracker::get_ion_cross_section_ptr() const
{
  return (DWORD_PTR)&m_fCrossSection;
}

inline void CTracker::set_ion_cross_section(double fCrossSect)
{
  m_fCrossSection = fCrossSect;
}

inline bool CTracker::get_vel_depend_flag() const
{
  return m_bVelDependent;
}

inline DWORD_PTR CTracker::get_vel_depend_flag_ptr() const
{
  return (DWORD_PTR)&m_bVelDependent;
}

inline void CTracker::set_vel_depend_flag(bool bFlag)
{
  m_bVelDependent = bFlag;
}

inline double CTracker::get_act_energy() const
{
  return m_fActEnergy;
}

inline DWORD_PTR CTracker::get_act_energy_ptr() const
{
  return (DWORD_PTR)&m_fActEnergy;
}

inline void CTracker::set_act_energy(double fEact)
{
  m_fActEnergy = fEact;
}

inline bool CTracker::get_enable_rf() const
{
  return m_bEnableRF;
}

inline DWORD_PTR CTracker::get_enable_rf_ptr() const
{
  return (DWORD_PTR)&m_bEnableRF;
}

inline void CTracker::set_enable_rf(bool bEnable)
{
  m_bEnableRF = bEnable;
}

inline double CTracker::get_rf_amplitude() const
{
  return m_fAmplRF;
}

inline DWORD_PTR CTracker::get_rf_amplitude_ptr() const
{
  return (DWORD_PTR)&m_fAmplRF;
}

inline void CTracker::set_rf_amplitude(double fAmpl)
{
  m_fAmplRF = fAmpl;
}

inline double CTracker::get_rf_frequency() const
{
  return m_fOmega / Const_2PI;
}

inline DWORD_PTR CTracker::get_rf_frequency_ptr() const
{
  return (DWORD_PTR)&m_fOmega;
}

inline void CTracker::set_rf_frequency(double fFreq)
{
  m_fOmega = Const_2PI * fFreq;
}

inline double CTracker::get_rf_Q00_ampl() const
{
  return m_fAmplRF_Q00;
}

inline DWORD_PTR CTracker::get_rf_Q00_ampl_ptr() const
{
  return (DWORD_PTR)&m_fAmplRF_Q00;
}

inline void CTracker::set_rf_Q00_ampl(double fAmpl)
{
  m_fAmplRF_Q00 = fAmpl;
}

inline double CTracker::get_rf_Q00_freq() const
{
  return m_fOmega_Q00 / Const_2PI;
}

inline DWORD_PTR CTracker::get_rf_Q00_freq_ptr() const
{
  return (DWORD_PTR)&m_fOmega_Q00;
}

inline void CTracker::set_rf_Q00_freq(double fFreq)
{
  m_fOmega_Q00 = Const_2PI * fFreq;
}

inline double CTracker::get_Q00_trans() const
{
  return m_fX_Q00;
}

inline DWORD_PTR CTracker::get_Q00_trans_ptr() const
{
  return (DWORD_PTR)&m_fX_Q00;
}

inline void CTracker::set_Q00_trans(double fX)
{
  m_fX_Q00 = fX;
}

inline double CTracker::get_rf_flatapole_ampl() const
{
  return m_fAmplRF_Q0;
}

inline DWORD_PTR CTracker::get_rf_flatapole_ampl_ptr() const
{
  return (DWORD_PTR)&m_fAmplRF_Q0;
}

inline void CTracker::set_rf_flatapole_ampl(double fAmpl)
{
  m_fAmplRF_Q0 = fAmpl;
}

inline double CTracker::get_rf_flatapole_freq() const
{
  return m_fOmega_Q0 / Const_2PI;
}

inline DWORD_PTR CTracker::get_rf_flatapole_freq_ptr() const
{
  return (DWORD_PTR)&m_fOmega_Q0;
}

inline void CTracker::set_rf_flatapole_freq(double fFreq)
{
  m_fOmega_Q0 = Const_2PI * fFreq;
}

inline double CTracker::get_flatapole_trans() const
{
  return m_fX_Q0;
}

inline DWORD_PTR CTracker::get_flatapole_trans_ptr() const
{
  return (DWORD_PTR)&m_fX_Q0;
}

inline void CTracker::set_flatapole_trans(double fX)
{
  m_fX_Q0 = fX;
}

inline void CTracker::get_ampl_freq(const Vector3D& vPos, double& fAmpl, double& fFreq) const
{
  if(vPos.x < m_fX_Q00) // in the funnel (S-Lens).
  {
    fAmpl = m_fAmplRF;
    fFreq = m_fOmega;
  }
  else if(vPos.x < m_fX_Q0) // in the Q00 region.
  {
    fAmpl = m_fAmplRF_Q00;
    fFreq = m_fOmega_Q00;
  }
  else // in the Q0 (flatapole) region.
  {
    fAmpl = m_fAmplRF_Q0;
    fFreq = m_fOmega_Q0;
  }
}

// Coulomb effects:
inline bool CTracker::get_enable_coulomb() const
{
  return m_bEnableCoulomb;
}

inline DWORD_PTR CTracker::get_enable_coulomb_ptr() const
{
  return (DWORD_PTR)&m_bEnableCoulomb;
}

inline void CTracker::set_enable_coulomb(bool bEnable)
{
  m_bEnableCoulomb = bEnable;
}

inline double CTracker::get_full_current() const
{
  return m_fFullCurrent;
}

inline DWORD_PTR CTracker::get_full_current_ptr() const
{
  return (DWORD_PTR)&m_fFullCurrent;
}

inline void CTracker::set_full_current(double fCurrent)
{
  m_fFullCurrent = fCurrent;
}

inline double CTracker::get_bunch_r0() const
{
  return m_fInitBunchRadius;
}

inline DWORD_PTR CTracker::get_bunch_r0_ptr() const
{
  return (DWORD_PTR)&m_fInitBunchRadius;
}

inline void CTracker::set_bunch_r0(double fR0)
{
  m_fInitBunchRadius = fR0;
}

inline bool CTracker::get_axial_symm() const
{
  return m_bAxialSymm;
}

inline DWORD_PTR CTracker::get_axial_symm_ptr() const
{
  return (DWORD_PTR)&m_bAxialSymm;
}

inline void CTracker::set_axial_symm(bool bAxialSymm)
{
  m_bAxialSymm = bAxialSymm;
}

inline UINT CTracker::get_iter_count() const
{
  return m_nIterCount;
}

inline DWORD_PTR CTracker::get_iter_count_ptr() const
{
  return (DWORD_PTR)&m_nIterCount;
}

inline void CTracker::set_iter_count(UINT nCount)
{
  m_nIterCount = nCount;
}

inline double CTracker::get_BH_dist_par() const
{
  return m_fBHDist;
}

inline DWORD_PTR CTracker::get_BH_dist_par_ptr() const
{
  return (DWORD_PTR)&m_fBHDist;
}

inline void CTracker::set_BH_dist_par(double fDist)
{
  m_fBHDist = fDist;
}

inline double CTracker::get_crit_radius() const
{
  return m_fBHCritR;
}
inline DWORD_PTR CTracker::get_crit_radius_ptr() const
{
  return (DWORD_PTR)&m_fBHCritR;
}
inline void CTracker::set_crit_radius(double fR)
{
  m_fBHCritR = fR;
}

inline UINT CTracker::get_max_rec_depth() const
{
  return m_nMaxRecDepth;
}

inline DWORD_PTR CTracker::get_max_rec_depth_ptr() const
{
  return (DWORD_PTR)&m_nMaxRecDepth;
}

inline void CTracker::set_max_rec_depth(UINT nDepth)
{
  m_nMaxRecDepth = nDepth;
}

inline bool CTracker::get_enable_quad_terms() const
{
  return m_bEnableQuadTerms;
}

inline DWORD_PTR CTracker::get_enable_quad_terms_ptr() const
{
  return (DWORD_PTR)&m_bEnableQuadTerms;
}

inline void CTracker::set_enable_quad_terms(bool bEnable)
{
  m_bEnableQuadTerms = bEnable;
}

inline bool CTracker::get_use_radial_coulomb() const
{
  return m_bUseRadialCoulomb;
}

inline DWORD_PTR CTracker::get_use_radial_coulomb_ptr() const
{
  return (DWORD_PTR)&m_bUseRadialCoulomb;
}

inline void CTracker::set_use_radial_coulomb(bool bEnable)
{
  m_bUseRadialCoulomb = bEnable;
}

inline double CTracker::get_radial_coulomb_trans() const
{
  return m_fRadialCoulombX;
}

inline DWORD_PTR CTracker::get_radial_coulomb_trans_ptr() const
{
  return (DWORD_PTR)&m_fRadialCoulombX;
}

inline void CTracker::set_radial_coulomb_trans(double fX)
{
  m_fRadialCoulombX = fX;
}

inline bool CTracker::get_use_pre_calc_coulomb() const
{
  return m_bUsePreCalcCoulomb;
}

inline DWORD_PTR CTracker::get_use_pre_calc_coulomb_ptr() const
{
  return (DWORD_PTR)&m_bUsePreCalcCoulomb;
}

inline const char* CTracker::get_pre_calc_clmb_file() const
{
  return m_sClmbDataFile.c_str();
}

inline DWORD_PTR CTracker::get_pre_calc_clmb_file_ptr() const
{
  return (DWORD_PTR)&m_sClmbDataFile;
}

inline void CTracker::set_pre_calc_clmb_file(const char* pName)
{
  m_sClmbDataFile = pName;
}

// Import of gas-dynamic data support:
inline CImportOpenFOAM& CTracker::get_importer()
{
  return m_Importer;
}

inline void CTracker::conv_to_cgs(float& fPress, float& fDens, float& fDynVisc, float& fThermCond, float& fCp,
  float& fVx, float& fVy, float& fVz, float& fEx, float& fEy, float& fEz, float& fRFEx, float& fRFEy, float& fRFEz)
{
  fPress *= (float)SI_to_CGS_Press;
  fDens *= (float)SI_to_CGS_Dens;
  fDynVisc *= (float)SI_to_CGS_DynVisc;
  fThermCond *= (float)SI_to_CGS_ThermCond;

  fEx *= (float)SI_to_CGS_ElecField;
  fEy *= (float)SI_to_CGS_ElecField;
  fEz *= (float)SI_to_CGS_ElecField;

  fRFEx *= (float)SI_to_CGS_ElecField;
  fRFEy *= (float)SI_to_CGS_ElecField;
  fRFEz *= (float)SI_to_CGS_ElecField;

  fVx *= (float)SI_to_CGS_Vel;
  fVy *= (float)SI_to_CGS_Vel;
  fVz *= (float)SI_to_CGS_Vel;

  fCp *= (float)SI_to_CGS_Cp;
}

inline double CTracker::get_Re(const Vector3D& vVel, double fDens, double fDynVisc, double fD) const
{
  return fDens * fD * vVel.length() / fDynVisc;
}

// ANSYS suggested expression leading to constant Cd at high Re. At low Re Stokes drag coefficient is obtained.
inline double CTracker::get_Cd(double fRe) const
{
  if(fRe < Const_Almost_Zero)
    return 0.;

  return max((24. / fRe) * (1. + 0.15 * exp(0.687 * log(fRe))), 0.44);
}

inline double CTracker::get_particle_mass(double fD) const
{
  return Const_One_Sixth * Const_PI * fD * fD * fD * m_fPartDens;
}

inline double CTracker::get_particle_diameter(double fMass) const
{
  return exp(Const_One_Third * log(fMass * 6. * Const_1_PI / m_fPartDens));
}

inline CElementsCollection& CTracker::get_elems()
{
  return m_vElems;
}

// Drawing support:
inline CNodesCollection& CTracker::get_nodes()
{
  return m_vNodes;
}

inline CTrackVector& CTracker::get_tracks()
{
  return m_Tracks;
}

inline CBox& CTracker::get_box()
{
  return m_Box;
}

inline Vector3D CTracker::get_center() const
{
  return m_Box.get_center();
}

inline bool CTracker::is_ready()
{
  return m_bReady;
}

inline bool CTracker::get_result_flag() const
{
  return m_bResult;
}

// Output stuff:
inline COutputEngine& CTracker::get_output_engine()
{
  return m_OutputEngine;
}

// Different integrators:
inline bool CTracker::get_use_old_integrator() const
{
  return m_bOldIntegrator;
}

inline DWORD_PTR CTracker::get_use_old_integrator_ptr() const
{
  return (DWORD_PTR)&m_bOldIntegrator;
}

inline void CTracker::set_use_old_integrator(bool bEnable)
{
  m_bOldIntegrator = bEnable;
}

inline bool CTracker::get_enable_ansys_field() const
{
  return m_bAnsysFields;
}

inline DWORD_PTR CTracker::get_enable_ansys_field_ptr() const
{
  return (DWORD_PTR)&m_bAnsysFields;
}

inline bool CTracker::get_enable_save_tracks() const
{
  return m_bSaveTracks;
}

inline DWORD_PTR CTracker::get_enable_save_tracks_ptr() const
{
  return (DWORD_PTR)&m_bSaveTracks;
}

inline void CTracker::set_enable_save_tracks(bool bEnable)
{
  m_bSaveTracks = bEnable;
}

inline CSpaceChargeDistrib& CTracker::get_space_charge_dist()
{
  return m_SpaceChargeDist;
}

inline CTransform& CTracker::get_transform()
{
  return m_Transform;
}

inline CFieldPtbCollection& CTracker::get_field_ptb()
{
  return m_vFieldPtbColl;
}

// Random diffusion (for ion type of particles only).
inline bool CTracker::get_enable_diffusion() const
{
  return m_bEnableDiffusion;
}

inline DWORD_PTR CTracker::get_enable_diffusion_ptr() const
{
  return (DWORD_PTR)&m_bEnableDiffusion;
}

inline long CTracker::get_random_seed() const
{
  return m_nRandomSeed;
}

inline DWORD_PTR CTracker::get_random_seed_ptr() const
{
  return (DWORD_PTR)&m_nRandomSeed;
}

inline void CTracker::set_random_seed(long nSeed)
{
  m_nRandomSeed = UINT(std::abs(nSeed));
}

inline CRandomDiffType CTracker::get_rand_diff_type() const
{
  return m_nRndDiffType;
}

inline DWORD_PTR CTracker::get_rand_diff_type_ptr() const
{
  return (DWORD_PTR)&m_nRndDiffType;
}

inline void CTracker::set_rand_diff_type(CRandomDiffType nType)
{
  m_nRndDiffType = nType;
}

}; // namespace EvaporatingParticle

#endif