#pragma once
#ifndef _BatchSim_
#define _BatchSim_

#include "vector3d.hpp"
#include <vector>

namespace EvaporatingParticle
{

class CTracker;
class CBarnesHut;
class CElectricFieldData;
class CTrackCalculator;
struct CNode3D;

struct CFieldParams
{
  double  fAmpl,  // scale, in V.
          fFreq;  // RF frequency in Hz.
};

struct CIonSimParams  // all 5 double parameters are in CGSE.
{
  CIonSimParams();

  CIonSimParams(double fM, double fQ, double fK0, double fCS, double fCurr)
    : fMass(fM), fCharge(fQ), fMob(fK0), fCrossSect(fCS), fFullCurr(fCurr), sLegend(" ")
  {
  }

  double  fMass,
          fCharge,
          fMob,       // at STP.
          fCrossSect,
          fFullCurr;

  std::vector<CFieldParams> vFieldPar;

// A string that is added to the end of any file name saying what parameter is being varied. 
// Note: Each string in the input batch file must end with the legend!
  CString sLegend;
};

typedef std::vector<CIonSimParams> CSimParamsVector;
//-------------------------------------------------------------------------------------------------
// CBatchSim - an auxiliary class for batch simulations. 
//-------------------------------------------------------------------------------------------------
class CBatchSim
{
public:
  CBatchSim();
  ~CBatchSim();

  bool                  get_enable() const;
  DWORD_PTR             get_enable_ptr() const;

  CString               get_filename() const;
  DWORD_PTR             get_filename_ptr() const;
  void                  set_filename(const CString& sName);

  bool                  get_use_sliding_aver() const;
  DWORD_PTR             get_use_sliding_aver_ptr() const;

  double                get_curr_incr_iter() const;
  DWORD_PTR             get_curr_incr_iter_ptr() const;
  void                  set_curr_incr_iter(double fCurrIncr);

  UINT                  get_aver_width() const;
  DWORD_PTR             get_aver_width_ptr() const;
  void                  set_aver_width(UINT nWidth);

// The following 6 functions are called from CPropertiesWnd::OnDoTracking()
  bool                  batch_calc_init();
  void                  batch_calc_relax();

  UINT                  get_batch_steps_count() const;
  bool                  prepare_step(UINT nStep);

  bool                  intermediate_output();  // calculators outputs and Coulomb field output.
  bool                  intermediate_save();    // intermediate project save with simulation results on the given step.

  void                  accum_clmb_field_in_node(CNode3D* pNode, CBarnesHut* pBHObj, CElectricFieldData* pData, UINT nIter);

// Streaming:
  void                  save(CArchive& archive);
  void                  load(CArchive& archive);

protected:
  void                  set_default();

  bool                  read_sim_params();

  void                  backup_sim_params();
  bool                  set_sim_params_to_tracker(const CIonSimParams& param);

  bool                  get_filenames(CString& sTransFile, CString& sFragmFile, CString& sIonTempFile, CString& sClmbFile);

// Find the CTrackCalculator calculator among the calculators collection:
  bool                  get_calc();
  void                  release_calc();

  void                  init_prev_clmb();
  void                  clear_prev_clmb();

// DEBUG
  void                  debug_output();
// END DEBUG

private:
  bool                  m_bEnable;

  double                m_fIncrPerIter;   // user-defined average current increment per iteration.
  
  CString               m_sFileName;    // the full path to the file containing a list of the ion parameters for batch simulations.

  UINT                  m_nStepsCount;  // count of strings in the simulation list; run-time, defined after reading the input file.

  CTracker*             m_pObj;
  CTrackCalculator*     m_pCalc;

  CSimParamsVector      m_vSimParams;     // an array of simulation parameters for batch calculations.
  CIonSimParams         m_BackUpParams;   // parameters of the CTracker object and the Fields Collection to be stored before batch simulations.

  CString               m_sTskFile,       // run-time project file name for intermediate output.
                        m_sLegend;        // any intermediate file name must end up with this legend, which is set manually by user in the input batch file.

// These variables are not intended for frequent usage. In the latest version of the CTracker the sliding average is not supported yet.
  bool                  m_bUseSlidingAverage;

  UINT                  m_nAverWidth;     // dimension of the Vector3F array for coulomb field averaging.
// These two arrays are allocated only if m_bUseSlidingAverage is true.
  Vector3F**            m_pPrevClmb;  // two-dimensional array of Coulomb field vector from previous iterations.
  float**               m_pPrevPhi;   // two-dimensional array of Coulomb potential (including mirror if any).

  double                m_fNormCoeff; // run-time, m_fNormCoeff = 1 / (1 + m_nAverWidth).
};

//-------------------------------------------------------------------------------------------------
// Inline implementation
//-------------------------------------------------------------------------------------------------
inline bool CBatchSim::get_enable() const
{
  return m_bEnable;
}

inline DWORD_PTR CBatchSim::get_enable_ptr() const
{
  return (DWORD_PTR)&m_bEnable;
}

inline CString CBatchSim::get_filename() const
{
  return m_sFileName;
}

inline DWORD_PTR CBatchSim::get_filename_ptr() const
{
  return (DWORD_PTR)&m_sFileName;
}

inline void CBatchSim::set_filename(const CString& sName)
{
  m_sFileName = sName;
}

inline UINT CBatchSim::get_batch_steps_count() const
{
  return m_nStepsCount; // this is a run-time variable, it is defined in CBatchSim::read_sim_params(). 
}

inline bool CBatchSim::get_use_sliding_aver() const
{
  return m_bUseSlidingAverage;
}

inline DWORD_PTR CBatchSim::get_use_sliding_aver_ptr() const
{
  return (DWORD_PTR)&m_bUseSlidingAverage;
}

inline double CBatchSim::get_curr_incr_iter() const
{
  return m_fIncrPerIter;
}

inline DWORD_PTR CBatchSim::get_curr_incr_iter_ptr() const
{
  return (DWORD_PTR)&m_fIncrPerIter;
}

inline UINT CBatchSim::get_aver_width() const
{
  return m_nAverWidth;
}

inline DWORD_PTR CBatchSim::get_aver_width_ptr() const
{
  return (DWORD_PTR)&m_nAverWidth;
}

inline void CBatchSim::set_aver_width(UINT nWidth)
{
  m_nAverWidth = nWidth;
  m_fNormCoeff = 1. / (1 + m_nAverWidth);
}

};  // namespace EvaporatingParticle

#endif
