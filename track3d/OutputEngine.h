
#pragma once

#include "CObject.h"
#include "Elements.h"
#include "TrackItem.h"
#include <vector>
#include <string>


namespace EvaporatingParticle
{

/**
 * CTrackItem was moved to a separate file by AC 25.06.2016
 */

//-------------------------------------------------------------------------------------------------
// COutputEngine - a class comprising all the output stuff. 
//-------------------------------------------------------------------------------------------------
class COutputEngine : public CObject
{
public:
  COutputEngine();
  virtual ~COutputEngine();

// User interface
  const char*       get_out_dir() const;
  DWORD_PTR         get_out_dir_ptr() const;
  bool              set_out_dir(const char* pPath);

  double            get_output_time_step() const;
  DWORD_PTR         get_output_time_step_ptr() const;
  void              set_output_time_step(double fOutTimeStep);

  bool              get_enable_file_output() const;
  DWORD_PTR         get_enable_file_output_ptr() const;
  void              set_enable_file_output(bool bEnable);

  bool              get_restrict_output() const;          // restrict the file output by only those tracks which have
  DWORD_PTR         get_restrict_output_ptr() const;      // passed to the Q00 region, m_bOnlyPassedQ00 flag.
  void              set_restrict_output(bool bRestrict);

// Output types:
  bool              get_enable_ens_by_radius() const;
  DWORD_PTR         get_enable_ens_by_radius_ptr() const;
  void              set_enable_ens_by_radius(bool bEnable);

  UINT              get_ens_by_radius_count() const;
  DWORD_PTR         get_ens_by_radius_count_ptr() const;
  void              set_ens_by_radius_count(UINT nCount);

  double            get_cross_section_x() const;
  DWORD_PTR         get_cross_section_x_ptr() const;
  void              set_cross_section_x(double fX);

  static std::string get_full_path(const char* pFullPath);    // returns the full path without the file name itself.
  static std::string get_file_name(const char* pFullPath);    // returns the file name without the path.
  static std::string get_base_name(const char* pFullPath);    // discards the extension only, keeps the rest path.

// Streams:
  void              save(CArchive& archive);
  void              load(CArchive& archive);

  bool              output_droplet_track(UINT nTrackIndex);

  void              output_ion_tracks();

  void              output_ion_ensemble(const CTrackVector& vEnsTracks, FILE* pStream, int nPercentStart, int nPercentEnd);
  UINT              get_max_ensemble_index() const;

  void              output_ensemble_by_ens_index();        // ion ensembles are built by ensemble index.
  void              output_ensemble_by_initial_radius();   // ion ensembles are built by initial radii of the ions.

  void              output_average_tracks();
  void              get_output_range(double& fArgMin, double& fArgMax);

  void              output_ion_current();
  void              output_axial_current(const char* cFileName);
  void              output_averaged_current(const char* cFileName);

// Averaging the output current over latest iterations:
  void              prepare_current_output(); // this function must be called just after the zero iteration to set up output ranges.
  void              add_current();

// Distribution of trajectories in a cross-section:
  void              out_cross_section();

protected:
  void              set_default();

private:
  std::string       m_sOutDir;

  double            m_fOutputTimeStep;  // s.

  bool              m_bEnableFileOutput,
                    m_bOnlyPassedQ00,   // output to files only those tracks which have passed to Q00 region.

// Build ensembles by the starting radii, not by ensemble indices. m_Src.m_nEnsembleSize is supposed to be unity if m_bByInitRadii is true.
                    m_bByInitRadii;

  UINT              m_nEnsByRadiusCount;  // count of ranges of initial radii.

// Distribution of trajectories in a cross-section:
  double            m_fCrssSctX;

  CAveragingEngine  m_AverEngine;
};

//-------------------------------------------------------------------------------------------------
// COutputEngine - inline implementation. 
//-------------------------------------------------------------------------------------------------
inline const char* COutputEngine::get_out_dir() const
{
  return m_sOutDir.c_str();
}

inline DWORD_PTR COutputEngine::get_out_dir_ptr() const
{
  return (DWORD_PTR)&m_sOutDir;
}

inline bool COutputEngine::set_out_dir(const char* pPath)
{
  m_sOutDir = pPath;
}

inline double COutputEngine::get_output_time_step() const
{
  return m_fOutputTimeStep;
}

inline DWORD_PTR COutputEngine::get_output_time_step_ptr() const
{
  return (DWORD_PTR)&m_fOutputTimeStep;
}

inline void COutputEngine::set_output_time_step(double fOutTimeStep)
{
  m_fOutputTimeStep = fOutTimeStep;
}

inline bool COutputEngine::get_enable_file_output() const
{
  return m_bEnableFileOutput;
}

inline DWORD_PTR COutputEngine::get_enable_file_output_ptr() const
{
  return (DWORD_PTR)&m_bEnableFileOutput;
}

inline void COutputEngine::set_enable_file_output(bool bEnable)
{
  m_bEnableFileOutput = bEnable;
}

inline bool COutputEngine::get_restrict_output() const
{
  return m_bOnlyPassedQ00;
}

inline DWORD_PTR COutputEngine::get_restrict_output_ptr() const
{
  return (DWORD_PTR)&m_bOnlyPassedQ00;
}

inline void COutputEngine::set_restrict_output(bool bRestrict)
{
  m_bOnlyPassedQ00 = bRestrict;
}

// File output types:
inline bool COutputEngine::get_enable_ens_by_radius() const
{
  return m_bByInitRadii;
}

inline DWORD_PTR COutputEngine::get_enable_ens_by_radius_ptr() const
{
  return (DWORD_PTR)&m_bByInitRadii;
}

inline void COutputEngine::set_enable_ens_by_radius(bool bEnable)
{
  m_bByInitRadii = bEnable;
}

inline UINT COutputEngine::get_ens_by_radius_count() const
{
  return m_nEnsByRadiusCount;
}

inline DWORD_PTR COutputEngine::get_ens_by_radius_count_ptr() const
{
  return (DWORD_PTR)&m_nEnsByRadiusCount;
}

inline void COutputEngine::set_ens_by_radius_count(UINT nCount)
{
  m_nEnsByRadiusCount = nCount;
}

inline double COutputEngine::get_cross_section_x() const
{
  return m_fCrssSctX;
}

inline DWORD_PTR COutputEngine::get_cross_section_x_ptr() const
{
  return (DWORD_PTR)&m_fCrssSctX;
}

inline void COutputEngine::set_cross_section_x(double fX)
{
  m_fCrssSctX = fX;
}

};  // namespace EvaporatingParticle