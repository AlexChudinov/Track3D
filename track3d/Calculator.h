
#pragma once

#include "CObject.h"
#include "Elements.h"
#include "ColorContour.h"

#include <vector>
#include <string>

namespace EvaporatingParticle
{

class CTracker;
//-------------------------------------------------------------------------------------------------
// CCalculator - the base class designed for calculations on arbitrary scalar or vector fields.
//-------------------------------------------------------------------------------------------------
class CCalculator : public CObject
{
public:
  CCalculator();
  virtual ~CCalculator();

  enum  // calculated variable.
  {
    clcMassFlow   = 0,
    clcVolumeFlow = 1,
    clcEnergyFlow = 2,
    clcAvePress   = 3,
    clcAveTemp    = 4,
    clcAveDens    = 5,
    clcAveVx      = 6,
    clcAveRe      = 7,  // Reynolds number averaged over the cross-section.
    clcCount      = 8
  };

  enum  // calculator type.
  {
    ctPlaneYZ        = 0,
    ctSelRegions     = 1,
    ctAlongLine      = 2,
    ctTrackCalc      = 3,
    ctTrackCrossSect = 4,
    ctCount          = 5
  };

  static CCalculator* create(int nType);
  static const char*  calc_name(int nType);

  virtual void        run() = 0;
  virtual void        do_calculate() = 0;
  virtual void        update() = 0;
  virtual void        clear() = 0;

  virtual int         type() const = 0;

  virtual const char* units() const;
  virtual const char* get_var_name(int nVar) const;

  virtual bool        get_update_flag() const;  // calculators which must re-calculate at every change return true.
  virtual int         calc_vars_count() const;

  void                invalidate();

  bool                get_enable() const;
  DWORD_PTR           get_enable_ptr() const;
  void                set_enable(bool bEnable);

  int                 get_clc_var_type() const;
  DWORD_PTR           get_clc_var_type_ptr() const;
  void                set_clc_var_type(int nClcVar);

  double              get_char_length() const;
  DWORD_PTR           get_char_length_ptr() const;
  void                set_char_length(double fLen);

  double              get_result() const;
  DWORD_PTR           get_result_ptr() const;

  const char*         get_name() const;

//-------------------------------------------------------------------------------------------------
// Streaming:
//-------------------------------------------------------------------------------------------------
  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

protected:
  bool                get_tracker_ptr();

  bool                intersect(const CRay& ray, CIntersectColl& results) const;

  bool                process_face(CFace* pFace, double& fSquare, double& fRes) const;

  void                normalize_result(double fSumSquare, double fSumRes);
  virtual void        convert_result();   // convertion to SI units (or other output units).

// User-defined:
  bool                m_bEnable;
  int                 m_nClcVar;

  double              m_fCharLength;   // specially for Reynolds number calculation.

  std::string         m_sCalcName;

// Run-time:
  double              m_fResult;

  bool                m_bClcVarChanged,
                      m_bReady;

  CTracker*           m_pObj;
};

//-------------------------------------------------------------------------------------------------
// CCalcCollection.
//-------------------------------------------------------------------------------------------------
class CCalcCollection : public std::vector<CCalculator*>
{
public:
  virtual ~CCalcCollection();

  void              clear_calcs();
  void              calculate();

  bool              sel_region_changed(CStringVector* pRegNames);

  void              save(CArchive& ar);
  void              load(CArchive& ar);
};

//-------------------------------------------------------------------------------------------------
// CPlaneYZCalculator.
//-------------------------------------------------------------------------------------------------
class CPlaneYZCalculator : public CCalculator
{
public:
  CPlaneYZCalculator();
  virtual ~CPlaneYZCalculator();

  virtual void      run();
  virtual void      do_calculate();
  virtual void      update();
  virtual void      clear();

  virtual int       type() const { return ctPlaneYZ; }

  double            get_plane_x() const;
  DWORD_PTR         get_plane_x_ptr() const;
  void              set_plane_x(double fX);

  double            get_mesh_step() const;
  DWORD_PTR         get_mesh_step_ptr() const;
  void              set_mesh_step(double fStep);

// Sequential calculations support. 
  double            get_start_pos() const;
  DWORD_PTR         get_start_pos_ptr() const;
  void              set_start_pos(double fX);

  double            get_end_pos() const;
  DWORD_PTR         get_end_pos_ptr() const;
  void              set_end_pos(double fX);

  UINT              get_seq_calc_count() const;
  DWORD_PTR         get_seq_calc_count_ptr() const;
  void              set_seq_calc_count(UINT nCount);

  const char*       get_filename() const;
  DWORD_PTR         get_filename_ptr() const;
  void              set_filename(const char* pName);

  void              do_sequence_calc();

//-------------------------------------------------------------------------------------------------
// Streaming:
//-------------------------------------------------------------------------------------------------
  virtual void      save(CArchive& ar);
  virtual void      load(CArchive& ar);

protected:
  void              set_default();

  bool              build_cs_mesh();

  void              triangulate(UINT nY, UINT nZ, double fdY, double fdZ, double fY0, double fZ0, CNode3D** pNodePtrs);

  void              process_quad(const Vector3D& v0,
                                 const Vector3D& v1,
                                 const Vector3D& v2,
                                 const Vector3D& v3,
                                 CNode3D*&       p0, 
                                 CNode3D*&       p1,
                                 CNode3D*&       p2,
                                 CNode3D*&       p3);

  void              process_triangle(const Vector3D& v0,
                                     const Vector3D& v1,
                                     const Vector3D& v2,
                                     CNode3D*&       p0, 
                                     CNode3D*&       p1, 
                                     CNode3D*&       p2);

  bool              move_node(const Vector3D& vOrig, Vector3D& vDest) const;

private:
  double            m_fPosX,        // x-coordinate of the cross-section plane.
                    m_fMeshSize;

  CNodesCollection  m_vNodes;
  CFacesCollection  m_vFaces;

// Sequential calculations support:
  double            m_fStartX,
                    m_fEndX;

  UINT              m_nSeqCalcCount;

  std::string       m_sOutputFile;
};

//-------------------------------------------------------------------------------------------------
// CSelectedRegionCalculator.
//-------------------------------------------------------------------------------------------------
class CSelectedRegionCalculator : public CCalculator
{
public:
  CSelectedRegionCalculator();
  virtual ~CSelectedRegionCalculator();

  virtual void        run();
  virtual void        do_calculate();
  virtual void        update();
  virtual void        clear();

  virtual int         type() const { return ctSelRegions; }

  CString             get_sel_reg_names() const;
  DWORD_PTR           get_sel_reg_names_ptr() const;

//-------------------------------------------------------------------------------------------------
// Streaming:
//-------------------------------------------------------------------------------------------------
  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

protected:
  bool                get_sel_regions(CRegionsCollection& vRegions, UINT& nFaceCount) const;

private:
  CStringVector       m_vSelRegNames;
};

//-------------------------------------------------------------------------------------------------
// CLineCalculator.
//-------------------------------------------------------------------------------------------------
class CLineCalculator : public CCalculator
{
public:
  CLineCalculator();
  virtual ~CLineCalculator();

  enum  // calculated variable.
  {
    lcPress   = 0,
    lcTemp    = 1,
    lcDens    = 2,
    lcVx      = 3,
    lcRe      = 4,
    lcEx      = 5,
    lcEy      = 6,
    lcEz      = 7,
    lcRFx     = 8,
    lcRFy     = 9,
    lcRFz     = 10,
    lcCount   = 11
  };

  virtual void        run();
  virtual void        do_calculate();
  virtual void        update();
  virtual void        clear();

  virtual int         type() const { return ctAlongLine; }

  virtual const char* units() const;
  virtual const char* get_var_name(int nVar) const;

  virtual bool        get_update_flag() const;  // calculators which must re-calculate at every change return true.
  virtual int         calc_vars_count() const;

  const char*         get_filename() const;
  DWORD_PTR           get_filename_ptr() const;
  void                set_filename(const char* pName);

  Vector3D            get_start() const;
  DWORD_PTR           get_start_ptr() const;
  void                set_start(const Vector3D& vPos);

  Vector3D            get_end() const;
  DWORD_PTR           get_end_ptr() const;
  void                set_end(const Vector3D& vPos);

  UINT                get_steps_count() const;
  DWORD_PTR           get_steps_count_ptr() const;
  void                set_steps_count(UINT nCount);

//-------------------------------------------------------------------------------------------------
// Streaming:
//-------------------------------------------------------------------------------------------------
  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

protected:
  void                set_default();
  void                assign_result(const CNode3D& node);

private:
  std::string         m_sOutputFile;

  UINT                m_nStepCount;

  Vector3D            m_vLineStart,
                      m_vLineEnd;
};

//-------------------------------------------------------------------------------------------------
// CTrackCalculator - a class for calculating averaged ion parameters by track data.
//-------------------------------------------------------------------------------------------------
class CTrackCalculator : public CCalculator
{
public:
  CTrackCalculator();
  virtual ~CTrackCalculator();

  enum  // calculated variable.
  {
    clcIonTemp    = 0,
    clcCurrent    = 1,
    clcTime       = 2,
    clcTrackCount = 3
  };

  virtual void        run();
  virtual void        do_calculate();
  virtual void        update();
  virtual void        clear();

  virtual int         type() const { return ctTrackCalc; }

  virtual const char* units() const;
  virtual const char* get_var_name(int nVar) const;

  virtual int         calc_vars_count() const;

  const char*         get_filename() const;
  DWORD_PTR           get_filename_ptr() const;
  void                set_filename(const char* pName);

  double              get_cs_pos() const;
  DWORD_PTR           get_cs_pos_ptr() const;
  void                set_cs_pos(double fX);

  double              get_start_x() const;
  DWORD_PTR           get_start_x_ptr() const;
  void                set_start_x(double fX);

  double              get_end_x() const;
  DWORD_PTR           get_end_x_ptr() const;
  void                set_end_x(double fX);

  UINT                get_cs_count() const;
  DWORD_PTR           get_cs_count_ptr() const;
  void                set_cs_count(UINT nCount);

  void                do_sequence_calc();

  virtual void        save(CArchive& ar);
  virtual void        load(CArchive& ar);

// Stationary ion temperature:
  static double       ion_temp(double fGasTemp, const Vector3D& vGasVel, const Vector3D& vIonVel);

protected:
  void                set_default();
  virtual void        convert_result();   // convertion to SI units (or other output units).

// Ion temperature calculation from scratch:
  bool                calc_gas_temp(const Vector3D& vIonPos, double& fGasTemp) const;
  bool                calc_ion_temp(const Vector3D& vIonPos, const Vector3D& vIonVel, double& fIonTemp, double& fGasTemp) const;
  bool                collect_elements(); // collect all the elements, which bounding boxes intersect with x = m_fPosCS plane.

  void                find_reached_last_cs();

private:
  std::string         m_sOutputFile;

  double              m_fPosCS,   // current cross-section position.
                      m_fStartPos,
                      m_fEndPos,
                      
                      m_fIonTempEq, // for the ion temperature the secondary result is the equlibrium ion temperature.
                      m_fGasTemp;

  UINT                m_nCrossSectCount;

// Run-time, for acceleration of the ion temperature calculation from scratch.
  CElementsCollection m_vElements;
// When the travelling time is being measured the averaging must take place over one and the same set of trajectories.
// This vector is of the tracks count size and it contains "true" for those tracks which have reached m_fEndPos. 
  std::vector<bool>   m_vReachedLastCS;
};

class CColoredCrossSection;
//-------------------------------------------------------------------------------------------------
// CTrackCrossSectionCalculator - a class for output ion parameters distributions in a given cross-section.
//-------------------------------------------------------------------------------------------------
class CTrackCrossSectionCalculator : public CCalculator
{
public:
  CTrackCrossSectionCalculator();
  virtual ~CTrackCrossSectionCalculator();

  virtual void          run();
  virtual void          do_calculate();
  virtual void          update();
  virtual void          clear();

  virtual int           type() const { return ctTrackCrossSect; }

  virtual const char*   units() const;
  virtual const char*   get_var_name(int nVar) const;

  const char*           get_filename() const;
  DWORD_PTR             get_filename_ptr() const;
  void                  set_filename(const char* pName);

  virtual void          save(CArchive& ar);
  virtual void          load(CArchive& ar);

  CColoredCrossSection* get_object();

// Cross-section statistics, run-time
  Vector3D              get_center() const;
  Vector3D              get_sigma() const;
  UINT                  get_count() const;

protected:
  void                  default_filename();
  void                  calculate_statistics();

private:
  CColoredCrossSection  m_Object;
  std::string           m_sOutputFile;

// Cross-section statistics, run-time
  Vector3D              m_vCenter,
                        m_vSigma;     // mean square deviations.
};

//-------------------------------------------------------------------------------------------------
// Inline implementation
//-------------------------------------------------------------------------------------------------
inline void CCalculator::invalidate()
{
  m_bReady = false;
}

inline bool CCalculator::get_enable() const
{
  return m_bEnable;
}

inline DWORD_PTR CCalculator::get_enable_ptr() const
{
  return (DWORD_PTR)&m_bEnable;
}

inline void CCalculator::set_enable(bool bEnable)
{
  m_bEnable = bEnable;
}

inline int CCalculator::get_clc_var_type() const
{
  return m_nClcVar;
}

inline DWORD_PTR CCalculator::get_clc_var_type_ptr() const
{
  return (DWORD_PTR)&m_nClcVar;
}

inline void CCalculator::set_clc_var_type(int nClcVar)
{
  if(m_nClcVar != nClcVar)
  {
    m_nClcVar = nClcVar;
    m_bClcVarChanged = true;
  }
}

inline double CCalculator::get_char_length() const
{
  return m_fCharLength;
}

inline DWORD_PTR CCalculator::get_char_length_ptr() const
{
  return (DWORD_PTR)&m_fCharLength;
}

inline void CCalculator::set_char_length(double fLen)
{
  m_fCharLength = fLen;
}

inline const char* CCalculator::get_name() const
{
  return m_sCalcName.c_str();
}

inline double CCalculator::get_result() const
{
  return m_fResult;
}

inline DWORD_PTR CCalculator::get_result_ptr() const
{
  return (DWORD_PTR)&m_fResult;
}

inline int CCalculator::calc_vars_count() const
{
  return CCalculator::clcCount;
}

//-------------------------------------------------------------------------------------------------
// Inline implementation: CPlaneYZCalculator
//-------------------------------------------------------------------------------------------------
inline double CPlaneYZCalculator::get_plane_x() const
{
  return m_fPosX;
}

inline DWORD_PTR CPlaneYZCalculator::get_plane_x_ptr() const
{
  return (DWORD_PTR)&m_fPosX;
}

inline void CPlaneYZCalculator::set_plane_x(double fX)
{
  if(m_fPosX != fX)
  {
    m_fPosX = fX;
    m_bReady = false;
  }
}

inline double CPlaneYZCalculator::get_mesh_step() const
{
  return m_fMeshSize;
}

inline DWORD_PTR CPlaneYZCalculator::get_mesh_step_ptr() const
{
  return (DWORD_PTR)&m_fMeshSize;
}

inline void CPlaneYZCalculator::set_mesh_step(double fStep)
{
  if(m_fMeshSize != fStep)
  {
    m_fMeshSize = fStep;
    m_bReady = false;
  }
}

// Sequential calculations support. 
inline double CPlaneYZCalculator::get_start_pos() const
{
  return m_fStartX;
}

inline DWORD_PTR CPlaneYZCalculator::get_start_pos_ptr() const
{
  return (DWORD_PTR)&m_fStartX;
}

inline void CPlaneYZCalculator::set_start_pos(double fX)
{
  m_fStartX = fX;
}

inline double CPlaneYZCalculator::get_end_pos() const
{
  return m_fEndX;
}

inline DWORD_PTR CPlaneYZCalculator::get_end_pos_ptr() const
{
  return (DWORD_PTR)&m_fEndX;
}

inline void CPlaneYZCalculator::set_end_pos(double fX)
{
  m_fEndX = fX;
}

inline UINT CPlaneYZCalculator::get_seq_calc_count() const
{
  return m_nSeqCalcCount;
}

inline DWORD_PTR CPlaneYZCalculator::get_seq_calc_count_ptr() const
{
  return (DWORD_PTR)&m_nSeqCalcCount;
}

inline void CPlaneYZCalculator::set_seq_calc_count(UINT nCount)
{
  m_nSeqCalcCount = nCount;
}

inline const char* CPlaneYZCalculator::get_filename() const
{
  return m_sOutputFile.c_str();
}

inline DWORD_PTR CPlaneYZCalculator::get_filename_ptr() const
{
  return (DWORD_PTR)&m_sOutputFile;
}

inline void CPlaneYZCalculator::set_filename(const char* pName)
{
  m_sOutputFile = pName;
}

//-------------------------------------------------------------------------------------------------
// Inline implementation: CSelectedRegionCalculator
//-------------------------------------------------------------------------------------------------
inline CString CSelectedRegionCalculator::get_sel_reg_names() const
{
  return CObject::compile_string(m_vSelRegNames);
}

inline DWORD_PTR CSelectedRegionCalculator::get_sel_reg_names_ptr() const
{
  return (DWORD_PTR)&m_vSelRegNames;
}

//-------------------------------------------------------------------------------------------------
// Inline implementation: CLineCalculator
//-------------------------------------------------------------------------------------------------
inline int CLineCalculator::calc_vars_count() const
{
  return CLineCalculator::lcCount;
}

inline const char* CLineCalculator::get_filename() const
{
  return m_sOutputFile.c_str();
}

inline DWORD_PTR CLineCalculator::get_filename_ptr() const
{
  return (DWORD_PTR)&m_sOutputFile;
}

inline void CLineCalculator::set_filename(const char* pName)
{
  m_sOutputFile = pName;
}

inline Vector3D CLineCalculator::get_start() const
{
  return m_vLineStart;
}

inline DWORD_PTR CLineCalculator::get_start_ptr() const
{
  return (DWORD_PTR)&m_vLineStart;
}

inline void CLineCalculator::set_start(const Vector3D& vPos)
{
  m_vLineStart = vPos;
}

inline Vector3D CLineCalculator::get_end() const
{
  return m_vLineEnd;
}

inline DWORD_PTR CLineCalculator::get_end_ptr() const
{
  return (DWORD_PTR)&m_vLineEnd;
}

inline void CLineCalculator::set_end(const Vector3D& vPos)
{
  m_vLineEnd = vPos;
}

inline UINT CLineCalculator::get_steps_count() const
{
  return m_nStepCount;
}

inline DWORD_PTR CLineCalculator::get_steps_count_ptr() const
{
  return (DWORD_PTR)&m_nStepCount;
}

inline void CLineCalculator::set_steps_count(UINT nCount)
{
  m_nStepCount = nCount;
}


//-------------------------------------------------------------------------------------------------
// Inline implementation: CTrackCalculator
//-------------------------------------------------------------------------------------------------
inline int CTrackCalculator::calc_vars_count() const
{
  return CTrackCalculator::clcTrackCount;
}

inline const char* CTrackCalculator::get_filename() const
{
  return m_sOutputFile.c_str();
}

inline DWORD_PTR CTrackCalculator::get_filename_ptr() const
{
  return (DWORD_PTR)&m_sOutputFile;
}

inline void CTrackCalculator::set_filename(const char* pName)
{
  m_sOutputFile = pName;
}

inline double CTrackCalculator::get_cs_pos() const
{
  return m_fPosCS;
}

inline DWORD_PTR CTrackCalculator::get_cs_pos_ptr() const
{
  return (DWORD_PTR)&m_fPosCS;
}

inline void CTrackCalculator::set_cs_pos(double fX)
{
  m_fPosCS = fX;
}

inline double CTrackCalculator::get_start_x() const
{
  return m_fStartPos;
}

inline DWORD_PTR CTrackCalculator::get_start_x_ptr() const
{
  return (DWORD_PTR)&m_fStartPos;
}

inline void CTrackCalculator::set_start_x(double fX)
{
  m_fStartPos = fX;
}

inline double CTrackCalculator::get_end_x() const
{
  return m_fEndPos;
}

inline DWORD_PTR CTrackCalculator::get_end_x_ptr() const
{
  return (DWORD_PTR)&m_fEndPos;
}

inline void CTrackCalculator::set_end_x(double fX)
{
  m_fEndPos = fX;
}

inline UINT CTrackCalculator::get_cs_count() const
{
  return m_nCrossSectCount;
}

inline DWORD_PTR CTrackCalculator::get_cs_count_ptr() const
{
  return (DWORD_PTR)&m_nCrossSectCount;
}

inline void CTrackCalculator::set_cs_count(UINT nCount)
{
  m_nCrossSectCount = nCount;
}

//-------------------------------------------------------------------------------------------------
// Inline implementation: CTrackCrossSectionCalculator
//-------------------------------------------------------------------------------------------------
inline const char* CTrackCrossSectionCalculator::get_filename() const
{
  return m_sOutputFile.c_str();
}

inline DWORD_PTR CTrackCrossSectionCalculator::get_filename_ptr() const
{
  return (DWORD_PTR)&m_sOutputFile;
}

inline void CTrackCrossSectionCalculator::set_filename(const char* pName)
{
  m_sOutputFile = pName;
}

inline CColoredCrossSection* CTrackCrossSectionCalculator::get_object()
{
  return &m_Object;
}

// Cross-section statistics, run-time
inline Vector3D CTrackCrossSectionCalculator::get_center() const
{
  return m_vCenter;
}

inline Vector3D CTrackCrossSectionCalculator::get_sigma() const
{
  return m_vSigma;
}

};  // namespace EvaporatingParticle
