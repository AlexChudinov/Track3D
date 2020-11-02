
#pragma once

#include "ColorImage.h"

namespace EvaporatingParticle
{

typedef std::vector<Vector3D> CVert3DColl;
//---------------------------------------------------------------------------------------
// CColorContour.
//---------------------------------------------------------------------------------------
class CColorContour : public CColorImage
{
public:
  CColorContour();
  virtual ~CColorContour();

  enum  // index of the variable to be plotted.
  {
    varPress    = 0,
    varDens     = 1,
    varTemp     = 2,
    varAbsVel   = 3,
    varVelX     = 4,
    varVelY     = 5,
    varVelZ     = 6,
    varAbsClmb  = 7,
    varClmbX    = 8,
    varClmbY    = 9,
    varClmbZ    = 10,
    varPhi      = 11,  // sum of potentials of all computed fields, which visualization is enabled.
    varCount    = 12
  };

  int                 get_var_index() const;
  DWORD_PTR           get_var_index_ptr() const;
  void                set_var_index(int nIndex);

  bool                get_enable_lines() const;
  DWORD_PTR           get_enable_lines_ptr() const;
  void                set_enable_lines(bool bEnable);

  double              get_min_val() const;
  DWORD_PTR           get_min_val_ptr() const;
  void                set_min_val(double fVal);

  double              get_max_val() const;
  DWORD_PTR           get_max_val_ptr() const;
  void                set_max_val(double fVal);

  virtual void        draw();

  virtual void        save(CArchive& archive);
  virtual void        load(CArchive& archive);

  virtual void        restore_user_range();

  virtual void        get_min_max();      // minimal and maximal values of the drawn variable over the cross-section...
  virtual void        get_glob_min_max(); // ... and over all the computational domain.

  virtual void        invalidate();

  const char*         get_var_name(int nVar) const;

protected:
  virtual void        set_default();

  virtual bool        build_draw_array();

  void                set_correct_phi_to_nodes();

  void                reorder_vertices(CFace* pFace, Vector3D* pFaceVert, double* pVal) const;

  void                get_face_values_array(CFace* pFace, double* pVal) const;

  double              get_node_value(const CNode3D& node) const;

  void                add_face(const Vector3D& v0, const Vector3D& v1, const Vector3D& v2, const RGB_Color& clr);
  void                add_edge(const Vector3D& v0, const Vector3D& v1);

  double              get_multiplier(bool bSI_to_CGS) const;
  
private:
  int                 m_nVar;
  bool                m_bDrawLines;

  double              m_vUserDefMin[varCount];
  double              m_vUserDefMax[varCount];

  CVert3DColl         m_vFaceVert;
  CColorVector        m_vFaceColors;

  CVert3DColl         m_vIsolines;
};


struct CTrackItem;
typedef std::vector<CVert3DColl> CTracksVert;
typedef std::vector<CColorVector> CTracksColors;
//---------------------------------------------------------------------------------------
// CColoredTracks.
//---------------------------------------------------------------------------------------
class CColoredTracks : public CColorImage
{
public:
  CColoredTracks();
  virtual ~CColoredTracks();

  int                 get_var_index() const;
  DWORD_PTR           get_var_index_ptr() const;
  void                set_var_index(int nIndex);

  double              get_min_val() const;
  DWORD_PTR           get_min_val_ptr() const;
  void                set_min_val(double fVal);

  double              get_max_val() const;
  DWORD_PTR           get_max_val_ptr() const;
  void                set_max_val(double fVal);

  enum  // index of variable the track will be colored by.
  {
    tcmDefault        = 0,
    tcmIonTemp        = 1,
    tcmIonEqTemp      = 2,
    tcmIonConc        = 3, // AC 12/08/2016 Show ion fragmentation percentage
    tcmStartRadiusXY  = 4, // To see how far from the central point at the starting plane the track has started.
    tcmStartRadiusYZ  = 5,
    tcmStartRadiusXZ  = 6,
    tcmStartX         = 7,
    tcmStartY         = 8,
    tcmStartZ         = 9,
    tcmCount          = 10
  };

  virtual void        draw();

  virtual void        save(CArchive& archive);
  virtual void        load(CArchive& archive);

  virtual void        restore_user_range();
  virtual void        get_min_max();

  virtual bool        build_draw_array();

  const char*         get_var_name(int nVar) const;

  static double       get_start_radius(size_t nTrackIndex, int nVarIndex = tcmDefault);

protected:
  virtual void        set_default();

  double              get_vert_value(const CTrackItem& item, size_t nTrackIndex) const;

  double              get_start_coord(size_t nTrackIndex) const;

  void                clear();

private:
  int                 m_nVarColorMap;

  CTracksVert         m_vTracksVert;
  CTracksColors       m_vTracksColors;

  double              m_vUserDefMin[tcmCount];
  double              m_vUserDefMax[tcmCount];
};

//---------------------------------------------------------------------------------------
// CColoredCrossSection - a class for output of the ion beam cross-section colored by 
//                        initial radius. The data could be taken and plotted by Origin.
//---------------------------------------------------------------------------------------
class CColoredCrossSection : public CColorImage
{
public:
  CColoredCrossSection();
  virtual ~CColoredCrossSection();

  enum
  {
    varIonTemp     = 0,
    varStartRadius = 1,
    varRadialField = 2,
    varCount       = 3
  };

  enum  // direction of the cross-section normal, X by default.
  {
    nrmDirX  = 0,
    nrmDirY  = 1,
    nrmDirZ  = 2
  };

  double              get_cross_sect_pos() const;
  DWORD_PTR           get_cross_sect_pos_ptr() const;
  void                set_cross_sect_pos(double fX);

  int                 get_var() const;
  DWORD_PTR           get_var_ptr() const;
  void                set_var(int nVar);

  double              get_min_val() const;
  DWORD_PTR           get_min_val_ptr() const;
  void                set_min_val(double fVal);

  double              get_max_val() const;
  DWORD_PTR           get_max_val_ptr() const;
  void                set_max_val(double fVal);

  void                set_norm_dir(int nDir);

  size_t              get_points_count() const;
  void                get_colored_point(size_t nIndex, Vector3D& vPoint, double& fVal, RGB_Color& clr) const;

  virtual void        draw();
  virtual void        restore_user_range();
  virtual void        get_min_max();

  const char*         get_var_name(int nVar) const;

  virtual void        save(CArchive& archive);
  virtual void        load(CArchive& archive);

protected:
  virtual void        set_default();
  virtual bool        build_draw_array();

  double              get_vert_value(const CTrackItem& item, size_t nTrackIndex) const;

  void                collect_intersections();
  bool                is_intersected(const Vector3D& vPrev, const Vector3D& vCurr);
  void                clear();

private:
  double              m_fCrossSectX;  // can be X, Y or Z depending on m_nDir.

  CVert3DColl         m_vPoints;
  CValueVector        m_vCrossSectValues;
  CColorVector        m_vCrossSectColors;

  int                 m_nVar,
                      m_nDir;

  double              m_vUserDefMin[varCount];
  double              m_vUserDefMax[varCount];
};

//---------------------------------------------------------------------------------------
// CColorContour - inline implementation.
//---------------------------------------------------------------------------------------
inline int CColorContour::get_var_index() const
{
  return m_nVar;
}

inline DWORD_PTR CColorContour::get_var_index_ptr() const
{
  return (DWORD_PTR)&m_nVar;
}

inline void CColorContour::set_var_index(int nIndex)
{
  if(m_nVar != nIndex)
  {
    m_nVar = nIndex;
    if(m_bUserDefRange && (m_nVar < varCount))
      restore_user_range();

    m_bReady = false;
  }
}

inline bool CColorContour::get_enable_lines() const
{
  return m_bDrawLines;
}

inline DWORD_PTR CColorContour::get_enable_lines_ptr() const
{
  return (DWORD_PTR)&m_bDrawLines;
}

inline void CColorContour::set_enable_lines(bool bEnable)
{
  if(m_bDrawLines != bEnable)
  {
    m_bDrawLines = bEnable;
    m_bReady = false;
  }
}

inline DWORD_PTR CColorContour::get_min_val_ptr() const
{
  return (DWORD_PTR)&m_fMinVal;
}

inline DWORD_PTR CColorContour::get_max_val_ptr() const
{
  return (DWORD_PTR)&m_fMaxVal;
}

//---------------------------------------------------------------------------------------
// CColoredTracks - inline implementation.
//---------------------------------------------------------------------------------------
inline int CColoredTracks::get_var_index() const
{
  return m_nVarColorMap;
}

inline DWORD_PTR CColoredTracks::get_var_index_ptr() const
{
  return (DWORD_PTR)&m_nVarColorMap;
}

inline void CColoredTracks::set_var_index(int nVar)
{
  if(m_nVarColorMap != nVar)
  {
    m_nVarColorMap = nVar;
    m_bReady = false;
  }
}

inline double CColoredTracks::get_min_val() const
{
  return m_fMinVal;
}

inline DWORD_PTR CColoredTracks::get_min_val_ptr() const
{
  return (DWORD_PTR)&m_fMinVal;
}

inline void CColoredTracks::set_min_val(double fVal)
{
  if(m_fMinVal != fVal)
  {
    m_fMinVal = fVal;
    m_bReady = false;
    if(m_bUserDefRange)
      m_vUserDefMin[m_nVarColorMap] = m_fMinVal;
  }
}

inline double CColoredTracks::get_max_val() const
{
  return m_fMaxVal;
}

inline DWORD_PTR CColoredTracks::get_max_val_ptr() const
{
  return (DWORD_PTR)&m_fMaxVal;
}

inline void CColoredTracks::set_max_val(double fVal)
{
  if(m_fMaxVal != fVal)
  {
    m_fMaxVal = fVal;
    m_bReady = false;
    if(m_bUserDefRange)
      m_vUserDefMax[m_nVarColorMap] = m_fMaxVal;
  }
}

//---------------------------------------------------------------------------------------
// CColoredCrossSection - inline implementation.
//---------------------------------------------------------------------------------------
inline double CColoredCrossSection::get_cross_sect_pos() const
{
  return m_fCrossSectX;
}

inline DWORD_PTR CColoredCrossSection::get_cross_sect_pos_ptr() const
{
  return (DWORD_PTR)&m_fCrossSectX;
}

inline void CColoredCrossSection::set_cross_sect_pos(double fX)
{
  if(m_fCrossSectX != fX)
  {
    m_fCrossSectX = fX;
    m_bReady = false;
  }
}

inline int CColoredCrossSection::get_var() const
{
  return m_nVar;
}

inline DWORD_PTR CColoredCrossSection::get_var_ptr() const
{
  return (DWORD_PTR)&m_nVar;
}

inline void CColoredCrossSection::set_var(int nVar)
{
  if(m_nVar != nVar)
  {
    m_nVar = nVar;
    m_bReady = false;
  }
}

inline double CColoredCrossSection::get_min_val() const
{
  return m_fMinVal;
}

inline DWORD_PTR CColoredCrossSection::get_min_val_ptr() const
{
  return (DWORD_PTR)&m_fMinVal;
}

inline void CColoredCrossSection::set_min_val(double fVal)
{
  if(m_fMinVal != fVal)
  {
    m_fMinVal = fVal;
    m_bReady = false;
    if(m_bUserDefRange)
      m_vUserDefMin[m_nVar] = m_fMinVal;
  }
}

inline double CColoredCrossSection::get_max_val() const
{
  return m_fMaxVal;
}

inline DWORD_PTR CColoredCrossSection::get_max_val_ptr() const
{
  return (DWORD_PTR)&m_fMaxVal;
}

inline void CColoredCrossSection::set_max_val(double fVal)
{
  if(m_fMaxVal != fVal)
  {
    m_fMaxVal = fVal;
    m_bReady = false;
    if(m_bUserDefRange)
      m_vUserDefMax[m_nVar] = m_fMaxVal;
  }
}

inline size_t CColoredCrossSection::get_points_count() const
{
  return m_vPoints.size();
}

inline void CColoredCrossSection::set_norm_dir(int nDir)
{
  m_nDir = nDir;
}

};  // namespace EvaporatingParticle