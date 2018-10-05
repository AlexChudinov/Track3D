
#pragma once

#include <GL/gl.h>
#include <gl/glu.h>

#include "ScreenCapture.h"
#include "ColorContour.h"

#include "vector2d.hpp"
#include "matrix3d.hpp"
#include <string>
#include <vector>

namespace EvaporatingParticle
{

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
struct CFaceVertex
{
  CFaceVertex()
    : x(0), y(0), z(0), nx(1), ny(0), nz(0)
  {
  }

  CFaceVertex(const Vector3D& pos, const Vector3D& norm)
    : x(pos.x), y(pos.y), z(pos.z), nx(norm.x), ny(norm.y), nz(norm.z)
  {
  }

  double x, y, z;
  double nx, ny, nz;
};

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
struct CEdgeVertex
{
  CEdgeVertex()
    : x(0), y(0), z(0)
  {
  }

  CEdgeVertex(const Vector3D& pos)
    : x(pos.x), y(pos.y), z(pos.z)
  {
  }

  double x, y, z;
};

typedef std::vector<CEdgeVertex> CEdgeVertexColl;

struct CRay;
struct CFace;
struct CNode3D;
struct CRegion;
class  CTracker;
class  CColorContour;

typedef std::vector<std::string> CNamesVector;
typedef std::vector<CColorContour*> CContourColl;
typedef std::vector<CFace*> CExternalFaces;

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
class CTrackDraw
{
public:
  CTrackDraw();
  ~CTrackDraw();

  void              draw();

  enum  // drawing mode.
  {
    dmNone        = 0,
    dmWire        = 1,
    dmFlatAndWire = 2,
    dmFlatOnly    = 3,
    dmCount       = 4
  };

  CTracker*         get_tracker() const;
  void              set_tracker(CTracker* pTr);

  int               get_draw_mode() const;
  DWORD_PTR         get_draw_mode_ptr() const;
  void              set_draw_mode(int nMode);

  const char*       get_draw_mode_name(int nMode) const;

  bool              get_enable_tracks() const;
  DWORD_PTR         get_enable_tracks_ptr() const;
  void              set_enable_tracks(bool bEnable);

  double            get_opacity() const;
  DWORD_PTR         get_opacity_ptr() const;
  void              set_opacity(double fAlpha);

  COLORREF          get_faces_color() const;
  DWORD_PTR         get_faces_color_ptr() const;
  void              set_faces_color(COLORREF clr);

  COLORREF          get_bkgr_color() const;
  DWORD_PTR         get_bkgr_color_ptr() const;
  void              set_bkgr_color(COLORREF clr);

  CNamesVector      get_hidden_reg_names() const;
  DWORD_PTR         get_hidden_reg_names_ptr() const;
  void              set_hidden_reg_names();

  bool              get_rot_center() const;
  DWORD_PTR         get_rot_center_ptr() const;
  void              set_rot_center(bool bEnable);

// Auxiliary lines (visualization of normals or Dirichlet cells)
  bool              get_enable_draw_norm() const;
  DWORD_PTR         get_enable_draw_norm_ptr() const;
  void              set_enable_draw_norm(bool bEnable);

  UINT              get_cell_index() const;
  DWORD_PTR         get_cell_index_ptr() const;
  void              set_cell_index(UINT nId);

// Contours and vector plots:
  void              add_contour();
  void              remove_contour(CColorContour* pItem);
  void              invalidate_contour(DWORD_PTR pRegNames);

  size_t            get_contours_count() const;
  CColorContour*    get_contour(size_t nIndex) const;

// Windows kitchen:
  HWND              get_window_handle();
  void              set_window_handle(HWND hwnd);

  bool              create_window_layer();

// Global settings:
  void              set_global();

// Coloring:
  void              set_lights();
  void              set_materials();

  void              invalidate_all();       // makes the draw object re-build vertex and normal arrays.
  void              invalidate_tracks();
  void              invalidate_geometry();
  void              invalidate_contours();
  void              invalidate_hidden();
  void              new_data();             // makes the draw object set the default projection.

// Move and Rotate toolbar support:
  enum
  {
    nContextMove    = 0,
    nContextRotX    = 1,
    nContextRotY    = 2,
    nContextRotZ    = 3,
  };

  int               get_context() const;
  void              set_context(int nContext);

// Mouse events support:
  void              on_mouse_move(const CPoint& point);
  void              on_mouse_wheel(short nDelta, const CPoint& point);

  void              on_left_button_down(const CPoint& point);
  void              on_left_button_up(const CPoint& point);

  void              on_right_button_up(const CPoint& point);

// Selecting regions support:
  bool              get_sel_flag() const;
  void              enter_sel_context(CNamesVector* pRegNames, bool bAllowSel = true);
  void              exit_sel_context(CNamesVector* pRegNames);
  void              hide_selected();

  bool              get_ctrl_pressed() const;
  void              set_ctrl_pressed(bool bCtrlKeyDown);

  CRegion*          get_region_under_cursor() const;
  void              set_region_under_cursor(CRegion* pReg);

  CString           get_hidden_names_str() const;  // merges the names in one single string using ", " as delimiters.

  void              show_all_regions(); // this is the only way to show hidden regions.
  void              set_visibility_status(CNamesVector* pRegNames, bool bVisible);

// Cross-sections of the calculators suppport:
  void              set_cross_sections_array(CExternalFaces* pFaces); // called from Calculator::update().

// Saving image support:
  bool              capture_image();
  bool              save_image(const CString& cFileName);
  bool              is_busy();

// Tracks color map:
  CColoredTracks&   get_colored_tracks();

// Streaming:
  void              save(CArchive& archive);
  void              load(CArchive& archive);

  void              clear();
  void              build_wireframe_array();  // public because it can be called from CTracker::load().

protected:
  void              set_projection();
  void              set_view(int nViewDir = dirNone);

  void              get_resolution(long& nx, long& ny) const;

  void              build_arrays();
  void              build_faces_array();
  void              build_norm_array();
  void              build_aux_array();

  void              draw_tracks();
  void              draw_geometry();
  void              draw_contours();
  void              draw_axes();

  void              draw_flat();
  void              draw_wire();
  void              draw_norm();

  void              draw_selected_faces();

// These two functions do wireframe drawing:
  void              draw_sel_region();
  void              draw_cross_sections();

  const char*       window_caption() const;

// Cursor coordinates in the status line:
  void              screen_to_world(const CPoint& point, Vector3D& world, bool bWorldDepth = true) const;
// Drawing progress in the status line:
  void              set_progress(const char* cJobName, int nPercent) const;

  CRay              get_view_dir(const CPoint& point) const;
  CRegion*          intersect(const CRay& ray) const;

  std::string       dbl_to_str(double val) const; // keeps three digits after decimal point.

// Text output support:
  void              get_inv_rot(double* pInvRotMtx);  // double pInvRotMtx[16] is assumed to be declared in the calling function.

// Check whether the mouse cursor is over some axis sign (X, -X, Y, -Y, Z or -Z).
  void              get_over_axis(const CPoint& point);
  void              draw_letter(double* pInvRot, const Vector3D& vWorldPos, int nLetterType);

private:
  std::vector<CFaceVertex> m_vFaceVert;
  std::vector<CFaceVertex> m_vSelFaceVert;    // selected faces (regions).

// Colored Tracks
  CColoredTracks    m_ColoredTracks;

// Contours and vector plots:
  CContourColl      m_vContours;

// Highlited region drawing:
  CRegion*          m_pRegUnderCursor;
  CEdgeVertexColl   m_vSelRegVert; // the mesh lines of the highlighted region will be drawn.

// Calculator cross-section(s) drawing:
  CEdgeVertexColl   m_vCrossSectVert;
// Wireframe lines:
  CEdgeVertexColl   m_vWireFrame;
// Auxiliary lines (visualization of normals or Dirichlet cells):
  CEdgeVertexColl   m_vAuxLines;

  int               m_nDrawMode;  // can be one of the following: dmNone, dmWire and dmFlat.
  bool              m_bDrawTracks,
                    m_bDrawNorm;  // optional auxiliary lines drawing.

  UINT              m_nDrawnCell; // index of the Dirichlet cell to visualize.

// Names of regions to be hidden:
  CNamesVector      m_vHiddenRegNames;

  double            m_fOpacity;
  COLORREF          m_Color;
  COLORREF          m_BkgrColor;
  float             m_fBkgrRed,
                    m_fBkgrGreen,
                    m_fBkgrBlue;

  bool              m_bRotCenter; // if true the rotation takes place around the center m_vCenter, otherwise around zero point.
  Vector3D          m_vCenter;    // run-time variable used generally in tests.

  CTracker*         m_pTracker;

  HWND              m_hWnd;
  HDC               m_hDC;
  HGLRC             m_hRC;

  double            m_fPix2cm,
                    m_fScrWidth,  // in cm.
                    m_fAxisLen,
                    m_fLetDist;

// Mouse events support:
  enum
  {
    nRegimeNone   = 0,
    nRegimeMove   = 1,
    nRegimeRotate = 2,
  };

  int               m_nRegime,
                    m_nContext;

  bool              m_bSelFlag,     // this flag becomes "true" in enter_sel_context(...) and "false" in exit_sel_context(...).
                    m_bCtrlPressed; // true if the Control key is pressed.

  CPoint            m_StartPoint;

  Vector2D          m_vShift;     // shift in the viewport coordinate system, in cm. 

  double            m_fScale,
                    m_fRotAngleX,
                    m_fRotAngleY,
                    m_fRotAngleZ;

// Run-time:
  bool              m_bWireframeReady,
                    m_bFacesReady,
                    m_bNormReady,
                    m_bNewData,
                    m_bBusy;

// Saving image support:
  CScreenImage      m_WndImage;

// Drawing axes support:
  double            m_pAxesTriad[18];
  double            m_pXYZ[54]; // 54 = 48 + 6, for X, Y, Z and for one sign "-" (minus).

  UINT              m_nOvrAxis;

  enum
  {
    dirNone   = 0,
    dirPlusX  = 1,
    dirMinusX = 2,
    dirPlusY  = 3,
    dirMinusY = 4,
    dirPlusZ  = 5,
    dirMinusZ = 6
  };

  friend class CDirichletTesselation;
};


inline const char* CTrackDraw::window_caption() const
{
  return "Tracks of Evaporating Particles";
}

inline HWND CTrackDraw::get_window_handle()
{
  return m_hWnd;
}

inline CTracker* CTrackDraw::get_tracker() const
{
  return m_pTracker;
}

inline int CTrackDraw::get_draw_mode() const
{
  return m_nDrawMode;
}

inline DWORD_PTR CTrackDraw::get_draw_mode_ptr() const
{
  return (DWORD_PTR)&m_nDrawMode;
}

inline void CTrackDraw::set_draw_mode(int nMode)
{
  m_nDrawMode = nMode;
}

inline bool CTrackDraw::get_enable_tracks() const
{
  return m_bDrawTracks;
}

inline DWORD_PTR CTrackDraw::get_enable_tracks_ptr() const
{
  return (DWORD_PTR)&m_bDrawTracks;
}

inline void CTrackDraw::set_enable_tracks(bool bEnable)
{
  m_bDrawTracks = bEnable;
}

inline double CTrackDraw::get_opacity() const
{
  return m_fOpacity;
}

inline DWORD_PTR CTrackDraw::get_opacity_ptr() const
{
  return (DWORD_PTR)&m_fOpacity;
}

inline void CTrackDraw::set_opacity(double fAlpha)
{
  m_fOpacity = fAlpha;
}

inline COLORREF CTrackDraw::get_faces_color() const
{
  return m_Color;
}

inline DWORD_PTR CTrackDraw::get_faces_color_ptr() const
{
  return (DWORD_PTR)&m_Color;
}

inline void CTrackDraw::set_faces_color(COLORREF clr)
{
  m_Color = clr;
}

inline COLORREF CTrackDraw::get_bkgr_color() const
{
  return m_BkgrColor;
}

inline DWORD_PTR CTrackDraw::get_bkgr_color_ptr() const
{
  return (DWORD_PTR)&m_BkgrColor;
}

inline void CTrackDraw::set_bkgr_color(COLORREF clr)
{
  m_BkgrColor = clr;
  float fc = 1.0f / 255;
  m_fBkgrRed = fc * GetRValue(clr);
  m_fBkgrGreen = fc * GetGValue(clr);
  m_fBkgrBlue = fc * GetBValue(clr);
}

inline CNamesVector CTrackDraw::get_hidden_reg_names() const
{
  return m_vHiddenRegNames;
}

inline DWORD_PTR CTrackDraw::get_hidden_reg_names_ptr() const
{
  return (DWORD_PTR)&m_vHiddenRegNames;
}

inline int CTrackDraw::get_context() const
{
  return m_nContext;
}

inline void CTrackDraw::set_context(int nContext)
{
  m_nContext = nContext;
}

inline void CTrackDraw::new_data()
{
  m_bNewData = true;
}

inline CRegion* CTrackDraw::get_region_under_cursor() const
{
  return m_pRegUnderCursor;
}

inline bool CTrackDraw::get_ctrl_pressed() const
{
  return m_bCtrlPressed;
}

inline void CTrackDraw::set_ctrl_pressed(bool bCtrlKeyDown)
{
  m_bCtrlPressed = bCtrlKeyDown;
}

inline bool CTrackDraw::get_rot_center() const
{
  return m_bRotCenter;
}

inline DWORD_PTR CTrackDraw::get_rot_center_ptr() const
{
  return (DWORD_PTR)&m_bRotCenter;
}

inline void CTrackDraw::set_rot_center(bool bEnable)
{
  m_bRotCenter = bEnable;
}

inline bool CTrackDraw::get_enable_draw_norm() const
{
  return m_bDrawNorm;
}

inline DWORD_PTR CTrackDraw::get_enable_draw_norm_ptr() const
{
  return (DWORD_PTR)&m_bDrawNorm;
}

inline void CTrackDraw::set_enable_draw_norm(bool bEnable)
{
  if(m_bDrawNorm != bEnable)
  {
    m_bDrawNorm = bEnable;
    m_bNormReady = false;
  }
}

inline UINT CTrackDraw::get_cell_index() const
{
  return m_nDrawnCell;
}

inline DWORD_PTR CTrackDraw::get_cell_index_ptr() const
{
  return (DWORD_PTR)&m_nDrawnCell;
}

inline void CTrackDraw::set_cell_index(UINT nId)
{
  if(m_nDrawnCell != nId)
  {
    m_nDrawnCell = nId;
    m_bNormReady = false;
  }
}

inline bool CTrackDraw::get_sel_flag() const
{
  return m_bSelFlag;
}

inline size_t CTrackDraw::get_contours_count() const
{
  return m_vContours.size();
}

inline CColorContour* CTrackDraw::get_contour(size_t nIndex) const
{
  return nIndex < m_vContours.size() ? m_vContours.at(nIndex) : NULL;
}

inline CColoredTracks& CTrackDraw::get_colored_tracks()
{
  return m_ColoredTracks;
}

inline void CTrackDraw::invalidate_all()
{
  invalidate_tracks();
  invalidate_geometry();
  invalidate_contours();
  invalidate_hidden();
}

inline void CTrackDraw::invalidate_tracks()
{
  m_ColoredTracks.invalidate();
}

inline void CTrackDraw::invalidate_geometry()
{
  m_bWireframeReady = false;
}

inline void CTrackDraw::invalidate_hidden()
{
  m_bFacesReady = false;
}

inline void CTrackDraw::invalidate_contours()
{
  size_t nCount = m_vContours.size();
  for(size_t i = 0; i < nCount; i++)
    m_vContours.at(i)->invalidate();
}

inline bool CTrackDraw::is_busy()
{
  return m_bBusy;
}

};  // namespace EvaporatingParticle