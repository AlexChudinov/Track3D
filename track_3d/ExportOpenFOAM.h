#pragma once
#ifndef _EXPORTOPENFOAM_
#define _EXPORTOPENFOAM_

#include "Tracker.hpp"

namespace EvaporatingParticle
{

enum
{
  bcWall  = 0,
  bcInlet = 1,
  bcPatch = 2,
  bcSymm  = 3
};

//---------------------------------------------------------------------------------------
// CBoundaryConditions
//---------------------------------------------------------------------------------------
struct CBoundaryConditions
{
  CBoundaryConditions(double fP = 200., double fT = 300., int nt = bcWall)
    : fPress(fP), fTemp(fT), nType(nt)
  {
  }

  double          fPress,   // Pa.
                  fTemp;    // K.

  int             nType;

  CStringVector   vRegNames;

  bool operator == (const CBoundaryConditions& bc)  { return (fPress == bc.fPress) && (fTemp == bc.fTemp) && (nType == bc.nType); }

  void            save(CArchive& archive);
  void            load(CArchive& archive);
};

typedef std::vector<CBoundaryConditions*> CBoundCondColl;

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
struct COpenFoamFace
{
  COpenFoamFace(size_t n0, size_t n1, size_t n2)
    : nNbrIndex(-1)
  {
    nodes[0] = n0;
    nodes[1] = n1;
    nodes[2] = n2;
    nNodesCount = 3;
  }

  COpenFoamFace(size_t n0, size_t n1, size_t n2, size_t n3)
    : nNbrIndex(-1)
  {
    nodes[0] = n0;
    nodes[1] = n1;
    nodes[2] = n2;
    nodes[3] = n3;
    nNodesCount = 4;
  }

  size_t                nodes[4];
  size_t                nNodesCount;  // can be 3 or 4.

  int                   nOwnerIndex;
  int                   nNbrIndex;

  Vector3D              vFaceCenter;

  bool operator == (const COpenFoamFace& face) const;
};

typedef std::vector<COpenFoamFace*> CExportFaceCollection;

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
struct CExportRegion
{
  CExportRegion(const char* pName)
    : sName(pName)
  {
    sTitle = sName;
  }

  CExportFaceCollection vFaces;
  std::string           sTitle;   // a new wall name of the type of "wall#1"; if the region is not a wall both names coincide.
  std::string           sName;    // an old ANSYS name to identify this region.
};

typedef std::vector<CExportRegion*> CPatchCollection;
typedef std::vector<Vector3D> CArrayVector3D;

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
class CExportOpenFOAM : public CObject
{
public:
  CExportOpenFOAM();
  ~CExportOpenFOAM();

  void                  set_default();
  bool                  do_export();

  const char*           get_path() const;
  DWORD_PTR             get_path_ptr() const;
  void                  set_path(const char* pPath);

  bool                  get_enable_bound_cond() const;
  DWORD_PTR             get_enable_bound_cond_ptr() const;
  void                  set_enable_bound_cond(bool bEnable);

  bool                  get_export_internal() const;
  DWORD_PTR             get_export_internal_ptr() const;
  void                  set_export_internal(bool bEnable);

  const char*           get_boundary_cond_file() const;
  DWORD_PTR             get_boundary_cond_file_ptr() const;
  void                  set_boundary_cond_file(const char* pPath);

  int                   get_boundary_object_symm() const;
  DWORD_PTR             get_boundary_object_symm_ptr() const;
  void                  set_boundary_object_symm(int nSymPlanes);

  double                get_coord_shift() const;
  DWORD_PTR             get_coord_shift_ptr() const;
  void                  set_coord_shift(double fShiftX);

  double                get_def_press() const;
  DWORD_PTR             get_def_press_ptr() const;
  void                  set_def_press(double fPress);

  double                get_def_temp() const;
  DWORD_PTR             get_def_temp_ptr() const;
  void                  set_def_temp(double fTemp);

  double                get_def_vx() const;
  DWORD_PTR             get_def_vx_ptr() const;
  void                  set_def_vx(double fVx);

  enum
  {
    euMeters      = 0,
    euCentimeters = 1,
    euMillimeters = 2
  };

  void                  add_bound_cond(CBoundaryConditions* pBC);
  void                  remove_bound_cond(CBoundaryConditions* pBC);

  size_t                get_bc_count() const;
  CBoundaryConditions*  get_bound_cond(size_t nIndex) const;

//-------------------------------------------------------------------------------------------------
// Streaming:
//-------------------------------------------------------------------------------------------------
  void                  save(CArchive& archive);
  void                  load(CArchive& archive);

protected:
  void                  prepare();
  void                  clear_bound_cond();
  void                  clear();

  bool                  export_vertices();
  bool                  export_faces();

  bool                  export_owners();
  bool                  export_neighbors();

  bool                  export_boundary();

  void                  print_bound_names();

  int                   get_bc_type(CExportRegion* pReg) const;
  CBoundaryConditions*  get_bc(CExportRegion* pReg) const;

  void                  collect_internal_faces();
  bool                  process_face(CElem3D* pElem, COpenFoamFace* pFace);
  void                  correct_normal(CElem3D* pOwnerElem, COpenFoamFace* pFace);

  COpenFoamFace*        identify_boundary_face(COpenFoamFace* pFace);

  int                   find_neighbour(CElem3D* pElem, COpenFoamFace* pFace);
  bool                  node_in_elem(CNode3D* pNode, CElem3D* pElem);

  bool                  read_2D_regions();
  void                  merge_2D_regions(); // merge all regions of the "wall" type going one after one.
  void                  merge(CExportRegion* pDest, CExportRegion* pSrc);

// When all the internal faces are collected, add the boundary faces to the end of the collection.
  void                  add_boundary_faces();

  void                  print_header(FILE* pStream);

// Internal data and boundary conditions.
  void                  backup_element_centers();
  void                  export_internal_and_boundary_data();
  void                  export_internal_data(CTracker* pObj, FILE* pFileT, FILE* pFileU, FILE* pFileN);
  void                  export_boundary_data(CTracker* pObj, FILE* pFileT, FILE* pFileU, FILE* pFileN);
  void                  write_uniform_internal_fields(FILE* pFileT, FILE* pFileU, FILE* pFileN);

  void                  print_header_T(FILE* pOutFile);
  void                  print_header_U(FILE* pOutFile);
  void                  print_header_N(FILE* pOutFile);

  void                  print_region_name(const char* pName, size_t nCount, bool bScalar, FILE* pOutFile);
  void                  print_symm_name(const char* pName, FILE* pOutFile);

  Vector3D              get_face_center(COpenFoamFace* pFace) const;
  Vector3D              get_elem_center(const CElem3D* pElem) const;

private:
  std::string           m_sPath;

  std::string           m_sBoundCondPath;

// Symmetry of the CTracker object, data from which are used as boundary conditions. User must set this value manually.
  int                   m_nSymPlanes;

  bool                  m_bEnableBoundCond,
// If this flag is true, internal data from ANSYS are exported. Otherwise, default values are written in the internal points.
                        m_bExportInternal;

  CNodesVector*         m_pNodes;
  CElementsCollection*  m_pElems;

// When the gas parameters at internal points or boundary conditions on inlet are exporting, the global CTracker object is
// re-initialized by the data from the source file. Data in m_pNodes and m_pElems became unusable. Meanwhile, the 3D locations
// of the centers of the elements are needed. These data are stored in m_vElemCenters array.
  CArrayVector3D        m_vElemCenters;

// Collection of pointers to boundary COpenFoamFace objects.
  CPatchCollection      m_Patches;

  CExportFaceCollection m_Faces;
  COpenFoamFace**       m_pAuxFaces;  // array of 6 * m_pElems.size() dimension for quick search through faces collection.

  int                   m_nBoundFaceCount;
  int                   m_nUnits;

// As a rule only a low-press part of the mesh is exported, so that there must be a shift.
// Note: if the mesh is exported as a whole, m_sBoundCondPath can be empty.
  double                m_fShiftX;

// Default internal and boundary values:
  double                m_fDefPress,
                        m_fDefTemp,
                        m_fDefDen,
                        m_fDefVx;

  CBoundCondColl        m_vBoundConditions;
};

typedef std::vector<Vector3D> CLocVector;
//-------------------------------------------------------------------------------------------------
// CExternalGridExport - export to grid nodes
//-------------------------------------------------------------------------------------------------
class CExternalGridExport : public CObject
{
public:
  CExternalGridExport()
    : CObject(), m_sInputFile(_T("")), m_sOutputFile(_T("")), m_bUseSI(false)
  {
  }

  enum  // output file format type (CSV or DAT).
  {
    fmtCSV    = 0,
    fmtDAT4   = 1,
    fmtDAT3   = 2,
    fmtCount  = 3
  };

  static const char*    get_format_name(int nFmtType);

  CString               get_input_file() const;
  DWORD_PTR             get_input_file_ptr() const;
  void                  set_input_file(const CString& sFile);

  CString               get_output_file() const;
  DWORD_PTR             get_output_file_ptr() const;
  void                  set_output_file(const CString& sFile);

  int                   get_format_type() const;
  DWORD_PTR             get_format_type_ptr() const;
  void                  set_format_type(int nType);

  bool                  get_use_SI_units() const;
  DWORD_PTR             get_use_SI_units_ptr() const;

  bool                  do_export();

protected:
  FILE*                 open_input_file() const;
  FILE*                 open_output_file() const;

  void                  read_locations(FILE* pInFile, CLocVector& vLoc);
  void                  interpolate_and_write_data(FILE* pOutFile, const CLocVector& vLoc);

  void                  convert_to_SI(CNode3D& node) const;
  
private:
  CString               m_sInputFile;
  CString               m_sOutputFile;
  bool                  m_bUseSI;
  int                   m_nFmtType;
};

//-------------------------------------------------------------------------------------------------
// Inline implementation. CExportOpenFOAM:
//-------------------------------------------------------------------------------------------------
inline const char* CExportOpenFOAM::get_path() const
{
  return m_sPath.c_str();
}

inline DWORD_PTR CExportOpenFOAM::get_path_ptr() const
{
  return (DWORD_PTR)&m_sPath;
}

inline void CExportOpenFOAM::set_path(const char* pPath)
{
  m_sPath = pPath;
}

inline bool CExportOpenFOAM::get_enable_bound_cond() const
{
  return m_bEnableBoundCond;
}

inline DWORD_PTR CExportOpenFOAM::get_enable_bound_cond_ptr() const
{
  return (DWORD_PTR)&m_bEnableBoundCond;
}

inline void CExportOpenFOAM::set_enable_bound_cond(bool bEnable)
{
  m_bEnableBoundCond = bEnable;
}

inline bool CExportOpenFOAM::get_export_internal() const
{
  return m_bExportInternal;
}

inline DWORD_PTR CExportOpenFOAM::get_export_internal_ptr() const
{
  return (DWORD_PTR)&m_bExportInternal;
}

inline void CExportOpenFOAM::set_export_internal(bool bEnable)
{
  m_bExportInternal = bEnable;
}

inline const char* CExportOpenFOAM::get_boundary_cond_file() const
{
  return m_sBoundCondPath.c_str();
}

inline DWORD_PTR CExportOpenFOAM::get_boundary_cond_file_ptr() const
{
  return (DWORD_PTR)&m_sBoundCondPath;
}

inline void CExportOpenFOAM::set_boundary_cond_file(const char* pPath)
{
  m_sBoundCondPath = pPath;
}

inline int CExportOpenFOAM::get_boundary_object_symm() const
{
  return m_nSymPlanes;
}

inline DWORD_PTR CExportOpenFOAM::get_boundary_object_symm_ptr() const
{
  return (DWORD_PTR)&m_nSymPlanes;
}

inline void CExportOpenFOAM::set_boundary_object_symm(int nSymPlanes)
{
  m_nSymPlanes = nSymPlanes;
}

inline double CExportOpenFOAM::get_coord_shift() const
{
  return 10 * m_fShiftX;
}

inline DWORD_PTR CExportOpenFOAM::get_coord_shift_ptr() const
{
  return (DWORD_PTR)&m_fShiftX;
}

inline void CExportOpenFOAM::set_coord_shift(double fShiftX)
{
  m_fShiftX = 0.1 * fShiftX;
}

inline double CExportOpenFOAM::get_def_press() const
{
  return m_fDefPress;
}

inline DWORD_PTR CExportOpenFOAM::get_def_press_ptr() const
{
  return (DWORD_PTR)&m_fDefPress;
}

inline void CExportOpenFOAM::set_def_press(double fPress)
{
  m_fDefPress = fPress;
  m_fDefDen = m_fDefTemp > Const_Almost_Zero ? m_fDefPress / (Const_Boltzmann_SI * m_fDefTemp) : 0;
}

inline double CExportOpenFOAM::get_def_temp() const
{
  return m_fDefTemp;
}

inline DWORD_PTR CExportOpenFOAM::get_def_temp_ptr() const
{
  return (DWORD_PTR)&m_fDefTemp;
}

inline void CExportOpenFOAM::set_def_temp(double fTemp)
{
  m_fDefTemp = fTemp;
  m_fDefDen = m_fDefTemp > Const_Almost_Zero ? m_fDefPress / (Const_Boltzmann_SI * m_fDefTemp) : 0;
}

inline double CExportOpenFOAM::get_def_vx() const
{
  return m_fDefVx;
}

inline DWORD_PTR CExportOpenFOAM::get_def_vx_ptr() const
{
  return (DWORD_PTR)&m_fDefVx;
}

inline void CExportOpenFOAM::set_def_vx(double fVx)
{
  m_fDefVx = fVx;
}

inline size_t CExportOpenFOAM::get_bc_count() const
{
  return m_vBoundConditions.size();
}

inline CBoundaryConditions* CExportOpenFOAM::get_bound_cond(size_t nIndex) const
{
  return nIndex < m_vBoundConditions.size() ? m_vBoundConditions.at(nIndex) : NULL;
}

// CExternalGridExport:
inline CString CExternalGridExport::get_input_file() const
{
  return m_sInputFile;
}

inline DWORD_PTR CExternalGridExport::get_input_file_ptr() const
{
  return (DWORD_PTR)&m_sInputFile;
}

inline void CExternalGridExport::set_input_file(const CString& sFile)
{
  m_sInputFile = sFile;
}

inline CString CExternalGridExport::get_output_file() const
{
  return m_sOutputFile;
}

inline DWORD_PTR CExternalGridExport::get_output_file_ptr() const
{
  return (DWORD_PTR)&m_sOutputFile;
}

inline void CExternalGridExport::set_output_file(const CString& sFile)
{
  m_sOutputFile = sFile;
}

inline bool CExternalGridExport::get_use_SI_units() const
{
  return m_bUseSI;
}

inline DWORD_PTR CExternalGridExport::get_use_SI_units_ptr() const
{
  return (DWORD_PTR)&m_bUseSI;
}

inline int CExternalGridExport::get_format_type() const
{
  return m_nFmtType;
}

inline DWORD_PTR CExternalGridExport::get_format_type_ptr() const
{
  return (DWORD_PTR)&m_nFmtType;
}

inline void CExternalGridExport::set_format_type(int nType)
{
  m_nFmtType = nType;
}

};  // EvaporatingParticle

#endif