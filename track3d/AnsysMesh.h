#pragma once

#ifndef _AnsysMesh_
#define _AnsysMesh_

#include "CObject.h"
#include "Elements.h"
#include <intsafe.h>

namespace EvaporatingParticle
{

struct CSymCorrData
{
  CSymCorrData()
    : nSymFlag(0), fNrmY(0), fNrmZ(0)
  {
  }

  int     nSymFlag;
  double  fNrmY,
          fNrmZ;
};

class CAnsysMesh : public CObject
{
public:
  CAnsysMesh(bool bAux = false);
  virtual ~CAnsysMesh();

  enum  // Symmetry planes
  {
    spNone  = 0,
    spXY    = 1,
    spXZ    = 2,
    spYZ    = 4,
    spAxial = 8
  };

// User interface
  const char*             get_filename() const;
  DWORD_PTR               get_filename_ptr() const;
  bool                    set_filename(const char* pName);

  int                     get_sym_plane() const;
  DWORD_PTR               get_sym_plane_ptr() const;
  void                    set_sym_plane(int nPlane);

  bool                    get_2D_flag() const;
  DWORD_PTR               get_2D_flag_ptr() const;

  bool                    get_convert_to_cgs() const;
  void                    set_convert_to_cgs(bool bEnable);

  bool                    read_data();  // main function for reading ANSYS data.
  bool                    is_ready();

  bool                    read_geometry();
  bool                    read_2D_regions();

  bool                    read_gasdyn_data(bool bFieldsOnly = false);

  bool                    get_need_read_ansys_field() const;
  void                    set_need_read_ansys_field();

protected:
  void                    add_tetra(CNode3D& p0, CNode3D& p1, CNode3D& p2, CNode3D& p3);
  void                    add_pyramid(CNode3D& p0, CNode3D& p1, CNode3D& p2, CNode3D& p3, CNode3D& p4);
  void                    add_wedge(CNode3D& p0, CNode3D& p1, CNode3D& p2, CNode3D& p3, CNode3D& p4, CNode3D& p5);
  void                    add_hexa(CNode3D& p0, CNode3D& p1, CNode3D& p2, CNode3D& p3, CNode3D& p4, CNode3D& p5, CNode3D& p6, CNode3D& p7);

  void                    bounding_box();
  
  void                    clear();

public:
// Reflect the particle's position, velocity and acceleration against the symmetry plane(s)if necessary.
// The function returns "true" if reflection has been done and the back reflection is needed. In this case set
// bForceReflect to "true" in the next call to ensure that the back reflection will occur.
//  bool                    sym_corr(Vector3D& vPos, Vector3D& vVel, Vector3D& vAccel, bool bForceReflect = false) const;
//  int                     sym_corr(Vector3D& vPos, Vector3D& vVel, Vector3D& vAccel, int nForceReflect = 0) const;

  CSymCorrData            sym_corr_forward(Vector3D& vPos, Vector3D& vVel) const;
  void                    sym_corr_back(Vector3D& vPos, Vector3D& vVel, Vector3D& vAccel, const CSymCorrData& data) const;

  CNodesVector&           get_nodes();
  CElementsCollection&    get_elems();
// If bExtReg == true, the regions created by cross-section planes are added. Otherwise m_vRegions is returned.
  CRegionsCollection&     get_regions(bool bExtReg = true);
  CBox&                   get_box();

  Vector3D                get_center() const;

// Mesh transformation:
  CTransform&             get_transform();

  static CRegion*         get_region(const std::string& sName);     // returns a region pointer or NULL if the name is not found.
  static int              get_region_id(const std::string& sName);  // returns a region index or -1 if the name is not found.

//-------------------------------------------------------------------------------------------------
// Streaming:
//-------------------------------------------------------------------------------------------------
  void                    save(CArchive& archive);
  void                    load(CArchive& archive);

  //AC: static functions to work with mesh
  static const CElem3D * find_global_elem(const CElem3D* elem, const Vector3D& pos);

  static const CElementsCollection& get_global_elements();

protected:
// Returns pointer to the element containing the input point or NULL if no element has been found.
//[AC] const added
	const CElem3D* find_elem(const CElem3D* pPrevElem, const Vector3D& vPos) const;

// Try the nearest neighbors of the element, where the previous location was found, including the element itself.
// Return pointer to the element containing the input point or NULL if no element has been found.
//[AC] const added
  const CElem3D* try_neighbors(const CElem3D* pElem, const Vector3D& vPos) const;

// Convertation to CGS:
  void                    conv_to_cgs(float& fPress, float& fDens, float& fDynVisc, float& fThermCond, float& fCp,
                            float& fVx, float& fVy, float& fVz, float& fEx, float& fEy, float& fEz, float& fRFEx, float& fRFEy, float& fRFEz);
// Calculators:
  void                    invalidate_calculators(); // called from read_data() and makes the calculators update their internal data.

// Termination:
  bool                    abort(FILE* pStream = NULL);

// Particle source management:
  virtual void            invalidate_src() {};

  std::string             m_sDataFile;

  CBox                    m_Box;          // bounding box.

  CNodesVector            m_vNodes;
  CElementsCollection     m_vElems;

  CRegionsCollection      m_vRegions;     // Regions containing triangular faces, for drawing only.
  CRegionsCollection      m_vExtRegions;  // m_vRegions + cross-section regions, run-time.

  bool                    m_bConv2CGS;    // a flag showing whether ANSYS data, which are normally in SI must be 
                                          // converted to the CGS system; this is always "true" so far.

  bool                    m_bMesh2D;      // a flag which is normally false but can be set to true if the mesh is 1 element wide.

  int                     m_nSymPlanes;

// Mesh transformation (applied in read_geometry()).
  CTransform              m_Transform;

// Run-time:
  bool                    m_bReady,       // this flag is set to "false" in set_filename(), to "true" in read_data().
                          m_bAux;         // this flag is set in the constructor and never changed afterwards; the only example - second step of import OpenFOAM.

  bool                    m_bNeedReadAnsysField;  // this flag is set to "true" every time when user checks "Use Ansys Fields" check-box.
};

//---------------------------------------------------------------------------------------
// Inline implementation
//---------------------------------------------------------------------------------------
inline const char* CAnsysMesh::get_filename() const
{
  return m_sDataFile.c_str();
}

inline DWORD_PTR CAnsysMesh::get_filename_ptr() const
{
  return (DWORD_PTR)&m_sDataFile;
}

inline bool CAnsysMesh::set_filename(const char* pName)
{
  if(strcmp(pName, m_sDataFile.c_str()) == 0)
    return false; // filename will not change, no need to set m_bReady flag to false. 

  m_sDataFile = pName;
  m_bReady = false; // force the data files to be re-read.

  invalidate_src(); // make the particle source be re-calculated.
  return true;
}

inline int CAnsysMesh::get_sym_plane() const
{
  return m_nSymPlanes;
}

inline DWORD_PTR CAnsysMesh::get_sym_plane_ptr() const
{
  return (DWORD_PTR)&m_nSymPlanes;
}

inline void CAnsysMesh::set_sym_plane(int nPlane)
{
  m_nSymPlanes = nPlane;
}

inline bool CAnsysMesh::get_2D_flag() const
{
  return m_bMesh2D;
}

inline DWORD_PTR CAnsysMesh::get_2D_flag_ptr() const
{
  return (DWORD_PTR)&m_bMesh2D;
}

inline bool CAnsysMesh::get_convert_to_cgs() const
{
  return m_bConv2CGS;
}

inline void CAnsysMesh::set_convert_to_cgs(bool bEnable)
{
  m_bConv2CGS = bEnable;
}

inline CNodesVector& CAnsysMesh::get_nodes()
{
  return m_vNodes;
}

inline CElementsCollection& CAnsysMesh::get_elems()
{
  return m_vElems;
}

inline CBox& CAnsysMesh::get_box()
{
  return m_Box;
}

inline Vector3D CAnsysMesh::get_center() const
{
  return m_Box.get_center();
}

inline CTransform& CAnsysMesh::get_transform()
{
  return m_Transform;
}

inline bool CAnsysMesh::is_ready()
{
  return m_bReady;
}

inline bool CAnsysMesh::get_need_read_ansys_field() const
{
  return m_bNeedReadAnsysField;
}

inline void CAnsysMesh::set_need_read_ansys_field()
{ 
  m_bNeedReadAnsysField = true;
}


};  // namespace EvaporatingParticle


#endif
