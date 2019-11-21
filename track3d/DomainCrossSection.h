
#pragma once

#include "Elements.h"

namespace EvaporatingParticle
{

struct CTetraEdge
{
  CTetraEdge(size_t i0, size_t i1, UINT excl)
    : n0(i0), n1(i1), nExcl(excl)
  {
  }

  size_t  n0, n1;   // indices of tetra nodes (0, 1, 2, 3) making this edge.
  UINT    nExcl;    // index of the edge, point on which can NOT be a neighbour of the point-intersection with this edge.
};

struct CPyramidEdge
{
  CPyramidEdge(size_t i0, size_t i1, UINT excl1, UINT excl2)
    : n0(i0), n1(i1), nExcl1(excl1), nExcl2(excl2)
  {
  }

  size_t  n0, n1;   // indices of pyramid nodes (0, 1, 2, 3, 4) making this edge.
  UINT    nExcl1,   // indices of the edges, points on which can NOT be a neighbour of the point-intersection with this edge.
          nExcl2;
};

struct CWedgeEdge
{
  CWedgeEdge(size_t i0, size_t i1, UINT excl1, UINT excl2, UINT excl3)
    : n0(i0), n1(i1), nExcl1(excl1), nExcl2(excl2), nExcl3(excl3)
  {
  }

  size_t  n0, n1;   // indices of pyramid nodes (0, 1, 2, 3, 4, 5) making this edge.
  UINT    nExcl1,   // indices of the edges, points on which can NOT be a neighbour of the point-intersection with this edge.
          nExcl2,
          nExcl3;
};

struct CHexaEdge
{
  CHexaEdge(size_t i0, size_t i1, UINT excl1, UINT excl2, UINT excl3, UINT excl4, UINT excl5)
    : n0(i0), n1(i1), nExcl1(excl1), nExcl2(excl2), nExcl3(excl3), nExcl4(excl4), nExcl5(excl5)
  {
  }

  size_t  n0, n1;   // indices of pyramid nodes (0, 1, 2, 3, 4, 5) making this edge.
  UINT    nExcl1,   // indices of the edges, points on which can NOT be a neighbour of the point-intersection with this edge.
          nExcl2,
          nExcl3,
          nExcl4,
          nExcl5;
};

const Vector3D cvDefOrigin(0, 0, 0);
const Vector3D cvDefNorm(0, 0, 1);

class CDomainCrossSection
{
public:
  CDomainCrossSection(const Vector3D& vOrigin = cvDefOrigin, const Vector3D& vNorm = cvDefNorm);
  ~CDomainCrossSection();

  enum  // type of the cross-section plane.
  {
    ptPlaneXY   = 0,
    ptPlaneXZ   = 1,
    ptPlaneYZ   = 2,
    ptCount     = 3
  };

  bool                get_enable() const;
  DWORD_PTR           get_enable_ptr() const;
  bool                set_enable(bool bEnable);

  int                 get_plane_type() const;
  DWORD_PTR           get_plane_type_ptr() const;
  bool                set_plane_type(int nType);

  Vector3D            get_plane_origin() const;
  DWORD_PTR           get_plane_origin_ptr() const;
  bool                set_plane_origin(const Vector3D& vPos);

  Vector3D            get_plane_norm() const;
  DWORD_PTR           get_plane_norm_ptr() const;
  bool                set_plane_norm(const Vector3D& vNorm);

  std::string         get_name() const;
  DWORD_PTR           get_name_ptr() const;
  void                set_name(const std::string& sName);

  void                invalidate();

  const char*         get_type_name(int nPlaneType) const;

  void                save(CArchive& ar);
  void                load(CArchive& ar);

  CRegion*            get_region(); 

protected:
  void                set_default();
  void                clear();

  void                build_mesh();
  bool                intersect_element(CFacesCollection& vFaces, CElem3D* pElem);

  bool                intersect_tetra(CFacesCollection& vFaces, CTetra* pElem);
  bool                intersect_piramid(CFacesCollection& vFaces, CPyramid* pElem);
  bool                intersect_wedge(CFacesCollection& vFaces, CWedge* pElem);
  bool                intersect_hexa(CFacesCollection& vFaces, CHexa* pElem);

// Linear interpolation between two nodes: ksi = 0 at p0, ksi = 1 at p1.
  CNode3D*            interpolate(const Vector3F& vPos, CNode3D* p0, CNode3D* p1, float ksi);

private:
  int                 m_nPlaneType;

// For compatibility with CColorContour the cross-section exists as a CRegion object.
  CRegion             m_Region;
  CNodesCollection    m_vNodes;
  CPlane              m_Plane;

// Run-time:
  bool                m_bReady;
};

//---------------------------------------------------------------------------------------
//  CCrossSectColl
//---------------------------------------------------------------------------------------
class CCrossSectColl : public std::vector<CDomainCrossSection*>
{
public:
  void              save(CArchive& ar);
  void              load(CArchive& ar);
};


inline bool CDomainCrossSection::get_enable() const
{
  return m_Region.bEnabled;
}

inline DWORD_PTR CDomainCrossSection::get_enable_ptr() const
{
  return (DWORD_PTR)&(m_Region.bEnabled);
}

inline bool CDomainCrossSection::set_enable(bool bEnable)
{
  if(m_Region.bEnabled != bEnable)
  {
    m_Region.bEnabled = bEnable;
    return true;
  }

  return false;
}

inline int CDomainCrossSection::get_plane_type() const
{
  return m_nPlaneType;
}

inline DWORD_PTR CDomainCrossSection::get_plane_type_ptr() const
{
  return (DWORD_PTR)&m_nPlaneType;
}

inline Vector3D CDomainCrossSection::get_plane_origin() const
{
  return m_Plane.pos;
}

inline DWORD_PTR CDomainCrossSection::get_plane_origin_ptr() const
{
  return (DWORD_PTR)&(m_Plane.pos);
}

inline bool CDomainCrossSection::set_plane_origin(const Vector3D& vPos)
{
  if((m_Plane.pos - vPos).length() > Const_Almost_Zero)
  {
    m_Plane.pos = vPos;
    m_bReady = false;
    return true;
  }

  return false;
}

inline Vector3D CDomainCrossSection::get_plane_norm() const
{
  return m_Plane.norm;
}

inline DWORD_PTR CDomainCrossSection::get_plane_norm_ptr() const
{
  return (DWORD_PTR)&(m_Plane.norm);
}

inline bool CDomainCrossSection::set_plane_norm(const Vector3D& vNorm)
{
  Vector3D vN = vNorm.normalized();
  if((m_Plane.norm - vN).length() > Const_Almost_Zero)
  {
    m_Plane.norm = vN;
    m_bReady = false;
    return true;
  }

  return false;
}

inline std::string CDomainCrossSection::get_name() const
{
  return m_Region.sName;
}

inline DWORD_PTR CDomainCrossSection::get_name_ptr() const
{
  return (DWORD_PTR)&(m_Region.sName);
}

inline void CDomainCrossSection::set_name(const std::string& sName)
{
  m_Region.sName = sName;
}

inline CRegion* CDomainCrossSection::get_region()
{
  if(!m_bReady)
    build_mesh();

  return m_bReady ? &m_Region : NULL;
}

inline void CDomainCrossSection::invalidate()
{
  m_bReady = false;
}

};  // namespace EvaporatingParticle.