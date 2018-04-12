#pragma once

#ifndef _DIRICHLET_TESSELATION_
#define _DIRICHLET_TESSELATION_

#include "CObject.h"
#include "Elements.h"

namespace EvaporatingParticle
{

//-------------------------------------------------------------------------------------------------
//  CDirichletCell - a finite volume surrounding a mesh node.
//-------------------------------------------------------------------------------------------------
struct CDirichletCell : public BlockAllocator<CDirichletCell>
{
  CDirichletCell(CNode3D* pNode);
  ~CDirichletCell();

  UINT        nCellNode;      // the central node of the cell.
  UINT        nFaceCount;     // the count of the neighbouring nodes.
  float*      pFaceSquare;    // array of nFaceCount size, contains the squares of polyangle faces, CGS.
  float       fVolume;        // volume of the cell, CGS.

  void        delete_cell();

  static void operator delete(void* ptr, size_t n);
};

typedef std::vector<Vector3D> CVertexColl;
typedef std::vector<CDirichletCell*> CDirichletCellColl;

class CAnsysMesh;
//-------------------------------------------------------------------------------------------------
//  CDirichletTesselation - a set of Dirichlet cells built around each inner node of the mesh.
//-------------------------------------------------------------------------------------------------
class CDirichletTesselation : public CObject
{
public:
  CDirichletTesselation(bool bTest = false);
  ~CDirichletTesselation();

  void                set_mesh(CAnsysMesh* pMesh);

// For visualization of the Dirichlet cells this function may be called from CTrackDraw::build_norm_array().
  void                build_cell_in_node(CNode3D* pNode, const CNodesCollection& vNodes);

protected:
  void                init();
  void                clear();

  bool                inside(const Vector3D& vPos, const CPlanesSet& vFaces) const;

// Find all vertices of the cell (all possible intersections of the planes forming the cell except for point outside the cell).
  void                collect_cell_vertices(const CPlanesSet& vPlanes, CVertexColl& vCellVert);

// Selects the vertices belonging to the plane (cell's face) and computes the square of the face.
  double              build_polygon_in_plane(const CPlane& plane, const CVertexColl& vCellVert) const;

// Intersection point (if exists) of three planes.
  bool                three_planes_intersect(const CPlane& a, const CPlane& b, const CPlane& c, Vector3D& vRes) const;

// Input:  unordered collection of vertices belonging to one and the same plane and forming a convex polygon.
// Output: ordered vertices of the polygon so that the polygon's region is always on the left-hand side if we go around.
  void                order_vert_in_plane(CVertexColl& poly, const CPlane& face) const;

// Input:  Vectors vA, vB and vNorm originate from one point. Vectors vA and vB are in the plane and vNorm is normal to that plane.
// Output: Angle between vA and vB counted in the direction from vA to vB counterclockwise.
  double              angle_0_360(const Vector3D& vA, const Vector3D& vB, const Vector3D& vNorm) const;

  void                cell_visualization(const CVertexColl& poly, const Vector3D& v0) const;

private:
  CAnsysMesh*         m_pMesh;
  CDirichletCellColl  m_Cells;
  bool                m_bTest;
};

//-------------------------------------------------------------------------------------------------
//  Inline implementation.
//-------------------------------------------------------------------------------------------------
inline void CDirichletCell::operator delete(void* ptr, size_t n)
{
  ((CDirichletCell*)ptr)->delete_cell(); 
}

inline void CDirichletTesselation::set_mesh(CAnsysMesh* pMesh)
{
  if(pMesh != m_pMesh)
    m_pMesh = pMesh;

  if(m_bTest)
    return;

  init();
}

};  // namespace EvaporatingParticle


#endif
