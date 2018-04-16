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
  CDirichletCell(CNode3D* pNode, bool bBoundNode = false);
  ~CDirichletCell();

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

// There is no range check in the following functions. The size of CDirichletTesselation::m_Cells is exactly m_pMesh::get_nodes().size().
  CDirichletCell*     get_cell(size_t nNodeId) const;

  bool                inner_cell(size_t nNodeId) const;

  double              get_cell_volume(size_t nNodeId) const;

// In the following 6 functions nNodeId is the absolute index of the cell's central node, but nNbrId is the index of the neighbour
// in the node's neighbours collection. 0 <= nNbrId < (Central Node).vNbrNodes.size().

// Returns the square of the nNbrId-th face of the nNodeId-th cell.
  double              get_cell_face_square(size_t nNodeId, size_t nNbrId) const;

// Distance between the central node of the nNodeId-th cell and its nNbrId-th neighbouring node.
  double              get_nodes_dist(size_t nNodeId, size_t nNbrId) const;

// Returns normal vector (normalized) at the nNbrId-th face of the nNodeId-th cell.
  Vector3D            get_norm(size_t nNodeId, size_t nNbrId) const;

// Returns normal vector (not normalized) at the nNbrId-th face of the nNodeId-th cell.
  Vector3D            get_norm_vector(size_t nNodeId, size_t nNbrId) const;

  double              get_coeff(size_t nNodeId, size_t nNbrId) const;

  size_t              get_abs_nbr_node_index(size_t nNodeId, size_t nNbrId) const;

  Vector3D            get_grad(size_t nNodeId, const std::vector<float>& vScalarField) const;

// For visualization of the Dirichlet cells this function may be called from CTrackDraw::build_norm_array().
  void                build_cell_in_node(CNode3D* pNode, const CNodesCollection& vNodes);

  bool                init();

protected:
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
  void                order_vert_in_plane(CVertexColl& poly, const Vector3D& vNorm, Vector3D& vFaceCenter) const;

// Input:  Vectors vA, vB and vNorm originate from one point. Vectors vA and vB are in the plane and vNorm is normal to that plane.
// Output: Angle between vA and vB counted in the direction from vA to vB counterclockwise.
  double              angle_0_360(const Vector3D& vA, const Vector3D& vB, const Vector3D& vNorm) const;

// Find the normal vector in the node and add one more plane to the planes collection of the boundary Dirichlet cell.
  void                additional_plane(CPlanesSet& vPlanes, CNode3D* pNode) const;

  void                cell_visualization(const CVertexColl& poly, const Vector3D& vFaceCenter) const;

private:
  CAnsysMesh*         m_pMesh;
  CDirichletCellColl  m_Cells;
  bool                m_bTest;

  bool                m_bReady; // run-time flag.
};

//-------------------------------------------------------------------------------------------------
//  Inline implementation CDirichletCell.
//-------------------------------------------------------------------------------------------------
inline void CDirichletCell::operator delete(void* ptr, size_t n)
{
  ((CDirichletCell*)ptr)->delete_cell(); 
}

//-------------------------------------------------------------------------------------------------
//  CDirichletTesselation.
//-------------------------------------------------------------------------------------------------
inline void CDirichletTesselation::set_mesh(CAnsysMesh* pMesh)
{
  if(pMesh != m_pMesh)
  {
    m_pMesh = pMesh;
    m_bReady = false;
  }
}

inline CDirichletCell* CDirichletTesselation::get_cell(size_t nNodeId) const
{
  return m_Cells.at(nNodeId);
}

inline bool CDirichletTesselation::inner_cell(size_t nNodeId) const
{
  return m_pMesh->get_nodes().at(nNodeId)->vNbrFaces.size() == 0;
}

inline double CDirichletTesselation::get_cell_volume(size_t nNodeId) const
{
  return m_Cells.at(nNodeId)->fVolume;
}

inline double CDirichletTesselation::get_cell_face_square(size_t nNodeId, size_t nNbrId) const
{
  return m_Cells.at(nNodeId)->pFaceSquare[nNbrId];
}

inline double CDirichletTesselation::get_nodes_dist(size_t nNodeId, size_t nNbrId) const
{
  return get_norm_vector(nNodeId, nNbrId).length();
}

inline Vector3D CDirichletTesselation::get_norm(size_t nNodeId, size_t nNbrId) const
{
  return get_norm_vector(nNodeId, nNbrId).normalized();
}

inline size_t CDirichletTesselation::get_abs_nbr_node_index(size_t nNodeId, size_t nNbrId) const
{
  return m_pMesh->get_nodes().at(nNodeId)->vNbrNodes.at(nNbrId);
}

inline double CDirichletTesselation::get_coeff(size_t nNodeId, size_t nNbrId) const
{
  return get_cell_face_square(nNodeId, nNbrId) / (get_nodes_dist(nNodeId, nNbrId) * get_cell_volume(nNodeId));
}

};  // namespace EvaporatingParticle


#endif
