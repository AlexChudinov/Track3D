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
  CDirichletCell();
  ~CDirichletCell();

  UINT        nFaceCount;     // in fact, this is the count of the neighbouring nodes.
  float*      pFaceSquare;    // for inner cells: array of nFaceCount size, contains the squares of polyangle faces, CGS; for boundary cells array of size 10.
  float*      pNbrDist;       // the distance between the central and neighbor nodes; the size of this array is nFaceCount for all cells.
  float       fVolume;        // volume of the inner cell, CGS; zero for boundary cells.

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

// Returns a neighbor vector (not normalized) drawn from the nNodeId-th vertex to its nNbrId-th neighbor.
  Vector3D            get_neighbor_vector(size_t nNodeId, size_t nNbrId) const;

// Returns normal vector at a boundary node. Do not mix with "get_norm" declared above, which merely normalizes the vector drawn from this node to a neighbor node. 
  Vector3D            get_bound_norm(const CNode3D& node) const;

// This condition is used for overdetermined system of equations. We take into account only those nodes, for which (e, vNorm) < Const_Almost_Zero.
  static bool         bound_deriv_calc_cond(const Vector3D& vNorm, const Vector3D& e);

  double              get_coeff(size_t nNodeId, size_t nNbrId) const;

  size_t              get_abs_nbr_node_index(size_t nNodeId, size_t nNbrId) const;

  Vector3D            get_grad(size_t nNodeId, const std::vector<float>& vScalarField) const;

// For visualization of the Dirichlet cells this function may be called from CTrackDraw::build_norm_array().
  CDirichletCell*     build_cell_in_node(const CNode3D& node, const CNodesVector& vNodes);

// The cell is supposed to be a boundary cell.
  Matrix3D            get_bound_cell_mtx(CDirichletCell* pCell) const;  // returns matrix C of the boundary cell.

  bool                init();

  void                invalidate();

protected:
  void                clear();

// For Dirichlet cells built around inner nodes, i.e. for strongly CONVEX cells.
// Input for the following two functions: vPos is the result of intersection of these three planes a, b and c.
  bool                inside(const CPlane& a, const CPlane& b, const CPlane& c, const Vector3D& vPos, const CPlanesSet& vFaces) const;

// Find all vertices of the cell (all possible intersections of the planes forming the cell except for points outside the cell).
// Input:  vInnerPlanes - collection of planes perpendicular to lines connecting this vertex and neighboring vertices;
// Output: vCellVert.
  void                collect_cell_vertices(const CPlanesSet& vInnerPlanes, CVertexColl& vCellVert) const;

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

  void                init_boundary_cell(CDirichletCell* pCell, const CNode3D& node, const CNodesVector& vNodes) const;

  Vector3D            get_bound_cell_grad(CDirichletCell* pCell, const CNode3D& node, const std::vector<float>& vScalarField) const;

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
  return m_pMesh->get_nodes().at(nNodeId).vNbrFaces.size() == 0;
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
  return get_neighbor_vector(nNodeId, nNbrId).length();
}

inline Vector3D CDirichletTesselation::get_norm(size_t nNodeId, size_t nNbrId) const
{
  return get_neighbor_vector(nNodeId, nNbrId).normalized();
}

inline size_t CDirichletTesselation::get_abs_nbr_node_index(size_t nNodeId, size_t nNbrId) const
{
  return m_pMesh->get_nodes().at(nNodeId).vNbrNodes.at(nNbrId);
}

inline double CDirichletTesselation::get_coeff(size_t nNodeId, size_t nNbrId) const
{
  return get_cell_face_square(nNodeId, nNbrId) / (get_nodes_dist(nNodeId, nNbrId) * get_cell_volume(nNodeId));
}

inline void CDirichletTesselation::invalidate()
{
  m_bReady = false;
}

};  // namespace EvaporatingParticle


#endif
