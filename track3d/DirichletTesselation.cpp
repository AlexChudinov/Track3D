#pragma once

#include "stdafx.h"

#include "math.h"
#include "float.h"
#include "AnsysMesh.h"

#include "DirichletTesselation.h"

#include <map>
#include <algorithm>

#include "ParticleTracking.h" // for the Dirichlet cells visualization.


namespace EvaporatingParticle
{

CDirichletCell::CDirichletCell(CNode3D* pNode, bool bBoundNode)
{
  nFaceCount = bBoundNode ? pNode->vNbrNodes.size() + 1 : pNode->vNbrNodes.size();
  pFaceSquare = new float[nFaceCount];
  for(UINT i = 0; i < nFaceCount; i++)
    pFaceSquare[i] = 0;

  fVolume = 0;
}

CDirichletCell::~CDirichletCell()
{
  delete[] pFaceSquare;
}

void CDirichletCell::delete_cell()
{ 
  BlockPool<CDirichletCell>::freeBlock(this);
}


//-------------------------------------------------------------------------------------------------
//  CDirichletTesselation - a set of Dirichlet cells built around each inner node of the mesh.
//-------------------------------------------------------------------------------------------------
CDirichletTesselation::CDirichletTesselation(bool bTest)
  : m_bTest(bTest), m_pMesh(NULL), m_bReady(false)
{
}

CDirichletTesselation::~CDirichletTesselation()
{
  clear();
}

Vector3D CDirichletTesselation::get_grad(size_t nNodeId, const std::vector<float>& vScalarField) const
{
  CNodesCollection& vNodes = m_pMesh->get_nodes();
  bool bOK = (vNodes.size() == m_Cells.size()) && (vNodes.size() == vScalarField.size()) && (nNodeId < vNodes.size());
  if(!bOK)
    return Vector3D(0, 0, 0);

  CNode3D* pNode = vNodes.at(nNodeId);      // the cell's central node.
  double fPhiC = vScalarField.at(nNodeId);  // the value of the potential at the cell's center.

  Vector3D vRes(0, 0, 0);
  CDirichletCell* pCell = m_Cells.at(nNodeId);
  double fVol = pCell->fVolume, fFaceSqr, fFacePhi;
  Vector3D vFaceNorm;

  size_t nNbrNodeIndex;
  size_t nNbrCount = vNodes.at(nNodeId)->vNbrNodes.size();
  for(size_t i = 0; i < nNbrCount; i++)
  {
    fFaceSqr = pCell->pFaceSquare[i];
    nNbrNodeIndex = pNode->vNbrNodes.at(i);
    fFacePhi = fPhiC + vScalarField.at(nNbrNodeIndex);
    vFaceNorm = get_norm(nNodeId, i);

    vRes += vFaceNorm * (fFacePhi * fFaceSqr);
  }

  vRes /= (2 * fVol);
  return vRes;
}

bool CDirichletTesselation::init()
{
  if(m_pMesh == NULL)
    return false;

  if(m_bReady)
    return true;

  set_job_name("Creating Dirichlet cells...");
  set_progress(0);

  clear();

  CNodesCollection& vNodes = m_pMesh->get_nodes();
  CNode3D* pNode = NULL;
  size_t nNodeCount = vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    pNode = vNodes.at(i);
    build_cell_in_node(pNode, vNodes);

    if(i % 100 == 0)
      set_progress(int(0.5 + 100 * i / nNodeCount));
    if(m_bTerminate)
      return false;
  }

  m_bReady = true;
  return true;
}

void CDirichletTesselation::clear()
{
  size_t nSize = m_Cells.size();
  for(size_t i = 0; i < nSize; i++)
    delete m_Cells.at(i);

  m_Cells.clear();
}

void CDirichletTesselation::build_cell_in_node(CNode3D* pNode, const CNodesCollection& vNodes)
{
  bool bBoundNode = pNode->vNbrFaces.size() > 0;
  CDirichletCell* pCell = new CDirichletCell(pNode, bBoundNode);

// 1. Collect all planes forming the Dirichlet cell around this node.
  size_t nPlanesCount = pCell->nFaceCount;
  CPlanesSet vPlanes;
  vPlanes.reserve(nPlanesCount);

  UINT nId;
  CNode3D* pNbrNode = NULL;
  Vector3D vC = pNode->pos, vOrg, vNorm, vNbr;  // cell center, plane origin and normal and neighbour cell center, respectively.
  size_t nNbrCount = pNode->vNbrNodes.size();   // for inner nodes nPlanesCount = nNbrCount, but for boundary nodes nPlanesCount = nNbrCount + 1.
  for(size_t i = 0; i < nNbrCount; i++)
  {
    nId = pNode->vNbrNodes.at(i);
    pNbrNode = vNodes.at(nId);
    vNbr = pNbrNode->pos;
    vOrg = 0.5 * (vC + vNbr);
    vNorm = (vNbr - vC).normalized();   // plane normal looks always out of the cell.
    CPlane face(vOrg, vNorm);
    vPlanes.push_back(face);
  }

  if(bBoundNode)
    additional_plane(vPlanes, pNode);

// 2. Collect all vertices of the Dirichlet cell as all possible intersections of triads of different planes.
//    Exclude the points lying outside the cell.
  CVertexColl vCellVert;
  collect_cell_vertices(vPlanes, vCellVert);

// 3. For each plane select vertices belonging to the plane and build the cell's polygon. Compute the polygon's square.
  double fS = 0;
  pCell->fVolume = 0;
  nPlanesCount = vPlanes.size();
  for(size_t j = 0; j < nPlanesCount; j++)
  {
    const CPlane& plane = vPlanes.at(j);
    fS = build_polygon_in_plane(plane, vCellVert); // the square of that face of the cell, which belongs to this plane.

    pCell->pFaceSquare[j] = fS;
    pCell->fVolume += Const_One_Third * fS * (vOrg - vC).length();  // V = 1/3 * S * h, a pyramid volume.
  }

  m_Cells.push_back(pCell);
}

void CDirichletTesselation::collect_cell_vertices(const CPlanesSet& vPlanes, CVertexColl& vCellVert)
{
  Vector3D vVert;
  size_t nPlanesCount = vPlanes.size();
  for(size_t i = 0; i < nPlanesCount; i++)
  {
    const CPlane& PlaneA = vPlanes.at(i);
    for(size_t j = 0; j < i; j++)
    {
      if(j == i)
        continue;

      const CPlane& PlaneB = vPlanes.at(j);
      for(size_t k = 0; k < j; k++)
      {
        if(k == j || k == i)
          continue;

        const CPlane& PlaneC = vPlanes.at(k);
        if(!three_planes_intersect(PlaneA, PlaneB, PlaneC, vVert))
          continue;

        if(!inside(vVert, vPlanes))
          continue;

        vCellVert.push_back(vVert);
      }
    }
  }
}

double CDirichletTesselation::build_polygon_in_plane(const CPlane& plane, const CVertexColl& vCellVert) const
{
  Vector3D vVert;
  CVertexColl vPoly;
  double fDist;

  size_t nVertCount = vCellVert.size();
// Find all vertices belonging to the face.
  for(size_t i = 0; i < nVertCount; i++)
  {
    vVert = vCellVert.at(i);
    fDist = fabs((vVert - plane.pos) & plane.norm);
    if(fDist > Const_Almost_Zero)
      continue;   // the vertex is out of the plane.

    vPoly.push_back(vVert);
  }

// At this point vPoly contains unordered vertices of the face polygon. What we need is to order them and find the polygon's square.
  Vector3D vC;
// Note: 1. After this call vPoly contains ordered vertices of the polygone.
//       2. Vector coordinates in vPoly are in the c.s. originating in the face center vC. 
  order_vert_in_plane(vPoly, plane.norm, vC);

  if(m_bTest)
    cell_visualization(vPoly, vC);

// When the vertices in the polygon are ordered, the square can be easily found:
  double fSquare = 0;
  size_t nPolyCount = vPoly.size();
  for(size_t k = 1; k < nPolyCount; k++)
    fSquare += 0.5 * (vPoly.at(k - 1) * vPoly.at(k)).length();

  return fSquare;
}

static const Vector3D scvNull(0, 0, 0);

void CDirichletTesselation::order_vert_in_plane(CVertexColl& poly, const Vector3D& vNorm, Vector3D& vC) const
{
  size_t nVertCount = poly.size();
  if(nVertCount == 0)
    return;

  vC = scvNull;
// For correct ordering and face square calculation the face center must be inside the face.
  for(size_t j = 0; j < nVertCount; j++)
    vC += poly.at(j);

  vC /= (double)nVertCount;

  double fAng = 0;
  Vector3D v0 = poly.at(0) - vC, v1;

  std::map<double, Vector3D> ordered_poly;
  ordered_poly.insert(std::pair<double, Vector3D>(fAng, v0));

  for(size_t i = 1; i < nVertCount; i++)
  {
    v1 = poly.at(i) - vC;
    fAng = angle_0_360(v0, v1, vNorm);
    ordered_poly.insert(std::pair<double, Vector3D>(fAng, v1));
  }

  poly.clear();
  std::map<double, Vector3D>::iterator iter;
  for(iter = ordered_poly.begin(); iter != ordered_poly.end(); iter++)
    poly.push_back(iter->second);
}

bool CDirichletTesselation::three_planes_intersect(const CPlane& faceA, const CPlane& faceB, const CPlane& faceC, Vector3D& vRes) const
{
  Matrix3D m(faceA.norm, faceB.norm, faceC.norm);
  m.transpose();
  double D = m.det();
  if(fabs(D) < FLT_MIN)
    return false;

// For better accuracy move the coordinates origin to the origin of faceA:
  Vector3D r1 = faceB.pos - faceA.pos;
  Vector3D r2 = faceC.pos - faceA.pos;

  Vector3D vRight(0, r1 & faceB.norm, r2 & faceC.norm); // the right part of the linear system.

// Cramer's rule:
  Matrix3D mr = m;
  mr.set_col(vRight, 0);
  vRes.x = mr.det() / D;

  mr = m;
  mr.set_col(vRight, 1);
  vRes.y = mr.det() / D;

  mr = m;
  mr.set_col(vRight, 2);
  vRes.z = mr.det() / D;

  vRes += faceA.pos; // go back to the absolute coordinates.
  return true;
}

bool CDirichletTesselation::inside(const Vector3D& vPos, const CPlanesSet& vPlanes) const
{
  size_t nFaceCount = vPlanes.size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    const CPlane& plane = vPlanes.at(i);
// The criterion is different from that of CPlane::inside() to make the points belonging to the plane satisfy it.
    if(((vPos - plane.pos) & plane.norm) > Const_Almost_Zero)
      return false;
  }

  return true;
}

double CDirichletTesselation::angle_0_360(const Vector3D& vA, const Vector3D& vB, const Vector3D& vNorm) const
{
  double fPhi = vA ^ vB;  // 0 <= fPhi <= Const_PI.
  double fProj = (vNorm & (vA * vB));
  return fProj > 0 ? fPhi : Const_2PI - fPhi;
}

void CDirichletTesselation::additional_plane(CPlanesSet& vPlanes, CNode3D* pNode) const
{
  Vector3D vNorm(0, 0, 0);
  const CRegionsCollection& vRegs = m_pMesh->get_regions();
  size_t nRegCount = vRegs.size(), nReg, nFace;

  size_t nNbrCount = pNode->vNbrFaces.size(); // count of neighboring FACES.
  for(size_t i = 0; i < nNbrCount; i++)
  {
    nReg = pNode->vNbrFaces.at(i).nReg;
    if(nReg >= nRegCount)
      continue;

    nFace = pNode->vNbrFaces.at(i).nFace;
    if(nFace >= vRegs.at(nReg)->vFaces.size())
      continue;

    vNorm += vRegs.at(nReg)->vFaces.at(nFace)->norm;
  }

  vNorm.normalize();
  CPlane plane(pNode->pos, vNorm);
  vPlanes.push_back(plane);
}

Vector3D CDirichletTesselation::get_norm_vector(size_t nNodeId, size_t nNbrId) const
{
  CNodesCollection& vNodes = m_pMesh->get_nodes();
  CNode3D* pNode = vNodes.at(nNodeId);
  Vector3D v0 = pNode->pos;
  size_t nNbrNodeIndex = pNode->vNbrNodes.at(nNbrId);
  Vector3D v1 = vNodes.at(nNbrNodeIndex)->pos;
  return v1 - v0;
}

void CDirichletTesselation::cell_visualization(const CVertexColl& poly, const Vector3D& vC) const
{
  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  size_t nCount = poly.size();
  for(size_t j = 0; j < nCount; j++)
  {
    pDrawObj->m_vAuxLines.push_back(CEdgeVertex(poly.at(j) + vC));
    size_t j1 = j < nCount - 1 ? j + 1 : 0;
    pDrawObj->m_vAuxLines.push_back(CEdgeVertex(poly.at(j1) + vC));
  }
}

};  // namespace EvaporatingParticle.
