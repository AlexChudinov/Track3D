#pragma once

#include "stdafx.h"

#include "math.h"
#include "float.h"
#include "AnsysMesh.h"

#include "DirichletTesselation.h"

#include <map>
#include <algorithm>

#include "../utilities/ParallelFor.h"
#include "ParticleTracking.h" // for the Dirichlet cells visualization.


namespace EvaporatingParticle
{

CDirichletCell::CDirichletCell()
{
  fVolume = 0;
}

CDirichletCell::~CDirichletCell()
{
  delete[] pFaceSquare;
  delete[] pNbrDist;
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

  CDirichletCell* pCell = m_Cells.at(nNodeId);

  CNode3D* pNode = vNodes.at(nNodeId);      // the cell's central node.
  bool bBoundNode = pNode->vNbrFaces.size() != 0;
  if(bBoundNode)
    return get_bound_cell_grad(pCell, pNode, vScalarField);

  double fPhiC = vScalarField.at(nNodeId);  // the value of the potential at the cell's center.

  Vector3D vRes(0, 0, 0);
  double fVol = pCell->fVolume, fFaceSqr, fFacePhi;
  if(fVol < Const_Almost_Zero)
    return Vector3D(0, 0, 0);

  Vector3D vFaceNorm;
  size_t nNbrNodeId;
  size_t nNbrCount = vNodes.at(nNodeId)->vNbrNodes.size();
  for(size_t i = 0; i < nNbrCount; i++)
  {
    fFaceSqr = pCell->pFaceSquare[i];
    nNbrNodeId = pNode->vNbrNodes.at(i);
    fFacePhi = fPhiC + vScalarField.at(nNbrNodeId);
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

  const CNodesCollection& vNodes = m_pMesh->get_nodes();

  //[AC 19042018] changed to parallel
  m_Cells.resize(vNodes.size());
  ThreadPool::splitInPar(vNodes.size(),
	  [&](size_t i) 
  {
	  CNode3D* pNode = vNodes[i];
	  m_Cells[i] = build_cell_in_node(pNode, vNodes);
  },
	  static_cast<CObject*>(this));

  m_bReady = true;
  return !get_terminate_flag();
}

void CDirichletTesselation::clear()
{
  size_t nSize = m_Cells.size();
  for(size_t i = 0; i < nSize; i++)
    delete m_Cells.at(i);

  m_Cells.clear();
}

CDirichletCell* CDirichletTesselation::build_cell_in_node(CNode3D* pNode, const CNodesCollection& vNodes)
{
  CDirichletCell* pCell = new CDirichletCell();

  bool bBoundNode = pNode->vNbrFaces.size() > 0;
  if(bBoundNode)
  {
    init_boundary_cell(pCell, pNode, vNodes);
    return pCell;
  }

// 1. Collect all planes forming the Dirichlet cell around this node.
  size_t nNbrCount = pNode->vNbrNodes.size();
  pCell->pNbrDist = new float[nNbrCount];
  pCell->nFaceCount = nNbrCount;

  size_t nPlanesCount = pCell->nFaceCount;
  CPlanesSet vInnerPlanes;
  vInnerPlanes.reserve(nPlanesCount);

  UINT nId;
  double fDist;
  CNode3D* pNbrNode = NULL;
  Vector3D vC = pNode->pos, vOrg, vNorm, vNbr;  // cell center, plane origin and normal and neighbour cell center, respectively.

  for(size_t i = 0; i < nNbrCount; i++)
  {
    nId = pNode->vNbrNodes.at(i);
    pNbrNode = vNodes.at(nId);
    vNbr = pNbrNode->pos;
    vOrg = 0.5 * (vC + vNbr);
    vNorm = vNbr - vC;          // plane normal looks always out of the cell.
    fDist = vNorm.length();
    pCell->pNbrDist[i] = fDist;
    vNorm /= fDist;

    CPlane face(vOrg, vNorm);
    vInnerPlanes.push_back(face);
  }

// 2. Collect all vertices of the Dirichlet cell as all possible intersections of triads of different planes.
//    Exclude the points lying outside the cell.
  CVertexColl vCellVert;
  collect_cell_vertices(vInnerPlanes, vCellVert);

// 3. For each plane select vertices belonging to the plane and build the cell's polygon. Compute the polygon's square.
  double fS, fH;
  pCell->fVolume = 0;
  pCell->pFaceSquare = new float[nNbrCount];
  for(size_t j = 0; j < nNbrCount; j++)
  {
    const CPlane& plane = vInnerPlanes.at(j);
    fS = build_polygon_in_plane(plane, vCellVert); // the square of that face of the cell, which belongs to this plane.
    fH = 0.5 * pCell->pNbrDist[j];

    pCell->pFaceSquare[j] = fS;
    pCell->fVolume += Const_One_Third * fS * fH;  // V = 1/3 * S * H, a pyramid volume.
  }

  return pCell;
}

void CDirichletTesselation::collect_cell_vertices(const CPlanesSet& vInnerPlanes, CVertexColl& vCellVert) const
{
  Vector3D vVert;
  size_t nPlanesCount = vInnerPlanes.size();
  for(size_t i = 0; i < nPlanesCount; i++)
  {
    const CPlane& PlaneA = vInnerPlanes.at(i);
    for(size_t j = 0; j < i; j++)
    {
      if(j == i)
        continue;

      const CPlane& PlaneB = vInnerPlanes.at(j);
      for(size_t k = 0; k < j; k++)
      {
        if(k == j || k == i)
          continue;

        const CPlane& PlaneC = vInnerPlanes.at(k);
        if(!three_planes_intersect(PlaneA, PlaneB, PlaneC, vVert))
          continue;

        if(!inside(PlaneA, PlaneB, PlaneC, vVert, vInnerPlanes))
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
  size_t nPolyCount = vPoly.size(), k1;
  for(size_t k = 0; k < nPolyCount; k++)
  {
    k1 = k < nPolyCount - 1 ? k + 1 : 0;
    fSquare += 0.5 * (vPoly.at(k) * vPoly.at(k1)).length();
  }

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

bool CDirichletTesselation::inside(const CPlane& a, const CPlane& b, const CPlane& c, const Vector3D& vPos, const CPlanesSet& vPlanes) const
{
  size_t nFaceCount = vPlanes.size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    const CPlane& plane = vPlanes.at(i);
    if(plane == a || plane == b || plane == c)
      continue;

    if(!plane.inside(vPos))
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

bool CDirichletTesselation::bound_deriv_calc_cond(const Vector3D& vNorm, const Vector3D& e)
{
  static const double fEps = 0.15;
  return (e & vNorm) < fEps;
}


void CDirichletTesselation::init_boundary_cell(CDirichletCell* pCell, CNode3D* pNode, const CNodesCollection& vNodes) const
{
  CNode3D* pNbrNode = NULL;
  Vector3D vC = pNode->pos, e;
  float Axx = 0, Axy = 0, Axz = 0, Ayy = 0, Ayz = 0, Azz = 0, fDet, fDist;

  size_t nNbrCount = pNode->vNbrNodes.size();

  pCell->nFaceCount = nNbrCount;
  pCell->pNbrDist = new float[nNbrCount];
  pCell->pFaceSquare = new float[10];      // Cxx, Cxy, Cxz, Cyy, Cyz, Czz, fDet, nx, ny, nz.

  Vector3D vNorm = get_bound_norm(pNode);
  pCell->pFaceSquare[7] = (float)vNorm.x;
  pCell->pFaceSquare[8] = (float)vNorm.y;
  pCell->pFaceSquare[9] = (float)vNorm.z;

// DEBUG
//  int nCount = 0;
// END DEBUG

  for(size_t i = 0; i < nNbrCount; i++)
  {
    pNbrNode = vNodes.at(pNode->vNbrNodes.at(i));
    e = pNbrNode->pos - vC;
    fDist = e.length();
    pCell->pNbrDist[i] = fDist;
    e /= fDist;

// Take into account only those neighbours, for which (vNorm, e) < some eps.
    if(!bound_deriv_calc_cond(vNorm, e))
      continue;

    Axx += e.x * e.x;
    Axy += e.x * e.y;
    Axz += e.x * e.z;
    Ayy += e.y * e.y;
    Ayz += e.y * e.z;
    Azz += e.z * e.z;

// DEBUG
//    nCount++;
// END DEBUG
  }

// DEBUG
//  if(nCount < 3)
//    AfxMessageBox("Too little items in the overdetermined system's sum.");
// END DEBUG

  Matrix3D m(Axx, Axy, Axz, Axy, Ayy, Ayz, Axz, Ayz, Azz);
  fDet = m.det();

  float Cxx = Ayy * Azz - Ayz * Ayz;
  float Cxy = Axz * Ayz - Axy * Azz;
  float Cxz = Axy * Ayz - Axz * Ayy;
  float Cyy = Axx * Azz - Axz * Axz;
  float Cyz = Axy * Axz - Axx * Ayz;
  float Czz = Axx * Ayy - Axy * Axy;

  pCell->pFaceSquare[0] = Cxx;
  pCell->pFaceSquare[1] = Cxy;
  pCell->pFaceSquare[2] = Cxz;
  pCell->pFaceSquare[3] = Cyy;
  pCell->pFaceSquare[4] = Cyz;
  pCell->pFaceSquare[5] = Czz;
  pCell->pFaceSquare[6] = fDet;
}

Vector3D CDirichletTesselation::get_bound_norm(CNode3D* pBoundNode) const
{
  size_t nNbrCount = pBoundNode->vNbrFaces.size();
  if(nNbrCount == 0)
    return scvNull;

  const CRegionsCollection& vRegs = m_pMesh->get_regions(false);
  size_t nRegCount = vRegs.size();

  UINT nReg, nFace;
  CFace* pFace = NULL;
  Vector3D vNorm(0, 0, 0);
  for(size_t i = 0; i < nNbrCount; i++)
  {
    nReg = pBoundNode->vNbrFaces.at(i).nReg;
    if(nReg >= nRegCount)
      continue;

    nFace = pBoundNode->vNbrFaces.at(i).nFace;
    if(nFace >= vRegs.at(nReg)->vFaces.size())
      continue;

    vNorm += vRegs.at(nReg)->vFaces.at(nFace)->norm;
  }

  vNorm.normalize();
  return vNorm;
}

Matrix3D CDirichletTesselation::get_bound_cell_mtx(CDirichletCell* pCell) const
{
  Matrix3D m(pCell->pFaceSquare[0], pCell->pFaceSquare[1], pCell->pFaceSquare[2],
             pCell->pFaceSquare[1], pCell->pFaceSquare[3], pCell->pFaceSquare[4],
             pCell->pFaceSquare[2], pCell->pFaceSquare[4], pCell->pFaceSquare[5]);

  return m;
}

Vector3D CDirichletTesselation::get_bound_cell_grad(CDirichletCell* pCell, CNode3D* pNode, const std::vector<float>& vScalarField) const
{
  float fPhi, fPhiC = vScalarField.at(pNode->nInd);
  Vector3D b(0, 0, 0), vNrm(pCell->pFaceSquare[7], pCell->pFaceSquare[8], pCell->pFaceSquare[9]), e;
  for(size_t i = 0; i < pCell->nFaceCount; i++)
  {
    fPhi = vScalarField.at(pNode->vNbrNodes.at(i));

    e = get_neighbor_vector(pNode->nInd, i) / pCell->pNbrDist[i];
    if(!bound_deriv_calc_cond(vNrm, e))
      continue;

    b += e * ((fPhi - fPhiC) / pCell->pNbrDist[i]);
  }

  Matrix3D m = get_bound_cell_mtx(pCell); // matrix C.
  Vector3D vG = m * b;
  vG /= pCell->pFaceSquare[6];  // determinant of matrix A.
  return vG;
}

Vector3D CDirichletTesselation::get_neighbor_vector(size_t nNodeId, size_t nNbrId) const
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
