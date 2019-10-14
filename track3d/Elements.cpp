
#include "stdafx.h"

#include <queue>
#include "math.h"
#include "float.h"
#include "ParticleTracking.h"
#include "Elements.h"

#include <algorithm>

namespace EvaporatingParticle
{

//------------------------------------------------------------------------------
// CRay for intersection of an arbitrary ray with the scene.
//------------------------------------------------------------------------------
CRay::CRay(const Vector3D& p, const Vector3D& d)
  : orig(p), dir(d)
{
  dir.normalize();
}

//------------------------------------------------------------------------------
// CBox the base class.
//------------------------------------------------------------------------------
bool CBox::intersect(const CRay& ray) const
{
  Vector3D vRes;
  for(int i = 0; i < 6; i++)
  {
    if(intersect_plane(ray, i, vRes))
      return true;
  }

  return false;
}

bool CBox::intersect_plane(const CRay& ray, int nface, Vector3D& vRes) const
{
  double ksi, x0, y0, z0;

  switch(nface)
  {
    case 0:
    case 1:
    {
// Planes x = vMin.x or x = vMax.x.
      if(fabs(ray.dir.x) < Const_Almost_Zero)
        return false;

      x0 = nface == 0 ? vMin.x : vMax.x;
      ksi = (x0 - ray.orig.x) / ray.dir.x;
      if(ksi < 0)
        return false;

      y0 = ray.orig.y + ksi * ray.dir.y;
      if(y0 < vMin.y || y0 > vMax.y)
        return false;

      z0 = ray.orig.z + ksi * ray.dir.z;
      if(z0 < vMin.z || z0 > vMax.z)
        return false;

      vRes = Vector3D(x0, y0, z0);
      return true;
    }
    case 2:
    case 3:
    {
// Planes y = vMin.y or y = vMax.y.
      if(fabs(ray.dir.y) < Const_Almost_Zero)
        return false;

      y0 = nface == 2 ? vMin.y : vMax.y;
      ksi = (y0 - ray.orig.y) / ray.dir.y;
      if(ksi < 0)
        return false;

      x0 = ray.orig.x + ksi * ray.dir.x;
      if(x0 < vMin.x || x0 > vMax.x)
        return false;

      z0 = ray.orig.z + ksi * ray.dir.z;
      if(z0 < vMin.z || z0 > vMax.z)
        return false;

      vRes = Vector3D(x0, y0, z0);
      return true;
    }
    case 4:
    case 5:
    {
// Planes z = vMin.z or z = vMax.z.
      if(fabs(ray.dir.z) < Const_Almost_Zero)
        return false;

      z0 = nface == 4 ? vMin.z : vMax.z;
      ksi = (z0 - ray.orig.z) / ray.dir.z;
      if(ksi < 0)
        return false;

      x0 = ray.orig.x + ksi * ray.dir.x;
      if(x0 < vMin.x || x0 > vMax.x)
        return false;

      y0 = ray.orig.y + ksi * ray.dir.y;
      if(y0 < vMin.y || y0 > vMax.y)
        return false;

      vRes = Vector3D(x0, y0, z0);
      return true;
    }
  }

  return false;
}

//-------------------------------------------------------------------------------------------------
// CPlane
//-------------------------------------------------------------------------------------------------
bool CPlane::intersect(const CRay& ray, Vector3D& pnt, double& ksi) const
{
  double dscr = norm & ray.dir;
  if(fabs(dscr) < Const_Almost_Zero)
    return false; // the ray is almost parallel to the plane.

  ksi = (norm & (pos - ray.orig)) / dscr;
  if(ksi < 0)
    return false; // the ray is directed away from the plane.

  pnt = ray.orig + ray.dir * ksi; // dir has always a unity length, so that ksi is the distance.
  return true;
}

bool CPlane::operator == (const CPlane& plane) const
{
  if((norm - plane.norm).sqlength() > Const_Almost_Zero)
    return false;

  Vector3D vD = pos - plane.pos;
  if(fabs(vD & norm) > Const_Almost_Zero)
    return false;

  return true;
}

//-------------------------------------------------------------------------------------------------
// CFace - a structure containing some run-time geometrical information in addition to pointers to
//         the three nodes of the triangle.
//-------------------------------------------------------------------------------------------------
CFace::CFace(CNode3D* pn0, CNode3D* pn1, CNode3D* pn2)
  : p0(pn0), p1(pn1), p2(pn2)
{
  init();
}

void CFace::init()
{
  Vector3D e1 = p1->pos - p0->pos;
  Vector3D e2 = p2->pos - p0->pos;
  norm = (e2 * e1).normalized();

  Vector3D vMin(p0->pos), vMax(p0->pos);
  Vector3D vNode[2] = { p1->pos, p2->pos };
  for(UINT i = 0; i < 2; i++)
  {
    if(vNode[i].x < vMin.x)
      vMin.x = vNode[i].x;
    if(vNode[i].x > vMax.x)
      vMax.x = vNode[i].x;
    if(vNode[i].y < vMin.y)
      vMin.y = vNode[i].y;
    if(vNode[i].y > vMax.y)
      vMax.y = vNode[i].y;
    if(vNode[i].z < vMin.z)
      vMin.z = vNode[i].z;
    if(vNode[i].z > vMax.z)
      vMax.z = vNode[i].z;
  }

  const double fTol = 1e-3;
  if(vMax.x - vMin.x < fTol)
  {
    vMax.x += fTol;
    vMin.x -= fTol;
  }
  if(vMax.y - vMin.y < fTol)
  {
    vMax.y += fTol;
    vMin.y -= fTol;
  }
  if(vMax.z - vMin.z < fTol)
  {
    vMax.z += fTol;
    vMin.z -= fTol;
  }

  box.vMin = vMin;
  box.vMax = vMax;
}

bool CFace::intersect(const CRay& ray, double& dist) const
{
  CPlane plane(p0->pos, norm);

  Vector3D vRes;
  if(!plane.intersect(ray, vRes, dist))
    return false;

  if(!box.inside(vRes))
    return false;

// All simple tests are passed. Now try to decompose the input point by the base vectors in the triangle.
  Vector3D e1 = p1->pos - p0->pos;
  Vector3D e2 = p2->pos - p0->pos;
  double e1sqr = e1 & e1;
  double e2sqr = e2 & e2;
  double e1e2 = e1 & e2;

  double det = e1sqr * e2sqr - e1e2 * e1e2;
  if(fabs(det) < Const_Almost_Zero)
    return false;

  double one_ovr_det = 1. / det;

  Vector3D vRel = vRes - p0->pos;
  double vpe1 = vRel & e1;
  double vpe2 = vRel & e2;

  double alpha = (vpe1 * e2sqr - vpe2 * e1e2) * one_ovr_det;
  if(alpha < 0 || alpha > 1)
    return false;

  double beta = (vpe2 * e1sqr - vpe1 * e1e2) * one_ovr_det;
  if(beta < 0 || beta > 1)
    return false;

  if(alpha + beta > 1)
    return false;

  return true;
}

bool CFace::is_wireframe_face() const
{
  return p0->is_wireframe_node() || p1->is_wireframe_node() || p2->is_wireframe_node();
}

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
void CRegion::bounding_box()
{
  Vector3D vMin(FLT_MAX, FLT_MAX, FLT_MAX), vMax(-FLT_MAX, -FLT_MAX, -FLT_MAX), v;
  size_t nFaceCount = vFaces.size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    v = vFaces.at(i)->box.vMin;
    if(v.x < vMin.x)
      vMin.x = v.x;
    if(v.y < vMin.y)
      vMin.y = v.y;
    if(v.z < vMin.z)
      vMin.z = v.z;

    v = vFaces.at(i)->box.vMax;
    if(v.x > vMax.x)
      vMax.x = v.x;
    if(v.y > vMax.y)
      vMax.y = v.y;
    if(v.z > vMax.z)
      vMax.z = v.z;
  }

  box.vMin = vMin;
  box.vMax = vMax;
}

bool CRegion::intersect(const CRay& ray, double& fDist, UINT& nFaceID) const
{
  if(!box.intersect(ray))
    return false;

  size_t nFaceCount = vFaces.size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    CFace* pFace = vFaces.at(i);
    if(pFace->intersect(ray, fDist))
    {
      nFaceID = i;
      return true;
    }
  }

  nFaceID = UINT_MAX;
  return false;
}

bool CRegion::intersect(const CRay& ray, CIntersectColl& results) const
{
  if(!box.intersect(ray))
    return false;

  double dist;
  Vector3D vRes;
  size_t nFaceCount = vFaces.size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    CFace* pFace = vFaces.at(i);
    if(pFace->intersect(ray, dist))
    {
      vRes = ray.orig + ray.dir * dist;   // TO DO: change the prototype of CFace::intersect(...) so that it would output the resulting point.
      CIntersectPoint point(vRes, dist);
      results.push_back(point);
    }
  }

  return results.size() > 0;
}

//-------------------------------------------------------------------------------------------------
// CTreeBox - is an element of the octo-tree of the CTracker.
//-------------------------------------------------------------------------------------------------
CTreeBox::~CTreeBox()
{
  delete_child_boxes();
  vElems.clear();
}

void CTreeBox::collect_elems_inside()
{
  if(pParent == NULL)
    return;

  CElem3D* pElem = NULL;
  CElementsCollection& vParentElems = pParent->vElems;
  size_t nParentElemCount = vParentElems.size(), nVertCount;
  for(size_t i = 0; i < nParentElemCount; i++)
  {
    pElem = vParentElems.at(i);
    nVertCount = pElem->get_node_count();
    for(size_t j = 0; j < nVertCount; j++)
    {
// If at least one vertex of the element is inside the box, include this element into the box collection:
      if(inside(pElem->get_node(j)->pos))
      {
        vElems.push_back(pElem);
        break;
      }
    }
  }
}

void CTreeBox::create_child_boxes()
{
  delete_child_boxes();

  Vector3D vDiag = 0.5 * (vMax - vMin);
  Vector3D vShiftX(vDiag.x, 0, 0);
  Vector3D vShiftY(0, vDiag.y, 0);
  Vector3D vShiftZ(0, 0, vDiag.z);

  pChild = new CTreeBox*[8];

  Vector3D v0;
  for(UINT i = 0; i < 8; i++)
  {
    switch(i)
    {
      case 0: v0 = vMin; break;
      case 1: v0 = vMin + vShiftX; break;
      case 2: v0 = vMin + vShiftY; break;
      case 3: v0 = vMin + vShiftX + vShiftY; break;
      case 4: v0 = vMin + vShiftZ; break;
      case 5: v0 = vMin + vShiftX + vShiftZ; break;
      case 6: v0 = vMin + vShiftY + vShiftZ; break;
      case 7: v0 = vMin + vShiftX + vShiftY + vShiftZ; break;
    }

    pChild[i] = new CTreeBox(v0, v0 + vDiag);
    pChild[i]->nLevel = nLevel + 1;
    pChild[i]->pParent = this;
  }
}

void CTreeBox::delete_child_boxes()
{
  if(pChild != NULL)
  {
    for(UINT i = 0; i < 8; i++)
      delete[] pChild[i];

    delete[] pChild;
    pChild = NULL;
  }
}

static const Vector3D vZero(0, 0, 0);
//-------------------------------------------------------------------------------------------------
// CNode3D
//-------------------------------------------------------------------------------------------------
CElementsCollection CNode3D::get_nbr_elems() const
{
  const CElementsCollection& vGlobElems = CParticleTrackingApp::Get()->GetTracker()->get_elems();

  size_t nNbrCount = vNbrElems.size();
  CElementsCollection vElems(nNbrCount);
  for(size_t i = 0; i < nNbrCount; i++)
    vElems[i] = vGlobElems[vNbrElems[i]];

  return vElems;
}

const CIndexVector & CNode3D::nbr_elems() const
{
	return vNbrElems;
}

bool CNode3D::is_wireframe_node() const
{
  size_t nNbrFacesCount = vNbrFaces.size();
  if(nNbrFacesCount == 0)
    return false;

  UINT nReg = vNbrFaces.at(0).nReg;
  for(size_t i = 1; i < nNbrFacesCount; i++)
  {
    if(vNbrFaces.at(i).nReg != nReg)
      return true;
  }

  return false;
}

//-------------------------------------------------------------------------------------------------
// CElem3D - a base class for all 3D elements - tetrahedra, pyramids, wedges and hexes. 
//-------------------------------------------------------------------------------------------------
void CElem3D::interpolate(const Vector3D& p, CNode3D& node) const
{
  double s, t, u;
  if(!param(p, s, t, u))  // compute parametric coordinates in the element by 3D position.
    return;

  node.pos = p;
  node.set_data(0, 0, 0, 0, 0, 0, vZero, vZero,vZero);  // this is just a container for interpolated data.

  double w;
  CNode3D* pNode = NULL;
  UINT nNodeCount = get_node_count();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    pNode = get_node(i);
    w = shape_func(s, t, u, i);
// Scalars:
    node.dens  += w * pNode->dens;
    node.press += w * pNode->press;
    node.temp  += w * pNode->temp;
    node.visc  += w * pNode->visc;
    node.cond  += w * pNode->cond;
    node.cp    += w * pNode->cp;
// Vectors:
    node.vel   += w * pNode->vel;
    node.field += w * pNode->field;
    node.rf    += w * pNode->rf;
  }
}

bool CElem3D::coeff(const Vector3D& p, double* pWeight) const
{
  double s, t, u;
  if(!param(p, s, t, u))  // compute parametric coordinates in the element by 3D position.
    return false;

  UINT nNodeCount = get_node_count();
  for(size_t i = 0; i < nNodeCount; i++)
    pWeight[i] = shape_func(s, t, u, i);

  return true;
}

void CElem3D::bounding_box()
{
  box.vMin = get_node(0)->pos;
  box.vMax = box.vMin;
  UINT nNodeCount = get_node_count();
  for(size_t i = 1; i < nNodeCount; i++)
  {
    CNode3D* pNode = get_node(i);
    if(pNode->pos.x < box.vMin.x)
      box.vMin.x = pNode->pos.x;
    if(pNode->pos.x > box.vMax.x)
      box.vMax.x = pNode->pos.x;

    if(pNode->pos.y < box.vMin.y)
      box.vMin.y = pNode->pos.y;
    if(pNode->pos.y > box.vMax.y)
      box.vMax.y = pNode->pos.y;

    if(pNode->pos.z < box.vMin.z)
      box.vMin.z = pNode->pos.z;
    if(pNode->pos.z > box.vMax.z)
      box.vMax.z = pNode->pos.z;
  }
}

// Newton method support.
// Discrepancy of positions: computed for a k-th iteration and given point.
Vector3D CElem3D::func(double s, double t, double u, const Vector3D& pos) const
{
  Vector3D vFunc;
  UINT nNodeCount = get_node_count();
  for(size_t i = 0; i < nNodeCount; i++)
    vFunc += shape_func(s, t, u, i) * get_node(i)->pos;

  vFunc -= pos;
  return vFunc;
}

//                            dFx/ds, dFx/dt, dFx/du
// Returned matrix contains:  dFy/ds, dFy/dt, dFy/du
//                            dFz/ds, dFz/dt, dFz/du
//
Matrix3D CElem3D::deriv(double s, double t, double u) const
{
  Vector3D dFds, dFdt, dFdu, vPos;
  UINT nNodeCount = get_node_count();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    vPos = get_node(i)->pos;
    dFds += shape_func_deriv(s, t, u, i, 0) * vPos;
    dFdt += shape_func_deriv(s, t, u, i, 1) * vPos;
    dFdu += shape_func_deriv(s, t, u, i, 2) * vPos;
  }

  return Matrix3D(dFds, dFdt, dFdu);
}

bool CElem3D::param(const Vector3D& p, double& s, double& t, double& u) const
{
  const double cfEps = 0.001;
  const UINT nMaxIterCount = 10;

// Initial approximation:
  elem_middle(s, t, u); // [MS] 02-04-2018, this point must be always inside the element; hopefully this will accelerate the convergence.

  Vector3D vF;
  Matrix3D mdF, mdFR;
  double fDet, fDetR, dKsi[3];

  for(UINT i = 0; i < nMaxIterCount; i++)
  {
// Vector of discrepancies:
    vF = func(s, t, u, p);

// Matrix of derivatives:
    mdF = deriv(s, t, u);
    fDet = mdF.det();
    if(fabs(fDet) < FLT_MIN)
      return false;

// Cramer's rule:
    for(UINT j = 0; j < 3; j++)
    {
      mdFR = mdF;
      mdFR.set_col(-vF, j);
      fDetR = mdFR.det();
      dKsi[j] = fDetR / fDet;
    }

    s += dKsi[0];
    if(s < 0)
      s = 0;
    if(s > 1)
      s = 1;

    t += dKsi[1];
    if(t < 0)
      t = 0;
    if(t > 1)
      t = 1;

    u += dKsi[2];
    if(u < 0)
      u = 0;
    if(u > 1)
      u = 1;

    if((fabs(dKsi[0]) < cfEps) && (fabs(dKsi[1]) < cfEps) && (fabs(dKsi[2]) < cfEps))
      break;
  }

  return true;
}

//-------------------------------------------------------------------------------------------------
//    CElemTetra - a simple object that consists of 4 tetrahedron planes (only planes, not nodes).
//-------------------------------------------------------------------------------------------------
bool CElemTetra::contains(const Vector3D& p) const
{
  for(size_t i = 0; i < 4; i++)
  {
    const CPlane& face = vFaces[i];
    if(!face.inside(p))
      return false;
  }

  return true;
}

void CElemTetra::elem_tetra(const Vector3D& p0, const Vector3D& p1, const Vector3D& p2, const Vector3D& p3)
{
  add_plane(0, p0, p1, p3);
  add_plane(1, p1, p2, p3);
  add_plane(2, p2, p0, p3);
  add_plane(3, p0, p2, p1);
}

void CElemTetra::add_plane(UINT nInd, const Vector3D& p0, const Vector3D& p1, const Vector3D& p2)
{
  Vector3D e1 = p1 - p0;
  Vector3D e2 = p2 - p0;
  Vector3D norm = (e1 * e2).normalized();
  CPlane plane(p0, norm);
  vFaces[nInd] = plane;
}

//-------------------------------------------------------------------------------------------------
//                      CTetra - a tetrahedron object of 4 nodes.
//-------------------------------------------------------------------------------------------------
CTetra::CTetra(const CNodesCollection& vNodes, UINT i0, UINT i1, UINT i2, UINT i3)
  : CElem3D(vNodes)
{
  vTetNodeIds[0] = i0;
  vTetNodeIds[1] = i1;
  vTetNodeIds[2] = i2;
  vTetNodeIds[3] = i3;
  valid = init();
}

CTetra::CTetra(const CNodesCollection& vNodes, const Vector3D& v0, const Vector3D& v1, const Vector3D& v2, const Vector3D& v3)
  : CElem3D(vNodes)
{
  add_plane(0, v0, v1, v3);
  add_plane(1, v1, v2, v3);
  add_plane(2, v2, v0, v3);
  add_plane(3, v0, v2, v1);

  Vector3D v, vVert[4] = { v0, v1, v2, v3 };
  box.vMin = vVert[0];
  box.vMax = vVert[0];
  for(size_t i = 1; i < 4; i++)
  {
    v = vVert[i];
    if(v.x < box.vMin.x)
      box.vMin.x = v.x;
    if(v.x > box.vMax.x)
      box.vMax.x = v.x;

    if(v.y < box.vMin.y)
      box.vMin.y = v.y;
    if(v.y > box.vMax.y)
      box.vMax.y = v.y;

    if(v.z < box.vMin.z)
      box.vMin.z = v.z;
    if(v.z > box.vMax.z)
      box.vMax.z = v.z;
  }
}

bool CTetra::init()
{
  Vector3D p[4] = { get_node(0)->pos, get_node(1)->pos, get_node(2)->pos, get_node(3)->pos };

  add_plane(0, p[0], p[1], p[3]);
  add_plane(1, p[1], p[2], p[3]);
  add_plane(2, p[2], p[0], p[3]);
  add_plane(3, p[0], p[2], p[1]);

  bounding_box();
  return true;
}

bool CTetra::inside(const Vector3D& pos) const
{
  if(!box.inside(pos))
    return false;

  if(!vFaces.contains(pos))
    return false;

  return true;
}

double CTetra::shape_func(double s, double t, double u, size_t node_id) const
{
  switch(node_id)
  {
    case 0: return 1.- s - t - u;
    case 1: return s;
    case 2: return t;
    case 3: return u;
  }

  return 0;
}

double CTetra::shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const
{
  switch(nfunc)
  {
    case 0: return -1.;
    case 1: return nvar == 0 ? 1. : 0.;
    case 2: return nvar == 1 ? 1. : 0.;
    case 3: return nvar == 2 ? 1. : 0.;
  }

  return 0.;
}

void CTetra::add_plane(UINT nInd, const Vector3D& p0, const Vector3D& p1, const Vector3D& p2)
{
  vFaces.add_plane(nInd, p0, p1, p2);
}

const UINT * CTetra::nodes() const
{
	return vTetNodeIds;
}
//[AC 27/03/2017] memory manager
void CTetra::deleteObj()
{
	BlockPool<CTetra>::freeBlock(this);
}
void CPyramid::deleteObj()
{
	BlockPool<CPyramid>::freeBlock(this);
}
void CWedge::deleteObj()
{
	BlockPool<CWedge>::freeBlock(this);
}
void CHexa::deleteObj()
{
	BlockPool<CHexa>::freeBlock(this);
}
//[/AC]

//-------------------------------------------------------------------------------------------------
//                          CPyramid - a pyramid object of 5 nodes.
//-------------------------------------------------------------------------------------------------
CPyramid::CPyramid(const CNodesCollection& vNodes, UINT i0, UINT i1, UINT i2, UINT i3, UINT i4)
  : CElem3D(vNodes) //, pTet1(NULL), pTet2(NULL)
{
  vPyrNodeIds[0] = i0;
	vPyrNodeIds[1] = i1;
  vPyrNodeIds[2] = i2;
  vPyrNodeIds[3] = i3;
  vPyrNodeIds[4] = i4;
  valid = init();
}

CPyramid::~CPyramid()
{
}

bool CPyramid::init()
{
  double fLen02 = (get_node(0)->pos - get_node(2)->pos).sqlength();
  double fLen13 = (get_node(1)->pos - get_node(3)->pos).sqlength();
  if(fLen02 < fLen13)
  {
    vSubTetra[0].elem_tetra(get_node(0)->pos, get_node(1)->pos, get_node(2)->pos, get_node(4)->pos);
    vSubTetra[1].elem_tetra(get_node(0)->pos, get_node(2)->pos, get_node(3)->pos, get_node(4)->pos);
  }
  else
  {
    vSubTetra[0].elem_tetra(get_node(1)->pos, get_node(3)->pos, get_node(0)->pos, get_node(4)->pos);
    vSubTetra[1].elem_tetra(get_node(1)->pos, get_node(2)->pos, get_node(3)->pos, get_node(4)->pos);
  }

  bounding_box(); 
  return true;
}

bool CPyramid::inside(const Vector3D& p) const
{
  if(!box.inside(p))
    return false;

  if(vSubTetra[0].contains(p) || vSubTetra[1].contains(p))
    return true;

  return false;
}

double CPyramid::shape_func(double s, double t, double u, size_t node_id) const
{
  switch(node_id)
  {
    case 0: return (1.0 - s) * (1.0 - t) * (1.0 - u);
    case 1: return s * (1.0 - t) * (1.0 - u);
    case 2: return s * t * (1.0 - u);
    case 3: return (1.0 - s) * t * (1.0 - u);
    case 4: return u;
  }

  return 0;
}

double CPyramid::shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const
{
  switch(nfunc)
  {
    case 0: return nvar == 0 ? -(1.0 - t) * (1.0 - u) : (nvar == 1 ? -(1.0 - s) * (1.0 - u) : -(1.0 - s) * (1.0 - t));
    case 1: return nvar == 0 ?  (1.0 - t) * (1.0 - u) : (nvar == 1 ? -s * (1.0 - u) : -s * (1.0 - t));
    case 2: return nvar == 0 ?  t * (1.0 - u) : (nvar == 1 ? s * (1.0 - u) : -s * t);
    case 3: return nvar == 0 ? -t * (1.0 - u) : (nvar == 1 ? (1.0 - s) * (1.0 - u) : -t * (1.0 - s));
    case 4: return nvar == 2 ? 1.0 : 0.0;
  }

  return 0.;
}

const UINT * CPyramid::nodes() const
{
	return vPyrNodeIds;
}

//-------------------------------------------------------------------------------------------------
//                         CWedge - a prismatic object of 6 nodes.
//-------------------------------------------------------------------------------------------------
CWedge::CWedge(const CNodesCollection& vNodes, UINT i0, UINT i1, UINT i2, UINT i3, UINT i4, UINT i5)
  : CElem3D(vNodes)
{
	vWdgNodeIds[0] = i0;
  vWdgNodeIds[1] = i1;
  vWdgNodeIds[2] = i2;
  vWdgNodeIds[3] = i3;
  vWdgNodeIds[4] = i4;
  vWdgNodeIds[5] = i5;
  valid = init();
}

CWedge::~CWedge()
{
}

bool CWedge::init()
{
  Vector3D vCenter = (get_node(0)->pos + get_node(1)->pos + get_node(2)->pos + get_node(3)->pos + get_node(4)->pos + get_node(5)->pos) / 6.0;

  vSubTetra[0].elem_tetra(get_node(0)->pos, get_node(1)->pos, get_node(2)->pos, vCenter);
  vSubTetra[1].elem_tetra(get_node(5)->pos, get_node(4)->pos, get_node(3)->pos, vCenter);

// Subdivide the non-flat quads by those diagonals that have minimal lengths:
  double fLen04 = (get_node(0)->pos - get_node(4)->pos).sqlength();
  double fLen13 = (get_node(1)->pos - get_node(3)->pos).sqlength();
  if(fLen04 < fLen13)
  {
    vSubTetra[2].elem_tetra(get_node(0)->pos, get_node(4)->pos, get_node(1)->pos, vCenter);
    vSubTetra[3].elem_tetra(get_node(0)->pos, get_node(3)->pos, get_node(4)->pos, vCenter);
  }
  else
  {
    vSubTetra[2].elem_tetra(get_node(1)->pos, get_node(3)->pos, get_node(4)->pos, vCenter);
    vSubTetra[3].elem_tetra(get_node(1)->pos, get_node(0)->pos, get_node(3)->pos, vCenter);
  }

  double fLen15 = (get_node(1)->pos - get_node(5)->pos).sqlength();
  double fLen24 = (get_node(2)->pos - get_node(4)->pos).sqlength();
  if(fLen15 < fLen24)
  {
    vSubTetra[4].elem_tetra(get_node(1)->pos, get_node(5)->pos, get_node(2)->pos, vCenter);
    vSubTetra[5].elem_tetra(get_node(1)->pos, get_node(4)->pos, get_node(5)->pos, vCenter);
  }
  else
  {
    vSubTetra[4].elem_tetra(get_node(2)->pos, get_node(4)->pos, get_node(5)->pos, vCenter);
    vSubTetra[5].elem_tetra(get_node(2)->pos, get_node(1)->pos, get_node(4)->pos, vCenter);
  }

  double fLen23 = (get_node(2)->pos - get_node(3)->pos).sqlength();
  double fLen05 = (get_node(0)->pos - get_node(5)->pos).sqlength();
  if(fLen23 < fLen05)
  {
    vSubTetra[6].elem_tetra(get_node(2)->pos, get_node(3)->pos, get_node(0)->pos, vCenter);
    vSubTetra[7].elem_tetra(get_node(2)->pos, get_node(5)->pos, get_node(3)->pos, vCenter);
  }
  else
  {
    vSubTetra[6].elem_tetra(get_node(0)->pos, get_node(5)->pos, get_node(3)->pos, vCenter);
    vSubTetra[7].elem_tetra(get_node(0)->pos, get_node(2)->pos, get_node(5)->pos, vCenter);
  }

  bounding_box();
  return true;
}

bool CWedge::inside(const Vector3D& p) const
{
  if(!box.inside(p))
    return false;

  for(size_t i = 0; i < 8; i++)
  {
    if(vSubTetra[i].contains(p))
      return true;
  }

  return false;
}

double CWedge::shape_func(double s, double t, double u, size_t node_id) const
{
  switch(node_id)
  {
    case 0: return (1.0 - s - t) * (1.0 - u);
    case 1: return s * (1.0 - u);
    case 2: return t * (1.0 - u);
    case 3: return (1.0 - s - t) * u;
    case 4: return s * u;
    case 5: return t * u;
  }

  return 0;
}

double CWedge::shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const
{
  switch(nfunc)
  {
    case 0: return (nvar == 0) || (nvar == 1) ? u - 1.0 : s + t - 1.0;
    case 1: return nvar == 0 ? 1.0 - u : (nvar == 1 ? 0.0 : -s);
    case 2: return nvar == 0 ? 0.0 : (nvar == 1 ? 1.0 - u : -t);
    case 3: return (nvar == 0) || (nvar == 1) ? -u : 1.0 - s - t;
    case 4: return nvar == 0 ?   u : (nvar == 1 ? 0.0 : s);
    case 5: return nvar == 0 ? 0.0 : (nvar == 1 ?   u : t);
  }

  return 0.;
}

const UINT * CWedge::nodes() const
{
	return vWdgNodeIds;
}

//-------------------------------------------------------------------------------------------------
//                        CHexa - a hexagonal object of 8 nodes.
//-------------------------------------------------------------------------------------------------
CHexa::CHexa(const CNodesCollection& vNodes, UINT i0, UINT i1, UINT i2, UINT i3, UINT i4, UINT i5, UINT i6, UINT i7)
  : CElem3D(vNodes)
{
	vHexNodeIds[0] = i0;
  vHexNodeIds[1] = i1;
  vHexNodeIds[2] = i2;
  vHexNodeIds[3] = i3;
  vHexNodeIds[4] = i4;
  vHexNodeIds[5] = i5;
  vHexNodeIds[6] = i6;
  vHexNodeIds[7] = i7;
  valid = init();
}

bool CHexa::init()
{
  Vector3D vCenter = Vector3D(0, 0, 0);
  for(size_t i = 0; i < 8; i++)
    vCenter += get_node(i)->pos;

  vCenter *= 0.125;

// Subdivide the non-flat quads by those diagonals that have minimal lengths:
  double fLen02 = (get_node(0)->pos - get_node(2)->pos).sqlength();
  double fLen13 = (get_node(1)->pos - get_node(3)->pos).sqlength();
  if(fLen02 < fLen13)
  {
    vSubTetra[0].elem_tetra(get_node(0)->pos, get_node(1)->pos, get_node(2)->pos, vCenter);
    vSubTetra[1].elem_tetra(get_node(0)->pos, get_node(2)->pos, get_node(3)->pos, vCenter);
  }
  else
  {
    vSubTetra[0].elem_tetra(get_node(0)->pos, get_node(1)->pos, get_node(3)->pos, vCenter);
    vSubTetra[1].elem_tetra(get_node(1)->pos, get_node(2)->pos, get_node(3)->pos, vCenter);
  }

  double fLen05 = (get_node(0)->pos - get_node(5)->pos).sqlength();
  double fLen14 = (get_node(1)->pos - get_node(4)->pos).sqlength();
  if(fLen05 < fLen14)
  {
    vSubTetra[2].elem_tetra(get_node(0)->pos, get_node(5)->pos, get_node(1)->pos, vCenter);
    vSubTetra[3].elem_tetra(get_node(0)->pos, get_node(4)->pos, get_node(5)->pos, vCenter);
  }
  else
  {
    vSubTetra[2].elem_tetra(get_node(0)->pos, get_node(4)->pos, get_node(1)->pos, vCenter);
    vSubTetra[3].elem_tetra(get_node(1)->pos, get_node(4)->pos, get_node(5)->pos, vCenter);
  }

  double fLen16 = (get_node(1)->pos - get_node(6)->pos).sqlength();
  double fLen25 = (get_node(2)->pos - get_node(5)->pos).sqlength();
  if(fLen16 < fLen25)
  {
    vSubTetra[4].elem_tetra(get_node(1)->pos, get_node(5)->pos, get_node(6)->pos, vCenter);
    vSubTetra[5].elem_tetra(get_node(1)->pos, get_node(6)->pos, get_node(2)->pos, vCenter);
  }
  else
  {
    vSubTetra[4].elem_tetra(get_node(1)->pos, get_node(5)->pos, get_node(2)->pos, vCenter);
    vSubTetra[5].elem_tetra(get_node(2)->pos, get_node(5)->pos, get_node(6)->pos, vCenter);
  }

  double fLen27 = (get_node(2)->pos - get_node(7)->pos).sqlength();
  double fLen36 = (get_node(3)->pos - get_node(6)->pos).sqlength();
  if(fLen27 < fLen36)
  {
    vSubTetra[6].elem_tetra(get_node(2)->pos, get_node(7)->pos, get_node(3)->pos, vCenter);
    vSubTetra[7].elem_tetra(get_node(6)->pos, get_node(7)->pos, get_node(2)->pos, vCenter);
  }
  else
  {
    vSubTetra[6].elem_tetra(get_node(2)->pos, get_node(6)->pos, get_node(3)->pos, vCenter);
    vSubTetra[7].elem_tetra(get_node(3)->pos, get_node(6)->pos, get_node(7)->pos, vCenter);
  }

  double fLen34 = (get_node(3)->pos - get_node(4)->pos).sqlength();
  double fLen70 = (get_node(7)->pos - get_node(0)->pos).sqlength();
  if(fLen34 < fLen70)
  {
    vSubTetra[8].elem_tetra(get_node(0)->pos, get_node(3)->pos, get_node(4)->pos, vCenter);
    vSubTetra[9].elem_tetra(get_node(4)->pos, get_node(3)->pos, get_node(7)->pos, vCenter);
  }
  else
  {
    vSubTetra[8].elem_tetra(get_node(3)->pos, get_node(7)->pos, get_node(0)->pos, vCenter);
    vSubTetra[9].elem_tetra(get_node(0)->pos, get_node(7)->pos, get_node(4)->pos, vCenter);
  }

  double fLen46 = (get_node(4)->pos - get_node(6)->pos).sqlength();
  double fLen75 = (get_node(7)->pos - get_node(5)->pos).sqlength();
  if(fLen46 < fLen75)
  {
    vSubTetra[10].elem_tetra(get_node(5)->pos, get_node(4)->pos, get_node(6)->pos, vCenter);
    vSubTetra[11].elem_tetra(get_node(6)->pos, get_node(4)->pos, get_node(7)->pos, vCenter);
  }
  else
  {
    vSubTetra[10].elem_tetra(get_node(4)->pos, get_node(7)->pos, get_node(5)->pos, vCenter);
    vSubTetra[11].elem_tetra(get_node(5)->pos, get_node(7)->pos, get_node(6)->pos, vCenter);
  }

  bounding_box();
  return true;
}

bool CHexa::inside(const Vector3D& p) const
{
  if(!box.inside(p))
    return false;

  for(size_t i = 0; i < 12; i++)
  {
    if(vSubTetra[i].contains(p))
      return true;
  }

  return false;
}

double CHexa::shape_func(double s, double t, double u, size_t node_id) const
{
  switch(node_id)
  {
    case 0: return (1.0 - s) * (1.0 - t) * (1.0 - u);
    case 1: return s * (1.0 - t) * (1.0 - u);
    case 2: return s * t * (1.0 - u);
    case 3: return (1.0 - s) * t * (1.0 - u);

    case 4: return (1.0 - s) * (1.0 - t) * u;
    case 5: return s * (1.0 - t) * u;
    case 6: return s * t * u;
    case 7: return (1.0 - s) * t * u;
  }

  return 0;
}

double CHexa::shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const
{
  switch(nfunc)
  {
    case 0: return nvar == 0 ? -(1.0 - t) * (1.0 - u) : (nvar == 1 ? -(1.0 - s) * (1.0 - u) : -(1.0 - s) * (1.0 - t));
    case 1: return nvar == 0 ?  (1.0 - t) * (1.0 - u) : (nvar == 1 ? -s * (1.0 - u) : -s * (1.0 - t));
    case 2: return nvar == 0 ?   t * (1.0 - u) : (nvar == 1 ? s * (1.0 - u) : -s * t);
    case 3: return nvar == 0 ?  -t * (1.0 - u) : (nvar == 1 ? (1.0 - s) * (1.0 - u) : -t * (1.0 - s));

    case 4: return nvar == 0 ?   (t - 1.0) * u : (nvar == 1 ? (s - 1.0) * u : (1.0 - s) * (1.0 - t));
    case 5: return nvar == 0 ?   (1.0 - t) * u : (nvar == 1 ? -s * u : s * (1.0 - t));
    case 6: return nvar == 0 ?   t * u : (nvar == 1 ? s * u : s * t);
    case 7: return nvar == 0 ?  -t * u : (nvar == 1 ? (1.0 - s) * u : (1.0 - s) * t);
  }

  return 0.;
}

const UINT * CHexa::nodes() const
{
	return vHexNodeIds;
}

//-------------------------------------------------------------------------------------------------
//  CSimpleStack
//-------------------------------------------------------------------------------------------------
CSimpleStack::CSimpleStack(UINT nSize)
  : m_pArr(NULL), m_nSize(nSize), m_nPos(0)
{
  init();
}

CSimpleStack::~CSimpleStack()
{
  if(m_pArr != NULL)
    delete[] m_pArr;
}

void CSimpleStack::add(double fVal)
{
  m_pArr[m_nPos] = fVal;
  m_nPos++;
  if(m_nPos == m_nSize)
    m_nPos = 0;
}

double CSimpleStack::get(UINT nOffset) const
{
  return m_pArr != NULL ? m_pArr[(m_nPos + nOffset) % m_nSize] : FLT_MAX;
}

double CSimpleStack::average(UINT nLen) const
{
  double fVal, fAver = 0;
  UINT nCount = 0, j = m_nPos;
  UINT nLength = nLen < m_nSize ? nLen : m_nSize;   // averaging interval.

  for(UINT i = 0; i < nLength; i++)
  {
    fVal = m_pArr[j];
    if(fVal != FLT_MAX)
    {
      fAver += fVal;
      nCount++;
    }

    if(j > 0)
      j--;
    else
      j = m_nSize - 1;
  }

  return nCount > 0 ? fAver / nCount : FLT_MAX;
}

void CSimpleStack::init()
{
  if(m_nSize > 0)
  {
    m_pArr = new double[m_nSize];
    for(UINT i = 0; i < m_nSize; i++)
      m_pArr[i] = FLT_MAX;
  }
}


//-------------------------------------------------------------------------------------------------
// CAverBin
//-------------------------------------------------------------------------------------------------
bool CAverBin::add_val(const CAverValue& vVal)
{
  bool bInside = (vVal.fX >= m_fArgMin) && (vVal.fX < m_fArgMax);
  if(!bInside)
    return false;

  double fCoeff = 1. / (m_nCount + 1); 
  m_vAver.fX = (m_nCount * m_vAver.fX + vVal.fX) * fCoeff;

  switch(m_nAverType)
  {
    case atAverage:
    {
      m_vAver.fY1 = (m_nCount * m_vAver.fY1 + vVal.fY1) * fCoeff;
      m_vAver.fY2 = (m_nCount * m_vAver.fY2 + vVal.fY2) * fCoeff;
      break;
    }
    case atMax:
    {
      if(vVal.fY1 > m_vAver.fY1)
        m_vAver.fY1 = vVal.fY1;
      if(vVal.fY2 > m_vAver.fY2)
        m_vAver.fY2 = vVal.fY2;
      break;
    }
  }

  m_nCount++;
  return true;
}

//-------------------------------------------------------------------------------------------------
// CAveragingEngine
//-------------------------------------------------------------------------------------------------
CAveragingEngine::CAveragingEngine()
{
  m_nAverType = atAverage;
}

CAveragingEngine::~CAveragingEngine()
{
  m_vAverBins.clear();
}

void CAveragingEngine::clear()
{
  m_vAverBins.clear();
}

void CAveragingEngine::set_range(double fArgMin, double fArgMax, UINT nBinsCount)
{
  m_fArgMin = fArgMin;
  m_fArgMax = fArgMax;
  m_fStep = (m_fArgMax - m_fArgMin) / nBinsCount;

  double fBinMin = m_fArgMin, fBinMax;
  for(UINT i = 0; i < nBinsCount; i++)
  {
    fBinMax = fBinMin + m_fStep;
    CAverBin bin;
    bin.set_type(m_nAverType);
    bin.set_bounds(fBinMin, fBinMax);
    m_vAverBins.push_back(bin);
    fBinMin = fBinMax;
  }
}

void CAveragingEngine::add_value(const CAverValue& vVal)
{
  double fArg = vVal.fX;
  if(fArg > m_fArgMax)
    return;
  
  if(fArg > m_fArgMax - Const_Almost_Zero)
    fArg = m_fArgMax - Const_Almost_Zero;

  size_t nBin = size_t(floor((fArg - m_fArgMin) / m_fStep));

  m_vAverBins.at(nBin).add_val(vVal);
}

};  // namespace EvaporatingParticle