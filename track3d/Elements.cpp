
#include "stdafx.h"

#include <queue>
#include "math.h"
#include "float.h"
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
  for(int i = 0; i < 6; i++)
  {
    if(intersect_plane(ray, i))
      return true;
  }

  return false;
}

bool CBox::intersect_plane(const CRay& ray, int nface) const
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

bool CRegion::intersect(const CRay& ray, double& dist) const
{
  if(!box.intersect(ray))
    return false;

  size_t nFaceCount = vFaces.size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    CFace* pFace = vFaces.at(i);
    if(pFace->intersect(ray, dist))
      return true;
  }

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
    nVertCount = pElem->vNodes.size();
    for(size_t j = 0; j < nVertCount; j++)
    {
// If at least one vertex of the element is inside the box, include this element into the box collection:
      if(inside(pElem->vNodes.at(j)->pos))
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
// CElem3D - a base class for all 3D elements - tetrahedra, pyramids, wedges and hexes. 
//-------------------------------------------------------------------------------------------------
bool CElem3D::inside(const Vector3D& pos) const
{
  if(!box.inside(pos))
    return false;

  Vector3D p;
  size_t nFacesCount = vFaces.size();
  for(size_t i = 0; i < nFacesCount; i++)
  {
    const CPlane& face = vFaces.at(i);
    if(!face.inside(pos))
      return false;
  }

  return true;
}

void CElem3D::interpolate(const Vector3D& p, CNode3D& node) const
{
  double s, t, u;
  if(!param(p, s, t, u))  // compute parametric coordinates in the element by 3D position.
    return;

  node.pos = p;
  node.set_data(0, 0, 0, 0, 0, 0, vZero, vZero,vZero);  // this is just a container for interpolated data.

  double w;
  CNode3D* pNode = NULL;
  size_t nNodesCount = vNodes.size();
  for(size_t i = 0; i < nNodesCount; i++)
  {
    pNode = vNodes.at(i);
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

  size_t nNodesCount = vNodes.size();
  for(size_t i = 0; i < nNodesCount; i++)
    pWeight[i] = shape_func(s, t, u, i);

  return true;
}

void CElem3D::bounding_box()
{
  box.vMin = vNodes.at(0)->pos;
  box.vMax = vNodes.at(0)->pos;
  size_t nNodesCount = vNodes.size();
  for(size_t i = 1; i < nNodesCount; i++)
  {
    CNode3D* pNode = vNodes.at(i);
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
  size_t nNodeCount = vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
    vFunc += shape_func(s, t, u, i) * vNodes.at(i)->pos;

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
  size_t nNodeCount = vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    vPos = vNodes.at(i)->pos;
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
  s = 0.49, t = 0.49, u = 0.49;

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

void CElem3D::add_plane(const Vector3D& p0, const Vector3D& p1, const Vector3D& p2)
{
  Vector3D e1 = p1 - p0;
  Vector3D e2 = p2 - p0;
  Vector3D norm = (e1 * e2).normalized();
  CPlane plane(p0, norm);
  vFaces.push_back(plane);
}

//-------------------------------------------------------------------------------------------------
//                      CTetra - a tetrahedron object of 4 nodes.
//-------------------------------------------------------------------------------------------------
CTetra::CTetra(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3)
{
	//Try if we can reduce memory consumption if we use explicit vector allocation
	vNodes.reserve(4);
	vNodes.push_back(p0);
  vNodes.push_back(p1);
  vNodes.push_back(p2);
  vNodes.push_back(p3);
  valid = init();
}

bool CTetra::init()
{
  vFaces.clear();
  Vector3D p[4] = { vNodes.at(0)->pos, vNodes.at(1)->pos, vNodes.at(2)->pos, vNodes.at(3)->pos };

  add_plane(p[0], p[1], p[3]);
  add_plane(p[1], p[2], p[3]);
  add_plane(p[2], p[0], p[3]);
  add_plane(p[0], p[2], p[1]);

  bounding_box();
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

//-------------------------------------------------------------------------------------------------
//                          CPyramid - a pyramid object of 5 nodes.
//-------------------------------------------------------------------------------------------------
CPyramid::CPyramid(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4)
  : pTet1(NULL), pTet2(NULL)
{
	vNodes.reserve(5);
  vNodes.push_back(p0);
  vNodes.push_back(p1);
  vNodes.push_back(p2);
  vNodes.push_back(p3);
  vNodes.push_back(p4);
  valid = init();
}

CPyramid::~CPyramid()
{
  if(pTet1 != NULL)
    delete pTet1;
  if(pTet2 != NULL)
    delete pTet2;
}

bool CPyramid::init()
{
  double fLen02 = (vNodes.at(0)->pos - vNodes.at(2)->pos).sqlength();
  double fLen13 = (vNodes.at(1)->pos - vNodes.at(3)->pos).sqlength();
  if(fLen02 < fLen13)
  {
    pTet1 = new CTetra(vNodes.at(0), vNodes.at(1), vNodes.at(2), vNodes.at(4));
    pTet2 = new CTetra(vNodes.at(0), vNodes.at(2), vNodes.at(3), vNodes.at(4));
  }
  else
  {
    pTet1 = new CTetra(vNodes.at(1), vNodes.at(3), vNodes.at(0), vNodes.at(4));
    pTet2 = new CTetra(vNodes.at(1), vNodes.at(2), vNodes.at(3), vNodes.at(4));
  }

  bounding_box(); 
  return true;
}

bool CPyramid::inside(const Vector3D& p) const
{
  if(!box.inside(p))
    return false;

  if(pTet1->inside(p) || pTet2->inside(p))
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

//-------------------------------------------------------------------------------------------------
//                         CWedge - a prismatic object of 6 nodes.
//-------------------------------------------------------------------------------------------------
CWedge::CWedge(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4, CNode3D* p5)
  : pCenter(NULL)
{
	vNodes.reserve(6);
  vNodes.push_back(p0);
  vNodes.push_back(p1);
  vNodes.push_back(p2);
  vNodes.push_back(p3);
  vNodes.push_back(p4);
  vNodes.push_back(p5);
  valid = init();
}

CWedge::~CWedge()
{
  size_t nSize = vTetras.size();
  if(nSize > 0)
  {
    for(size_t i = 0; i < nSize; i++)
      delete vTetras.at(i);
  }

  vTetras.clear();

  if(pCenter != NULL)
    delete pCenter;
}

bool CWedge::init()
{
  Vector3D vC = (vNodes.at(0)->pos + vNodes.at(1)->pos + vNodes.at(2)->pos + vNodes.at(3)->pos + vNodes.at(4)->pos + vNodes.at(5)->pos) / 6.0;
  pCenter = new CNode3D(vC);
  vTetras.reserve(8);
  vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(1), vNodes.at(2), pCenter));
  vTetras.push_back(new CTetra(vNodes.at(5), vNodes.at(4), vNodes.at(3), pCenter));

// Subdivide the non-flat quads by those diagonals that have minimal lengths:
  double fLen04 = (vNodes.at(0)->pos - vNodes.at(4)->pos).sqlength();
  double fLen13 = (vNodes.at(1)->pos - vNodes.at(3)->pos).sqlength();
  if(fLen04 < fLen13)
  {
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(4), vNodes.at(1), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(3), vNodes.at(4), pCenter));
  }
  else
  {
    vTetras.push_back(new CTetra(vNodes.at(1), vNodes.at(3), vNodes.at(4), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(1), vNodes.at(0), vNodes.at(3), pCenter));
  }

  double fLen15 = (vNodes.at(1)->pos - vNodes.at(5)->pos).sqlength();
  double fLen24 = (vNodes.at(2)->pos - vNodes.at(4)->pos).sqlength();
  if(fLen15 < fLen24)
  {
    vTetras.push_back(new CTetra(vNodes.at(1), vNodes.at(5), vNodes.at(2), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(1), vNodes.at(4), vNodes.at(5), pCenter));
  }
  else
  {
    vTetras.push_back(new CTetra(vNodes.at(2), vNodes.at(4), vNodes.at(5), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(2), vNodes.at(1), vNodes.at(4), pCenter));
  }

  double fLen23 = (vNodes.at(2)->pos - vNodes.at(3)->pos).sqlength();
  double fLen05 = (vNodes.at(0)->pos - vNodes.at(5)->pos).sqlength();
  if(fLen23 < fLen05)
  {
    vTetras.push_back(new CTetra(vNodes.at(2), vNodes.at(3), vNodes.at(0), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(2), vNodes.at(5), vNodes.at(3), pCenter));
  }
  else
  {
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(5), vNodes.at(3), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(2), vNodes.at(5), pCenter));
  }

  bounding_box();
  return true;
}

bool CWedge::inside(const Vector3D& p) const
{
  if(!box.inside(p))
    return false;

  size_t nSize = vTetras.size();
  for(size_t i = 0; i < nSize; i++)
  {
    if(vTetras.at(i)->inside(p))
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

//-------------------------------------------------------------------------------------------------
//                        CHexa - a hexagonal object of 8 nodes.
//-------------------------------------------------------------------------------------------------
CHexa::CHexa(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4, CNode3D* p5, CNode3D* p6, CNode3D* p7)
{
	vNodes.reserve(8);
  vNodes.push_back(p0);
  vNodes.push_back(p1);
  vNodes.push_back(p2);
  vNodes.push_back(p3);
  vNodes.push_back(p4);
  vNodes.push_back(p5);
  vNodes.push_back(p6);
  vNodes.push_back(p7);
  valid = init();
}

bool CHexa::init()
{
  Vector3D vC(0, 0, 0);
  size_t nNodeCount = vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
    vC += vNodes.at(i)->pos;

  vC /= (double)nNodeCount;
  pCenter = new CNode3D(vC);
  vTetras.reserve(12);

// Subdivide the non-flat quads by those diagonals that have minimal lengths:
  double fLen02 = (vNodes.at(0)->pos - vNodes.at(2)->pos).sqlength();
  double fLen13 = (vNodes.at(1)->pos - vNodes.at(3)->pos).sqlength();
  if(fLen02 < fLen13)
  {
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(1), vNodes.at(2), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(2), vNodes.at(3), pCenter));
  }
  else
  {
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(1), vNodes.at(3), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(1), vNodes.at(2), vNodes.at(3), pCenter));
  }

  double fLen05 = (vNodes.at(0)->pos - vNodes.at(5)->pos).sqlength();
  double fLen14 = (vNodes.at(1)->pos - vNodes.at(4)->pos).sqlength();
  if(fLen05 < fLen14)
  {
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(5), vNodes.at(1), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(4), vNodes.at(5), pCenter));
  }
  else
  {
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(4), vNodes.at(1), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(1), vNodes.at(4), vNodes.at(5), pCenter));
  }

  double fLen16 = (vNodes.at(1)->pos - vNodes.at(6)->pos).sqlength();
  double fLen25 = (vNodes.at(2)->pos - vNodes.at(5)->pos).sqlength();
  if(fLen16 < fLen25)
  {
    vTetras.push_back(new CTetra(vNodes.at(1), vNodes.at(5), vNodes.at(6), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(1), vNodes.at(6), vNodes.at(2), pCenter));
  }
  else
  {
    vTetras.push_back(new CTetra(vNodes.at(1), vNodes.at(5), vNodes.at(2), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(2), vNodes.at(5), vNodes.at(6), pCenter));
  }

  double fLen27 = (vNodes.at(2)->pos - vNodes.at(7)->pos).sqlength();
  double fLen36 = (vNodes.at(3)->pos - vNodes.at(6)->pos).sqlength();
  if(fLen27 < fLen36)
  {
    vTetras.push_back(new CTetra(vNodes.at(2), vNodes.at(7), vNodes.at(3), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(6), vNodes.at(7), vNodes.at(2), pCenter));
  }
  else
  {
    vTetras.push_back(new CTetra(vNodes.at(2), vNodes.at(6), vNodes.at(3), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(3), vNodes.at(6), vNodes.at(7), pCenter));
  }

  double fLen34 = (vNodes.at(3)->pos - vNodes.at(4)->pos).sqlength();
  double fLen70 = (vNodes.at(7)->pos - vNodes.at(0)->pos).sqlength();
  if(fLen34 < fLen70)
  {
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(3), vNodes.at(4), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(4), vNodes.at(3), vNodes.at(7), pCenter));
  }
  else
  {
    vTetras.push_back(new CTetra(vNodes.at(3), vNodes.at(7), vNodes.at(0), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(0), vNodes.at(7), vNodes.at(4), pCenter));
  }

  double fLen46 = (vNodes.at(4)->pos - vNodes.at(6)->pos).sqlength();
  double fLen75 = (vNodes.at(7)->pos - vNodes.at(5)->pos).sqlength();
  if(fLen46 < fLen75)
  {
    vTetras.push_back(new CTetra(vNodes.at(5), vNodes.at(4), vNodes.at(6), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(6), vNodes.at(4), vNodes.at(7), pCenter));
  }
  else
  {
    vTetras.push_back(new CTetra(vNodes.at(4), vNodes.at(7), vNodes.at(5), pCenter));
    vTetras.push_back(new CTetra(vNodes.at(5), vNodes.at(7), vNodes.at(6), pCenter));
  }
/*
  CNode3D* pBaseNode = vNodes.at(0);
  Vector3D e1 = vNodes.at(1)->pos - pBaseNode->pos;
  Vector3D e2 = vNodes.at(3)->pos - pBaseNode->pos;
  Vector3D e3 = vNodes.at(4)->pos - pBaseNode->pos;

  vFaces.clear();

// Planes defining the faces:
  Vector3D norm = (e2 * e1).normalized();
  CPlane plane1(pBaseNode->pos, norm);
  vFaces.push_back(plane1);

  norm = (e1 * e3).normalized();
  CPlane plane2(pBaseNode->pos, norm);
  vFaces.push_back(plane2);

  norm = (e3 * e2).normalized();
  CPlane plane3(pBaseNode->pos, norm);
  vFaces.push_back(plane3);

  pBaseNode = vNodes.at(7);
  e1 = vNodes.at(6)->pos - pBaseNode->pos;
  e2 = vNodes.at(4)->pos - pBaseNode->pos;
  e3 = vNodes.at(3)->pos - pBaseNode->pos;

  norm = (e2 * e1).normalized();
  CPlane plane4(pBaseNode->pos, norm);
  vFaces.push_back(plane4);

  norm = (e3 * e2).normalized();
  CPlane plane5(pBaseNode->pos, norm);
  vFaces.push_back(plane5);

  norm = (e1 * e3).normalized();
  CPlane plane6(pBaseNode->pos, norm);
  vFaces.push_back(plane6);
*/
  bounding_box();
  return true;
}

bool CHexa::inside(const Vector3D& p) const
{
  if(!box.inside(p))
    return false;

  size_t nSize = vTetras.size();
  for(size_t i = 0; i < nSize; i++)
  {
    if(vTetras.at(i)->inside(p))
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