
#include "stdafx.h"

#include "Primitives.h"
#include "ParticleTracking.h"

namespace EvaporatingParticle
{

//-------------------------------------------------------------------------------------------------
// CPrimitive
//-------------------------------------------------------------------------------------------------
void CPrimitive::set_default()
{
  m_bReady = false;
}

bool CPrimitive::get_regions(const CStringVector& vRegNames)
{
  m_vRegions.clear();
  size_t nRegCount = vRegNames.size();
  if(nRegCount == 0)
    return false;

  m_vRegions.reserve(nRegCount);
  for(size_t i = 0; i < nRegCount; i++)
  {
    CRegion* pReg = CAnsysMesh::get_region(vRegNames[i]);
    if(pReg != NULL)
      m_vRegions.push_back(pReg);
  }

  return m_vRegions.size() > 0;
}

bool CPrimitive::collect_reg_nodes()
{
  m_vRegNodes.clear();
  CFace* pFace = NULL;

  size_t nRegCount = m_vRegions.size();
  if(nRegCount == 0)
    return false;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  size_t nNodesCount = pObj->get_nodes().size();
  if(nNodesCount == 0)
    return false;

  std::vector<bool> vbSel(nNodesCount, false);

  for(size_t j = 0; j < nRegCount; j++)
  {
    CRegion* pReg = m_vRegions.at(j);
    size_t nFacesCount = pReg->vFaces.size();
    if(nFacesCount == 0)
      continue;

    for(size_t i = 0; i < nFacesCount; i++)
    {
// Nodes in m_vRegNodes must be unique.
      pFace = pReg->vFaces.at(i);
      if(!vbSel[pFace->p0->nInd])
      {
        m_vRegNodes.push_back(pFace->p0);
        vbSel[pFace->p0->nInd] = true;
      }
      if(!vbSel[pFace->p1->nInd])
      {
        m_vRegNodes.push_back(pFace->p1);
        vbSel[pFace->p1->nInd] = true;
      }
      if(!vbSel[pFace->p2->nInd])
      {
        m_vRegNodes.push_back(pFace->p2);
        vbSel[pFace->p2->nInd] = true;
      }
    }
  }

  vbSel.clear();
  m_vRegNodes.shrink_to_fit();
  return m_vRegNodes.size() > 3;
}


//-------------------------------------------------------------------------------------------------
// CPrimitiveBox
//-------------------------------------------------------------------------------------------------
CPrimitiveBox::CPrimitiveBox(const CStringVector& vRegNames, double fHeight)
{
  set_default();
  set_height(fHeight);
  if(get_regions(vRegNames))
    calc_box_params();
}

CPrimitiveBox::CPrimitiveBox(const Vector3D& vOrig, const Vector3D& vX, const Vector3D& vY, const Vector3D& vZ)
{
  create_box(vOrig, vX, vY, vZ);
}

Vector2D CPrimitiveBox::model_coord(const Vector3D& vPos) const
{
// It is assumed that contains(vPos) has returned "true".
  Vector3D v0 = vPos - m_vOrig;
  v0.z = 0;
  return Vector2D((v0 & m_vLocX), (v0 & m_vLocY));
}

bool CPrimitiveBox::contains(const Vector3D& vPos) const
{
  if(!m_bReady)
    return false;

  for(UINT i = 0; i < 6; i++)
  {
    if(!m_Planes[i].inside(vPos))
      return false;
  }

  return true;
}

void CPrimitiveBox::calc_box_params()
{
  calc_loc_triad();
  if(!collect_reg_nodes())
    return;

  Vector3D vYw(0, 1, 0), vZw(0, 0, 1);
  double fAngle = m_vLocY ^ vYw;          // fAngle must be in the range [0, pi].
  double fBeta = ((vYw * m_vLocY) & vZw);
  if(fBeta < -Const_Almost_Zero)
    fAngle = -fAngle;

  Matrix3D wrld2loc = Matrix3D::rot(vZw, fAngle);

// Build a usual bounding box in the local c.s.
  Vector3D vMin(FLT_MAX, FLT_MAX, FLT_MAX), vMax = -vMin, vPos;

  size_t nRegNodesCount = m_vRegNodes.size();
  for(size_t i = 0; i < nRegNodesCount; i++)
  {
    vPos = wrld2loc * m_vRegNodes.at(i)->pos;
    if(vPos.x < vMin.x)
      vMin.x = vPos.x;
    if(vPos.x > vMax.x)
      vMax.x = vPos.x;

    if(vPos.y < vMin.y)
      vMin.y = vPos.y;
    if(vPos.y > vMax.y)
      vMax.y = vPos.y;

    if(vPos.z < vMin.z)
      vMin.z = vPos.z;
    if(vPos.z > vMax.z)
      vMax.z = vPos.z;
  }

// Return back to the world c.s.:
  Matrix3D loc2wrld = Matrix3D::rot(vZw, -fAngle);

  Vector3D vOrig = loc2wrld * vMin;
  Vector3D vOppos = loc2wrld * vMax;

  double fdX = (vOppos - vOrig) & m_vLocX;
  Vector3D vX = fdX * m_vLocX;

  double fdZ = (vOppos - vOrig) & m_vLocZ;
  Vector3D vZ = fdZ * m_vLocZ;

  Vector3D vY = m_fHeight * m_vLocY;

  m_vOrig = vOrig;
  create_box(m_vOrig, vX, vY, vZ);
  m_bReady = true;
}

void CPrimitiveBox::calc_loc_triad()
{
  m_vLocZ = Vector3D(0, 0, 1);

// m_vLocY is the average of -pFace->norm over all faces.
  CFace* pFace = NULL;
  CRegion* pReg = NULL;
  Vector3D vAverNorm(0, 0, 0);
  size_t nRegCount = m_vRegions.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    pReg = m_vRegions.at(i);
    size_t nFaceCount = pReg->vFaces.size();
    for(size_t j = 0; j < nFaceCount; j++)
    {
      pFace = pReg->vFaces.at(j);
      vAverNorm += pFace->norm;
    }
  }

  m_vLocY = -vAverNorm.normalized();
  m_vLocX = m_vLocY * m_vLocZ;
}

void CPrimitiveBox::create_box(const Vector3D& vOrig, const Vector3D& vX, const Vector3D& vY, const Vector3D& vZ)
{
  Vector3D vOrigin = vOrig;
  m_Planes[0].pos = vOrigin;
  m_Planes[0].norm = -m_vLocZ;
  m_Planes[1].pos = vOrigin;
  m_Planes[1].norm = -m_vLocY;
  m_Planes[2].pos = vOrigin;
  m_Planes[2].norm = -m_vLocX;
  vOrigin = vOrig + vX + vY + vZ;
  m_Planes[3].pos = vOrigin;
  m_Planes[3].norm = m_vLocZ;
  m_Planes[4].pos = vOrigin;
  m_Planes[4].norm = m_vLocY;
  m_Planes[5].pos = vOrigin;
  m_Planes[5].norm = m_vLocX;
}

Vector3D CPrimitiveBox::get_loc_normal(const Vector3D& vPos) const
{
  return m_vLocY;
}

//-------------------------------------------------------------------------------------------------
// CPrimCylSector - a sectorial volume between two cylinders.
//-------------------------------------------------------------------------------------------------
CPrimCylSector::CPrimCylSector()
  : m_vCylAxis(Vector3D(0, 0, 0), Vector3D(0, 0, 1))
{
  set_default();
}

CPrimCylSector::CPrimCylSector(const CStringVector& vRegNames, double fX0, double fY0, double fR, double fHeight)
  : m_vCylAxis(Vector3D(0, 0, 0), Vector3D(0, 0, 1))
{
  set_default();
  set_height(fHeight);
  if(get_regions(vRegNames) && collect_reg_nodes())
  {
    assign(fX0, fY0, fR);
    m_bReady = true;
  }
//    calc_cyl_params();
}

Vector2D CPrimCylSector::model_coord(const Vector3D& vPos) const
{
// It is assumed that contains(vPos) has returned "true".
  Vector3D vP = vPos - m_vCylAxis.orig;
  vP.z = 0;
  double fR = vP.length();
  vP /= fR;

  double fArg = (m_vOrigDir & vP);
  if(fArg > 1.0)
    fArg = 1.0;
  if(fArg < -1.0)
    fArg = -1.0;

  double fPhi = acos(fArg);

// Note that we return here the "horizontal" coordinate as if the point is lying on the surface.
  return Vector2D(m_fSmallR * fPhi, fR - m_fSmallR);
}

Vector3D CPrimCylSector::get_loc_normal(const Vector3D& vPos) const
{
  Vector3D v = vPos - m_vCylAxis.orig;
  v.z = 0;
  return v.normalized();
}

bool CPrimCylSector::contains(const Vector3D& vPos) const
{
  if(!m_bReady)
    return false;

  if(vPos.z < m_vCylAxis.orig.z || vPos.z > m_fEndZ)
    return false;

  if(!m_vSectBounds[0].inside(vPos) || !m_vSectBounds[1].inside(vPos))
    return false;

  Vector3D v = vPos - m_vCylAxis.orig;
  v.z = 0;
  double fR = v.length();
  if(fR < m_fSmallR || fR > m_fBigR)
    return false;

  return true;
}

void CPrimCylSector::assign(double x0, double y0, double r)
{
  float fMinZ = FLT_MAX;
  float fMaxZ = -FLT_MAX;

  Vector3F vZ(0, 0, 1), vO(x0, y0, 0), vP;
  Vector3F vP0 = m_vRegNodes.at(0)->pos - vO; // arbitrary point on the cylinder surface.
  vP0.z = 0;
  vP0.normalize();

  float fDotMin = FLT_MAX;
  float fDotMax = -FLT_MAX;
  Vector3F vPmin, vPmax;    // normalized vectors with min and max dot products (vZ, [vP , vP0]), see the code below.

  float fDot;
  CNode3D* pNode = NULL;
  size_t N = m_vRegNodes.size();
  for(size_t i = 0; i < N; i++)
  {
    pNode = m_vRegNodes.at(i);

    vP = pNode->pos - vO;
    vP.z = 0;
    vP.normalize();

    fDot = vZ & (vP * vP0);
    if(fDot < fDotMin)
    {
      fDotMin = fDot;
      vPmin = vP;
    }
    if(fDot > fDotMax)
    {
      fDotMax = fDot;
      vPmax = vP;
    }

    if(pNode->pos.z < fMinZ)
      fMinZ = pNode->pos.z;
    if(pNode->pos.z > fMaxZ)
      fMaxZ = pNode->pos.z;
  }

  m_vCylAxis.dir = Vector3D(0, 0, 1);
  m_vCylAxis.orig = Vector3D(x0, y0, fMinZ);
  m_fEndZ = fMaxZ;
  m_fSmallR = r;
  m_fBigR = m_fSmallR + m_fHeight;

// Defining two restricting planes.
  m_vOrigDir = vPmin;
  m_vEndDir = vPmax;
  m_vSectBounds[0].pos = m_vCylAxis.orig;
  m_vSectBounds[0].norm = vZ * vPmin;
  m_vSectBounds[1].pos = m_vCylAxis.orig;
  m_vSectBounds[1].norm = -vZ * vPmax;
}

/*
void CPrimCylSector::calc_cyl_params()
{
  if(!collect_reg_nodes())
    return;

  calc_sums();

  m_bReady = calc_center_and_radius();
}

void CPrimCylSector::calc_sums()
{
  m_fSX3 = 0;
  m_fSY3 = 0;
  m_fSXY2 = 0;
  m_fSYX2 = 0;
  m_fSX = 0;
  m_fSY = 0;
  m_fSXY = 0;
  m_fSX2 = 0;
  m_fSY2 = 0;

  float fx2, fy2;
  Vector3F vPos;
  size_t nRegNodesCount = m_vRegNodes.size();
  for(size_t i = 0; i < nRegNodesCount; i++)
  {
    vPos = m_vRegNodes.at(i)->pos;

    m_fSX += vPos.x;
    m_fSY += vPos.y;

    m_fSXY += vPos.x * vPos.y;

    fx2 = vPos.x * vPos.x;
    fy2 = vPos.y * vPos.y;

    m_fSX2 += fx2;
    m_fSY2 += fy2;

    m_fSX3 += fx2 * vPos.x;
    m_fSY3 += fy2 * vPos.y;

    m_fSXY2 += vPos.x * fy2;
    m_fSYX2 += vPos.y * fx2;
  }
}

void CPrimCylSector::grad_F(double a, double b, double r, double& dFda, double& dFdb, double& dFdr) const
{
  double a2 = a * a;
  double b2 = b * b;
  double ab = a * b;
  double r2 = r * r;
  size_t N = m_vRegNodes.size();

  dFda = -4 * (m_fSX3 + m_fSXY2) + (12 * m_fSX2 + 4 * m_fSY2) * a - 12 * m_fSX * a2 
        + 4 * N * a * a2 - 4 * N * a * r2 + 8 * m_fSXY * b - 4 * m_fSX * b2 + 4 * m_fSX * r2 
        + 4 * N * a * b2 - 8 * m_fSY * ab;

  dFdb = -4 * (m_fSY3 + m_fSYX2) + (12 * m_fSY2 + 4 * m_fSX2) * b - 12 * m_fSY * b2
        + 4 * N * b * b2 - 4 * N * b * r2 + 8 * m_fSXY * a - 4 * m_fSY * a2 + 4 * m_fSY * r2
        + 4 * N * b * a2 - 8 * m_fSX * ab;

  dFdr = 4 * N * r * r2 - 4 * (m_fSX2 + m_fSY2) * r + 8 * m_fSX * a * r + 8 * m_fSY * b * r
        - 4 * N * a2 * r - 4 * N * b2 * r;
}

bool CPrimCylSector::calc_center_and_radius()
{
  double a1, b1, R1;
  if(!get_initial_guess(a1, b1, R1))
    return false;

// Test output.
  FILE* pStream = open_out_file();
  if(pStream != NULL)
    test_output(pStream, 0, a1, b1, R1, 0);

  double alpha = initial_estimate_alpha(a1, b1, R1);
  double beta = 0.5;

  double dFda, dFdb, dFdr;
  grad_F(a1, b1, R1, dFda, dFdb, dFdr);
  if(criterion(dFda, dFdb, dFdr))
  {
    assign(a1, b1, R1);
    if(pStream != NULL)
      fclose(pStream);

    return true;
  }

  double a2 = a1 - alpha * dFda;
  double b2 = b1 - alpha * dFdb;
  double R2 = R1 - alpha * dFdr;

  double a, b, R;
  UINT nMaxStepsCount = 300, nStep = 1;
  while(nStep < nMaxStepsCount)
  {
    if(pStream != NULL)
      test_output(pStream, nStep, a1, b1, R1, sqrt(dFda * dFda + dFdb * dFdb + dFdr * dFdr));

    grad_F(a2, b2, R2, dFda, dFdb, dFdr);
    if(criterion(dFda, dFdb, dFdr))
      break;

    a = a2 - alpha * dFda + beta * (a2 - a1);
    b = b2 - alpha * dFdb + beta * (b2 - b1);
    R = R2 - alpha * dFdr + beta * (R2 - R1);

    a1 = a2;
    b1 = b2;
    R1 = R2;

    a2 = a;
    b2 = b;
    R2 = R;

    nStep++;
  }

  if(nStep == nMaxStepsCount)
  {
    AfxMessageBox(_T("CPrimCylSector::calc_center_and_radius(): The minimum can not be found by 100 steps."));
  }

  assign(a2, b2, R2);
  if(pStream != NULL)
    fclose(pStream);

  return true;
}

double CPrimCylSector::initial_estimate_alpha(double x, double y, double r) const
{
  double dx = 0.01, dy = 0.01, dr = 0.1;
  double d = sqrt(dx * dx + dy * dy + dr * dr);

  double dFda1, dFdb1, dFdr1, dFda2, dFdb2, dFdr2;
  grad_F(x, y, r, dFda1, dFdb1, dFdr1);
  grad_F(x + dx, y + dy, r + dr, dFda2, dFdb2, dFdr2);

  double dGradFx = dFda2 - dFda1;
  double dGradFy = dFdb2 - dFdb1;
  double dGradFr = dFdr2 - dFdr1;
  double dGradF = sqrt(dGradFx * dGradFx + dGradFy * dGradFy + dGradFr * dGradFr);
  if(dGradF < Const_Almost_Zero)
    dGradF = 0.001;

  return d / dGradF;
}

bool CPrimCylSector::criterion(double dFdx, double dFdy, double dFdr) const
{
  const double fEps = 0.5;
  double fAbsGrad = sqrt(dFdx * dFdx + dFdy * dFdy + dFdr * dFdr);
  return fAbsGrad < fEps;
}

bool CPrimCylSector::get_initial_guess(double& a, double& b, double& r) const
{
  CFace* pFace1 = NULL; // the first face in the collection.
  CFace* pFace2 = NULL; // the face which (x, y) distance from the first face is maximal.
  float fDist, fMaxDist = 0;
  Vector3D v, v1, v2;   // centers of the corresponding faces.
  size_t nRegCount = m_vRegions.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    CRegion* pReg = m_vRegions.at(i);
    size_t nFaceCount = pReg->vFaces.size();
    if(nFaceCount == 0)
      continue;

    if(i == 0)
    {
      pFace1 = pReg->vFaces.at(0);
      v1 = (pFace1->p0->pos + pFace1->p1->pos + pFace1->p2->pos) / 3;
      v1.z = 0;
    }

    CFace* pFace = NULL;
    for(size_t j = 0; j < nFaceCount; j++)
    {
      pFace = pReg->vFaces.at(j);
      v = (pFace->p0->pos + pFace->p1->pos + pFace->p2->pos) / 3;
      v.z = 0;

      fDist = (v - v1).length();
      if(fDist > fMaxDist)
      {
        fMaxDist = fDist;
        pFace2 = pFace;
        v2 = v;
      }
    }
  }

// Intersection point (if exists) of two normals to faces 1 and 2 must give the initial guess
  double fDot = pFace1->norm & pFace2->norm;
  double D = 1 - fDot * fDot;
  if(D < Const_Almost_Zero)
    return false;   // n1 and n2 are parallel, there is no intersection.

  Vector3D dv = v1 - v2;
  double fDot1 = dv & pFace1->norm;
  double fDot2 = dv & pFace2->norm;
  r = (fDot * fDot2 - fDot1) / D;
  Vector3D v3 = v1 + r * pFace1->norm;

  a = v3.x;
  b = v3.y;
  return true;
}

FILE* CPrimCylSector::open_out_file()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  std::string cPath = COutputEngine::get_full_path(pObj->get_filename());
  std::string cName("Calc_Center_and_Radius");
  std::string cExt(".csv");
  std::string cFileName = cPath + cName + cExt;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if((nErr != 0) || (pStream == 0))
    return NULL;

  return pStream;
}

void CPrimCylSector::test_output(FILE* pStream, UINT nStep, double a, double b, double r, double grad)
{
  fprintf_s(pStream, "%d, %10.7lf, %10.7lf, %10.7lf, %10.6lf\n", nStep, a, b, r, grad);
}
*/

//-------------------------------------------------------------------------------------------------
// CEllipticalCylSector - a sectorial volume between two ellipses.
//-------------------------------------------------------------------------------------------------
CEllipticalCylSector::CEllipticalCylSector(const CStringVector& vRegNames, const CEllipseData& data, double fHeight, double fDelta)
  : m_fD(fDelta)
{
  set_height(fHeight);
  if(get_regions(vRegNames) && collect_reg_nodes())
  {
    assign(data);
    m_bReady = prepare_ellipse();
  }
}

Vector2D CEllipticalCylSector::model_coord(const Vector3D& vPos) const
{
  double h = get_h(vPos);
  Vector3D vNorm = get_loc_normal(vPos);

// Find the point on the basic ellipse (a, b) by drawing a normal from vPos onto its surface:
  double x = vPos.x - vNorm.x * h;
  clamp(x);

  int nInterv = 0;
  double len = m_LenSpl.get(x, 0, nInterv);
  return Vector2D(len, h);
}

Vector3D CEllipticalCylSector::get_loc_normal(const Vector3D& vPos) const
{
// The normal vector at a point (x1, y1) of a modified ellipse with semi-axes a' = a + h, b' = b + h:
  double nx = vPos.x - m_vCylAxis.orig.x;

  double h = get_h(vPos);
  double r = m_fA / m_fB;

  double coeff = r * (r + 2 * (1 - r) * h / m_fB);
  double ny = (vPos.y - m_vCylAxis.orig.y) * coeff;

  Vector3D vNorm(nx, ny, 0);
  return vNorm.normalized();
}

double CEllipticalCylSector::get_h(const Vector3D& vPos) const
{
  double fDX = (vPos.x - m_vCylAxis.orig.x) / m_fA;
  double fDY = (vPos.y - m_vCylAxis.orig.y) / m_fB;

  double fDX2 = fDX * fDX;
  double fDY2 = fDY * fDY;

// Approximation of the 2-nd order on h/m_fA << 1 and h/m_fB << 1. 
// Hint: vPos is supposed to lie on another elliptic surface with a' = a + h and b' = b + h. The height h can then be found approximately as
  double h = 0.5 * (fDX2 + fDY2 - 1) / (fDX2 / m_fA + fDY2 / m_fB);

  return h;
}

bool CEllipticalCylSector::contains(const Vector3D& vPos) const
{
  if(!m_bReady)
    return false;

  if(vPos.z < m_vCylAxis.orig.z || vPos.z > m_fEndZ)
    return false;

  if(!m_vSectBounds[0].inside(vPos) || !m_vSectBounds[1].inside(vPos))
    return false;

// Draw a virtual ellipse with semi-axes greater by m_fHeight. If the point is outside this ellipse, return false.
  double a = m_fA + m_fHeight;
  double b = m_fB + m_fHeight;
  if(ellipse(vPos.x, vPos.y, a, b) > 1)
    return false;

  return true;
}

double CEllipticalCylSector::ellipse(double x, double y, double a, double b) const
{
  double dx = (x - m_vCylAxis.orig.x) / a;
  double dy = (y - m_vCylAxis.orig.y) / b;

  return dx * dx + dy * dy;
}

bool CEllipticalCylSector::prepare_ellipse()
{
  double fXmin, fXmax;
  if(!get_x_ranges(fXmin, fXmax))
    return false;

  UINT nSplNodesCount = 11;  // 10 intervals.
  double fDX = (fXmax - fXmin) / (nSplNodesCount - 1);

  m_LenSpl.set_boundary_type(CubicSpline::btFree);
  m_LenSpl.set_nodes_count(nSplNodesCount);
  m_LenSpl.set_range(fXmin, fXmax, true);   // equidistant points.
  m_LenSpl.set_node_val(0, 0.0);

  double fX1 = fXmin, fX2 = fX1 + fDX, fLenX2 = 0;
  for(UINT i = 1; i < nSplNodesCount; i++)
  {
    fLenX2 += fabs(ellipse_arc_length(fX1, fX2));
    m_LenSpl.set_node_val(i, fLenX2);

    fX1 = fX2;
    fX2 = fX1 + fDX;
  }

  m_LenSpl.prepare();   // calculate all the spline coefficients.

// DEBUG output
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  std::string cPath = COutputEngine::get_full_path(pObj->get_filename());
  std::string cName("elliptical_arc_length");
  std::string cExt(".csv");
  std::string cFileName = cPath + cName + cExt;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if((nErr == 0) && (pStream != 0))
  {
    fputs(" x(mm),  length(mm)\n\n", pStream);
    int nInterv = 0;
    UINT nCount = 100;
    double y, x = fXmin, dx = (fXmax - fXmin) / (nCount - 1);
    for(UINT j = 0; j < nCount; j++)
    {
      if(x < m_LenSpl.get_x(0))
        x = m_LenSpl.get_x(0);
      if(x > m_LenSpl.get_x(m_LenSpl.get_nodes_count() - 1))
        x = m_LenSpl.get_x(m_LenSpl.get_nodes_count() - 1);

      y = m_LenSpl.get(x, 0, nInterv);
      fprintf_s(pStream, "%7.3lf, %7.3lf\n", 10*x, 10*y);   // x and y in millimeters.
      x += dx;
    }

    fclose(pStream);
  }

  cName = std::string("two_ellipses");
  cFileName = cPath + cName + cExt;
  nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if((nErr == 0) && (pStream != 0))
  {
    fputs(" x(mm),    y1(mm),    y2(mm),    x3(mm),    y3(mm)\n\n", pStream);
    UINT nCount = 100;
    double y1, y2, x = fXmin, dx = (fXmax - fXmin) / (nCount - 1);
    double dx1, dx2, h;
    Vector3D v, norm;
    for(UINT j = 0; j < nCount; j++)
    {
      
      dx2 = (x - m_vCylAxis.orig.x) / (m_fA + m_fHeight);
      y2 = m_vCylAxis.orig.y + (m_fB + m_fHeight) * sqrt(1.0 - dx2 * dx2);
      v.x = x;
      v.y = y2;
      h = get_h(v);
      norm = get_loc_normal(v);
      v -= (h * norm);

      dx1 = (x - m_vCylAxis.orig.x) / m_fA;
      y1 = m_vCylAxis.orig.y + m_fB * sqrt(1.0 - dx1 * dx1);
      
      fprintf_s(pStream, "%7.3lf, %10.6lf, %10.6lf, %10.6lf, %10.6lf\n", 10*x, 10*y1, 10*y2, 10*v.x, 10*v.y);   // x and y in millimeters.
      x += dx;
    }

    fclose(pStream);
  }

  return true;
}

bool CEllipticalCylSector::get_x_ranges(double& fXmin, double& fXmax) const
{
  fXmin = FLT_MAX;
  fXmax = -FLT_MAX;
  Vector3F vPos;
  size_t N = m_vRegNodes.size();
  for(size_t i = 0; i < N; i++)
  {
    vPos = m_vRegNodes.at(i)->pos;
    if(vPos.x < fXmin)
      fXmin = vPos.x;
    if(vPos.x > fXmax)
      fXmax = vPos.x;
  }

  return (fXmin != FLT_MAX) && (fXmax != -FLT_MAX);
}

double CEllipticalCylSector::ellipse_arc_length(double x1, double x2) const
{
  double fArg1 = (x1 - m_vCylAxis.orig.x) / m_fA;
  double fArg2 = (x2 - m_vCylAxis.orig.x) / m_fA;
  if((fabs(fArg1) > 1) || (fabs(fArg2) > 1))
    return 0;

  double fAlp1 = acos(fArg1);
  double fAlp2 = acos(fArg2);

  double par[2] = { m_fA, m_fB };

  return CMath::integr(fAlp1, fAlp2, &ellipse_elem_arc_length, 50, par);
}

double CEllipticalCylSector::ellipse_elem_arc_length(double alpha, double* par)
{
  double a = par[0];
  double b = par[1];

  double a2 = a * a;
  double b2 = b * b;

  double fCosAlp = cos(alpha);
  double f = a * sqrt(1.0 - fCosAlp * fCosAlp * (1.0 - b2 / a2));
  return f;
}

void CEllipticalCylSector::calc_sums()
{
}

bool CEllipticalCylSector::calc_center_and_radius()
{
/*
  CEllipseData data1;
  if(!get_initial_guess(data1))
    return false;

  collect_opt_nodes(0.1);

  double fEps2 = get_F(data1);

// Test output.
  FILE* pStream = open_out_file();
  if(pStream != NULL)
  {
    fputs("N,  x0 (mm),    y0 (mm),     a (mm),     b (mm),        eps\n", pStream);
    test_output(pStream, 0, data1, fEps2);
  }

  CEllipseData data2;
  UINT nMaxStepsCount = 300, nStep = 1, nBigStep = 0;
  while(nStep < nMaxStepsCount)
  {
    UINT nIterDone;
    if(nStep < 10)
    {
      nIterDone = one_iter_steepest_descent(data1, data2, fEps2);
      if(nIterDone == 0)
        break;
    }
    else
    {
      nIterDone = one_iter_Newton(data1, data2, fEps2);
      if(nIterDone == 0)
      {
        if(nBigStep > 5)
          break;

        nBigStep++;
        nStep = 1;
      }
    }

    if(pStream != NULL)
      test_output(pStream, nStep, data2, fEps2);

    if(criterion(fEps2))
      break;

    data1 = data2;
    nStep++;
  }

  if(nStep == nMaxStepsCount)
  {
    AfxMessageBox(_T("CEllipticalCylSector::calc_center_and_radius(): The minimum can not be found in the course of 300 steps."));
  }

  assign(data2);
  if(pStream != NULL)
    fclose(pStream);
*/
  return true;
}

/*
UINT CEllipticalCylSector::one_iter_Newton(const CEllipseData& start, CEllipseData& end, double& fMinF) const
{
  CGaussEliminationSolver gauss(4, 5);
  if(!grad2_F(start, gauss.get_matrix()))
    return 0;

  int nRes = gauss.solve();
  if(nRes != CGaussEliminationSolver::resOK)
    return 0;

  CEllipseData incr;
// Note that after this call incr will contain the increment of start leading to the minimum.
  get_solution(gauss.get_matrix(), incr);

  end = start + incr;
  return min_along_line(start, end, fMinF);
}

UINT CEllipticalCylSector::one_iter_steepest_descent(const CEllipseData& start, CEllipseData& end, double& fMinF) const
{
  CEllipseData grad;
  grad_F(start, grad);
  double abs_grad = sqrt(grad.x0 * grad.x0 + grad.y0 * grad.y0 + grad.a * grad.a + grad.b * grad.b);
  if(abs_grad < Const_Almost_Zero)
    return 0;

  grad /= abs_grad;   // demensionless direction of the steepest ascent.
  double step = 1e-4; // 5 mcm in cm.
  grad *= step;       // 5 mcm step along the steepest descent.

  UINT nStepsCount = 0;
  double f1 = get_F(start), f2 = -FLT_MAX;
  CEllipseData curr = start, prev = start;
  while(true)
  {
    curr -= grad;
    nStepsCount++;
    f2 = get_F(curr);
    if(f2 > f1)
    {
      end = prev;
      fMinF = f1;
      break;
    }

    f1 = f2;
    prev = curr;
  }

  return nStepsCount;
}

double CEllipticalCylSector::get_F(const CEllipseData& data) const
{
  double x0 = data.x0;
  double y0 = data.y0;
  double a = data.a;
  double b = data.b;

  return get_F(x0, y0, a, b);
}

double CEllipticalCylSector::get_F(double x0, double y0, double a, double b) const
{
  Vector3F vPos;
  double fX, fY, fEps, fFunc = 0;
  size_t N = m_vOptNodes.size();
  for(size_t i = 0; i < N; i++)
  {
    vPos = m_vOptNodes.at(i);
    fX = (vPos.x - x0) / a;
    fY = (vPos.y - y0) / b;

    fEps = fX * fX + fY * fY - 1;
    fFunc += fEps * fEps;
  }

  return fFunc;
}

void CEllipticalCylSector::grad_F(const CEllipseData& data, CEllipseData& grad) const
{
  double x = data.x0;
  double y = data.y0;
  double a = data.a;
  double b = data.b;

  double d = 0.5 / m_fD;

  double fFxp = get_F(x + m_fD, y, a, b);
  double fFxm = get_F(x - m_fD, y, a, b);  
  grad.x0 = d * (fFxp - fFxm);

  double fFyp = get_F(x, y + m_fD, a, b);
  double fFym = get_F(x, y - m_fD, a, b);
  grad.y0 = d * (fFyp - fFym);

  double fFap = get_F(x, y, a + m_fD, b);
  double fFam = get_F(x, y, a - m_fD, b);
  grad.a = d * (fFap - fFam);

  double fFbp = get_F(x, y, a, b + m_fD);
  double fFbm = get_F(x, y, a, b - m_fD);
  grad.b = d * (fFbp - fFbm);
}

bool CEllipticalCylSector::grad2_F(const CEllipseData& data, CMatrix_NxM& mtx) const
{
  UINT Nr = mtx.get_row_count();
  UINT Nc = mtx.get_col_count();
  if((Nr != 4) || (Nc != 5))
    return false;

  double x = data.x0;
  double y = data.y0;
  double a = data.a;
  double b = data.b;

  double d = 0.5 / m_fD;
  double d2 = d * d;

  double fFxp = get_F(x + m_fD, y, a, b);
  double fFxm = get_F(x - m_fD, y, a, b);  
  double dFdx = d * (fFxp - fFxm);
  mtx.set(0, 4, -dFdx);

  double fFyp = get_F(x, y + m_fD, a, b);
  double fFym = get_F(x, y - m_fD, a, b);
  double dFdy = d * (fFyp - fFym);
  mtx.set(1, 4, -dFdy);

  double fFap = get_F(x, y, a + m_fD, b);
  double fFam = get_F(x, y, a - m_fD, b);
  double dFda = d * (fFap - fFam);
  mtx.set(2, 4, -dFda);

  double fFbp = get_F(x, y, a, b + m_fD);
  double fFbm = get_F(x, y, a, b - m_fD);
  double dFdb = d * (fFbp - fFbm);
  mtx.set(3, 4, -dFdb);

  double fF = get_F(x, y, a, b);

  double d2Fdx2 = 4 * d2 * (fFxp - 2 * fF + fFxm);
  mtx.set(0, 0, d2Fdx2);

  double d2Fdxdy = d2 * (get_F(x + m_fD, y + m_fD, a, b) - get_F(x - m_fD, y + m_fD, a, b) - get_F(x + m_fD, y - m_fD, a, b) + get_F(x - m_fD, y - m_fD, a, b));
  mtx.set(0, 1, d2Fdxdy);
  mtx.set(1, 0, d2Fdxdy);

  double d2Fdxda = d2 * (get_F(x + m_fD, y, a + m_fD, b) - get_F(x - m_fD, y, a + m_fD, b) - get_F(x + m_fD, y, a - m_fD, b) + get_F(x - m_fD, y, a - m_fD, b));
  mtx.set(0, 2, d2Fdxda);
  mtx.set(2, 0, d2Fdxda);

  double d2Fdxdb = d2 * (get_F(x + m_fD, y, a, b + m_fD) - get_F(x - m_fD, y, a, b + m_fD) - get_F(x + m_fD, y, a, b - m_fD) + get_F(x - m_fD, y, a, b - m_fD));
  mtx.set(0, 3, d2Fdxdb);
  mtx.set(3, 0, d2Fdxdb);

  double d2Fdy2 = 4 * d2 * (fFyp - 2 * fF + fFym);
  mtx.set(1, 1, d2Fdy2);

  double d2Fdyda = d2 * (get_F(x, y + m_fD, a + m_fD, b) - get_F(x, y - m_fD, a + m_fD, b) - get_F(x, y + m_fD, a - m_fD, b) + get_F(x, y - m_fD, a - m_fD, b));
  mtx.set(1, 2, d2Fdyda);
  mtx.set(2, 1, d2Fdyda);

  double d2Fdydb = d2 * (get_F(x, y + m_fD, a, b + m_fD) - get_F(x, y - m_fD, a, b + m_fD) - get_F(x, y + m_fD, a, b - m_fD) + get_F(x, y - m_fD, a, b - m_fD));
  mtx.set(1, 3, d2Fdydb);
  mtx.set(3, 1, d2Fdydb);

  double d2Fda2 = 4 * d2 * (fFap - 2 * fF + fFam);
  mtx.set(2, 2, d2Fda2);

  double d2Fdadb = d2 * (get_F(x, y, a + m_fD, b + m_fD) - get_F(x, y, a - m_fD, b + m_fD) - get_F(x, y, a + m_fD, b - m_fD) + get_F(x, y, a - m_fD, b - m_fD));
  mtx.set(2, 3, d2Fdadb);
  mtx.set(3, 2, d2Fdadb);

  double d2Fdb2 = 4 * d2 * (fFbp - 2 * fF + fFbm);
  mtx.set(3, 3, d2Fdb2);

  return true;
}

UINT CEllipticalCylSector::min_along_line(const CEllipseData& beg, CEllipseData& end, double& fMinF) const
{
  double f1 = get_F(beg);

  UINT nIterCount = 500;
  CEllipseData d = (end - beg) / nIterCount;

  double f2;
  end = beg;
  CEllipseData prev = beg;
  UINT nIterDone = 0;
  for(UINT i = 0; i < nIterCount; i++)
  {
    end += d;
    f2 = get_F(end);
    fMinF = f2;
    if(f2 > f1)
    {
      end = prev;
      fMinF = f1;
      break;
    }

    nIterDone++;
    prev = end;
    f1 = f2;
  }

  return nIterDone;
}

bool CEllipticalCylSector::get_initial_guess(CEllipseData& data) const
{
  if(!CPrimCylSector::get_initial_guess(data.x0, data.y0, data.a))
    return false;

  data.b = data.a;

  return true;
}

void CEllipticalCylSector::collect_opt_nodes(double fdZ)
{
  m_vOptNodes.clear();

  Vector3F vPos;
  size_t N = m_vRegNodes.size();
  for(size_t i = 0; i < N; i++)
  {
    vPos = m_vRegNodes.at(i)->pos;
    if((vPos.z < -fdZ) || (vPos.z > fdZ))
      continue;

    m_vOptNodes.push_back(m_vRegNodes.at(i)->pos);
  }

  m_vOptNodes.shrink_to_fit();
}

bool CEllipticalCylSector::get_solution(const CMatrix_NxM& mtx, CEllipseData& data) const
{
  if(mtx.get_row_count() != 4 || mtx.get_col_count() != 5)
    return false;

  data.x0 = mtx.get(0, 4);
  data.y0 = mtx.get(1, 4);
  data.a  = mtx.get(2, 4);
  data.b  = mtx.get(3, 4);

  return true;
}

static const double scfEps = 1e-8;

bool CEllipticalCylSector::criterion(double fEps2) const
{
  return fEps2 < scfEps;
}

void CEllipticalCylSector::test_output(FILE* pStream, UINT nStep, const CEllipseData& data, double eps2)
{
  fprintf_s(pStream, "%d, %10.7lf, %10.7lf, %10.7lf, %10.7lf, %10.7lf\n", nStep, 10*data.x0, 10*data.y0, 10*data.a, 10*data.b, eps2);
}
*/
void CEllipticalCylSector::assign(const CEllipseData& data)
{
  CPrimCylSector::assign(data.x0, data.y0, data.a);

  m_fA = data.a;
  m_fB = data.b;
}

};  // namespace EvaporatingParticle
