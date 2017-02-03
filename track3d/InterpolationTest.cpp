
#include "stdafx.h"

#include "InterpolationTest.h"

namespace EvaporatingParticle
{

CInterpolationTest::CInterpolationTest()
  : m_pTetra(NULL), m_pPyramid(NULL), m_pWedge(NULL), m_pHexa(NULL)
{
  set_default();
  create_objects();
  do_test();
}

CInterpolationTest::~CInterpolationTest()
{
  if(m_pTetra != NULL)
    delete m_pTetra;

  if(m_pPyramid != NULL)
    delete m_pPyramid;

  if(m_pWedge != NULL)
    delete m_pWedge;

  if(m_pHexa != NULL)
    delete m_pHexa;
}

void CInterpolationTest::set_default()
{
  Vector3D vX(1, 0, 0), vY(0, 1, 0), vZ(0, 0, 1);

  m_vTetraVert[0] = Vector3D(1, 2, 3);
  m_vTetraVert[1] = m_vTetraVert[0] + 4 * vX + vY;
  m_vTetraVert[2] = m_vTetraVert[1] + 3 * vY - 2 * vX;
  m_vTetraVert[3] = (m_vTetraVert[0] + m_vTetraVert[1] + m_vTetraVert[2]) / 3.0 + 4 * vZ;

  m_vPyramidVert[0] = m_vTetraVert[0];
  m_vPyramidVert[1] = m_vPyramidVert[0] + 4 * vX - vY;
  m_vPyramidVert[2] = m_vPyramidVert[1] + 3 * vY + 2 * vX;
  m_vPyramidVert[3] = m_vPyramidVert[2] + 2 * vY - 3 * vX;
  m_vPyramidVert[4] = (m_vPyramidVert[0] + m_vPyramidVert[1] + m_vPyramidVert[2] + m_vPyramidVert[3]) / 4.0 + 4 * vZ + 2 * vX;

  m_vWedgeVert[0] = m_vTetraVert[0];
  m_vWedgeVert[1] = m_vTetraVert[1];
  m_vWedgeVert[2] = m_vTetraVert[2];
  m_vWedgeVert[3] = m_vWedgeVert[0] + 4 * vZ + 2 * vX;
  m_vWedgeVert[4] = m_vWedgeVert[1] + 4 * vZ + 2 * vX;
  m_vWedgeVert[5] = m_vWedgeVert[2] + 4 * vZ + 2 * vX;

  m_vHexaVert[0] = m_vPyramidVert[0];
  m_vHexaVert[1] = m_vPyramidVert[1];
  m_vHexaVert[2] = m_vPyramidVert[3];
  m_vHexaVert[3] = m_vPyramidVert[2];
  m_vHexaVert[4] = m_vHexaVert[0] + 4 * vZ;
  m_vHexaVert[5] = m_vHexaVert[1] + 4 * vZ;
  m_vHexaVert[6] = m_vHexaVert[2] + 4 * vZ;
  m_vHexaVert[7] = m_vHexaVert[3] + 4 * vZ;
}

void CInterpolationTest::do_test()
{
  test_func();
  test_param_at_nodes();
  test_inverse();

  test_inside_outside();
}

void CInterpolationTest::create_objects()
{
  CNode3D* p0 = new CNode3D(m_vTetraVert[0]);
  CNode3D* p1 = new CNode3D(m_vTetraVert[1]);
  CNode3D* p2 = new CNode3D(m_vTetraVert[2]);
  CNode3D* p3 = new CNode3D(m_vTetraVert[3]);
  m_pTetra = new CTetra(p0, p1, p2, p3);

  p0 = new CNode3D(m_vPyramidVert[0]);
  p1 = new CNode3D(m_vPyramidVert[1]);
  p2 = new CNode3D(m_vPyramidVert[2]);
  p3 = new CNode3D(m_vPyramidVert[3]);
  CNode3D* p4 = new CNode3D(m_vPyramidVert[4]);
  m_pPyramid = new CPyramid(p0, p1, p2, p3, p4);

  p0 = new CNode3D(m_vWedgeVert[0]);
  p1 = new CNode3D(m_vWedgeVert[1]);
  p2 = new CNode3D(m_vWedgeVert[2]);
  p3 = new CNode3D(m_vWedgeVert[3]);
  p4 = new CNode3D(m_vWedgeVert[4]);
  CNode3D* p5 = new CNode3D(m_vWedgeVert[5]);
  m_pWedge = new CWedge(p0, p1, p2, p3, p4, p5);

  p0 = new CNode3D(m_vHexaVert[0]);
  p1 = new CNode3D(m_vHexaVert[1]);
  p2 = new CNode3D(m_vHexaVert[2]);
  p3 = new CNode3D(m_vHexaVert[3]);
  p4 = new CNode3D(m_vHexaVert[4]);
  p5 = new CNode3D(m_vHexaVert[5]);
  CNode3D* p6 = new CNode3D(m_vHexaVert[6]);
  CNode3D* p7 = new CNode3D(m_vHexaVert[7]);
  m_pHexa = new CHexa(p0, p1, p2, p3, p4, p5, p6, p7);
}

void CInterpolationTest::test_func()
{
// All vRes must be exactly zero.
  Vector3D vRes = m_pTetra->func(0, 0, 0, m_vTetraVert[0]);
  vRes = m_pTetra->func(1, 0, 0, m_vTetraVert[1]);
  vRes = m_pTetra->func(0, 1, 0, m_vTetraVert[2]);
  vRes = m_pTetra->func(0, 0, 1, m_vTetraVert[3]);

  vRes = m_pPyramid->func(0, 0, 0, m_vPyramidVert[0]);
  vRes = m_pPyramid->func(1, 0, 0, m_vPyramidVert[1]);
  vRes = m_pPyramid->func(0, 1, 0, m_vPyramidVert[3]);
  vRes = m_pPyramid->func(1, 1, 0, m_vPyramidVert[2]);
  vRes = m_pPyramid->func(0, 0, 1, m_vPyramidVert[4]);

  vRes = m_pWedge->func(0, 0, 0, m_vWedgeVert[0]);
  vRes = m_pWedge->func(1, 0, 0, m_vWedgeVert[1]);
  vRes = m_pWedge->func(0, 1, 0, m_vWedgeVert[2]);
  vRes = m_pWedge->func(1, 0, 1, m_vWedgeVert[4]);
  vRes = m_pWedge->func(0, 1, 1, m_vWedgeVert[5]);
  vRes = m_pWedge->func(0, 0, 1, m_vWedgeVert[3]);

  vRes = m_pHexa->func(0, 0, 0, m_vHexaVert[0]);
  vRes = m_pHexa->func(1, 0, 0, m_vHexaVert[1]);
  vRes = m_pHexa->func(0, 1, 0, m_vHexaVert[2]);
  vRes = m_pHexa->func(1, 1, 0, m_vHexaVert[3]);
  vRes = m_pHexa->func(0, 0, 1, m_vHexaVert[4]);
  vRes = m_pHexa->func(1, 0, 1, m_vHexaVert[5]);
  vRes = m_pHexa->func(0, 1, 1, m_vHexaVert[6]);
  vRes = m_pHexa->func(1, 1, 1, m_vHexaVert[7]);
}

void CInterpolationTest::test_param_at_nodes()
{
  double s, t, u;
  m_pTetra->param(m_vTetraVert[0], s, t, u);  // (0, 0, 0)
  m_pTetra->param(m_vTetraVert[1], s, t, u);  // (1, 0, 0)
  m_pTetra->param(m_vTetraVert[2], s, t, u);  // (0, 1, 0)
  m_pTetra->param(m_vTetraVert[3], s, t, u);  // (0, 0, 1)

  m_pPyramid->param(m_vPyramidVert[0], s, t, u);  // (0, 0, 0)
  m_pPyramid->param(m_vPyramidVert[1], s, t, u);  // (1, 0, 0)
  m_pPyramid->param(m_vPyramidVert[2], s, t, u);  // (1, 1, 0)
  m_pPyramid->param(m_vPyramidVert[3], s, t, u);  // (0, 1, 0)
  m_pPyramid->param(m_vPyramidVert[4], s, t, u);  // (0, 0, 1)

  m_pWedge->param(m_vWedgeVert[0], s, t, u);  // (0, 0, 0)
  m_pWedge->param(m_vWedgeVert[1], s, t, u);  // (1, 0, 0)
  m_pWedge->param(m_vWedgeVert[2], s, t, u);  // (0, 1, 0)
  m_pWedge->param(m_vWedgeVert[3], s, t, u);  // (0, 0, 1)
  m_pWedge->param(m_vWedgeVert[4], s, t, u);  // (1, 0, 1)
  m_pWedge->param(m_vWedgeVert[5], s, t, u);  // (0, 1, 1)

  m_pHexa->param(m_vHexaVert[0], s, t, u);  // (0, 0, 0)
  m_pHexa->param(m_vHexaVert[1], s, t, u);  // (1, 0, 0)
  m_pHexa->param(m_vHexaVert[2], s, t, u);  // (0, 1, 0)
  m_pHexa->param(m_vHexaVert[3], s, t, u);  // (1, 1, 0)
  m_pHexa->param(m_vHexaVert[4], s, t, u);  // (0, 0, 1)
  m_pHexa->param(m_vHexaVert[5], s, t, u);  // (1, 0, 1)
  m_pHexa->param(m_vHexaVert[6], s, t, u);  // (0, 1, 1)
  m_pHexa->param(m_vHexaVert[7], s, t, u);  // (1, 1, 1)
}

void CInterpolationTest::apply_rotation(CElem3D* pElem, const Matrix3D& m)
{
  if(pElem == NULL)
    return;

  size_t nNodeCount = pElem->vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    CNode3D* pNode = pElem->vNodes.at(i);
    pNode->pos = m * pNode->pos;
  }

  pElem->init();
}

void CInterpolationTest::test_inverse()
{
  UINT i;
// Tetrahedron
  double pTetWeight[4] = { 2, 3, 3, 2 };
  double fSumWeight = 0;
  for(i = 0; i < 4; i++)
    fSumWeight += pTetWeight[i];

  Vector3D pos;
  for(i = 0; i < 4; i++)
  {
    pTetWeight[i] /= fSumWeight;
    pos += pTetWeight[i] * m_vTetraVert[i];     // must be inside the tetrahedron.
  }

  double s, t, u;
  m_pTetra->param(pos, s, t, u);
  Vector3D vRes = m_pTetra->func(s, t, u, pos);   // must be zero.

  Vector3D vDir(1, 2, 3);
  Matrix3D mRot = Matrix3D::rot(vDir, 170., false);
// A small test for Matrix3D:
  double fDet = mRot.det();             // must be unity.
  Matrix3D mRotInv = mRot.transposed();
  Matrix3D nUnity = mRot * mRotInv;     // must be unit matrix.

// Rotate the coordinate system...
  double srot, trot, urot;
  pos = mRot * pos;
  apply_rotation(m_pTetra, mRot);

// ... and again, find parametric coordinated of pos.
  //Commented AC 02/07/2016
  //CElem3D* pElem = m_pTetra->inside(pos);       // must be inside.
  m_pTetra->param(pos, srot, trot, urot);       // must be srot == s, trot == t, urot == u.
  vRes = m_pTetra->func(srot, trot, urot, pos); // must be zero.

// Pyramid
  double pPyrWeight[5] = { 1, 1, 2, 2, 5 };
  fSumWeight = 0;
  for(i = 0; i < 5; i++)
    fSumWeight += pPyrWeight[i];

  pos = Vector3D();
  for(i = 0; i < 5; i++)
  {
    pPyrWeight[i] /= fSumWeight;
    pos += pPyrWeight[i] * m_vPyramidVert[i];     // must be inside the pyramid.
  }

  m_pPyramid->param(pos, s, t, u);
  vRes = m_pPyramid->func(s, t, u, pos);   // must be zero.

// Rotate the coordinate system...
  pos = mRot * pos;
  apply_rotation(m_pPyramid, mRot);

// ... and again, find parametric coordinated of pos.
//Commented AC 02/07/2016
  //pElem = m_pPyramid->inside(pos);       // must be inside.
  m_pPyramid->param(pos, srot, trot, urot);       // must be srot == s, trot == t, urot == u.
  vRes = m_pPyramid->func(srot, trot, urot, pos); // must be zero.

// Wedge
  double pWedgeWeight[6] = { 6, 5, 4, 3, 2, 1 };
  fSumWeight = 0;
  for(i = 0; i < 6; i++)
    fSumWeight += pWedgeWeight[i];

  pos = Vector3D();
  for(i = 0; i < 6; i++)
  {
    pWedgeWeight[i] /= fSumWeight;
    pos += pWedgeWeight[i] * m_vWedgeVert[i];     // must be inside the pyramid.
  }

  m_pWedge->param(pos, s, t, u);
  vRes = m_pWedge->func(s, t, u, pos);   // must be zero.

// Rotate the coordinate system...
  pos = mRot * pos;
  apply_rotation(m_pWedge, mRot);

// ... and again, find parametric coordinated of pos.
//Commented AC 02/07/2016
  //pElem = m_pWedge->inside(pos);       // must be inside.
  m_pWedge->param(pos, srot, trot, urot);       // must be srot == s, trot == t, urot == u.
  vRes = m_pWedge->func(srot, trot, urot, pos); // must be zero.

// Hexahedron
  double pHexaWeight[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
  fSumWeight = 0;
  for(i = 0; i < 8; i++)
    fSumWeight += pHexaWeight[i];

  pos = Vector3D();
  for(i = 0; i < 8; i++)
  {
    pHexaWeight[i] /= fSumWeight;
    pos += pHexaWeight[i] * m_vHexaVert[i];     // must be inside the pyramid.
  }

  m_pHexa->param(pos, s, t, u);
  vRes = m_pHexa->func(s, t, u, pos);   // must be zero.

// Rotate the coordinate system...
  pos = mRot * pos;
  apply_rotation(m_pHexa, mRot);

// ... and again, find parametric coordinated of pos.
//Commented AC 02/07/2016
  //pElem = m_pHexa->inside(pos);       // must be inside.
  m_pHexa->param(pos, srot, trot, urot);       // must be srot == s, trot == t, urot == u.
  vRes = m_pHexa->func(srot, trot, urot, pos); // must be zero.
}

void CInterpolationTest::test_inside_outside()
{
  double fEps = 0.001;
// Tetrahedron.
  Vector3D vC;
  for(UINT i = 0; i < 4; i++)
    vC += 0.25 * m_vTetraVert[i];     // must be inside the tetrahedron.

  Vector3D vDirToInside = (vC - m_vTetraVert[0]).normalized();
  Vector3D vInside = m_vTetraVert[0] + fEps * vDirToInside;
  //Commented AC 02/07/2016
  //CElem3D* pElem = m_pTetra->inside(vInside);

  Vector3D vOutside = m_vTetraVert[0] - fEps * vDirToInside;
  //Commented AC 02/07/2016
  //pElem = m_pTetra->inside(vOutside);

// Pyramid.
  vC = Vector3D();
  for(UINT i = 0; i < 5; i++)
    vC += 0.2 * m_vPyramidVert[i];     // must be inside the pyramid.

  vDirToInside = (vC - m_vPyramidVert[4]).normalized();
  vInside = m_vPyramidVert[4] + fEps * vDirToInside;
  //Commented AC 02/07/2016
  //pElem = m_pPyramid->inside(vInside);

  vOutside = m_vPyramidVert[4] - fEps * vDirToInside;
  //Commented AC 02/07/2016
  //pElem = m_pPyramid->inside(vOutside);

// Wedge.
  vC = Vector3D();
  for(UINT i = 0; i < 6; i++)
    vC += 0.16666667 * m_vWedgeVert[i]; // must be inside the wedge.

  vDirToInside = (vC - m_vWedgeVert[3]).normalized();
  vInside = m_vWedgeVert[3] + fEps * vDirToInside;
  //Commented AC 02/07/2016
  //pElem = m_pWedge->inside(vInside);

  vOutside = m_vWedgeVert[3] - fEps * vDirToInside;
  //Commented AC 02/07/2016
  //pElem = m_pWedge->inside(vOutside);

// Hexahedron
  vC = Vector3D();
  for(UINT i = 0; i < 8; i++)
    vC += 0.125 * m_vHexaVert[i]; // must be inside the wedge.

  vDirToInside = (vC - m_vHexaVert[7]).normalized();
  vInside = m_vHexaVert[7] + fEps * vDirToInside;
  //Commented AC 02/07/2016
  //pElem = m_pHexa->inside(vInside);

  vOutside = m_vHexaVert[7] - fEps * vDirToInside;
  //Commented AC 02/07/2016
  //pElem = m_pHexa->inside(vOutside);
}

//static CInterpolationTest sTest;

}; // namespace EvaporatingParticle