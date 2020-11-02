
#include "stdafx.h"

#include "DrawTrack.h"
#include "Tracker.hpp"

#include "float.h"
#include "constant.hpp"

#include "ExecutionDialog.h"
#include "ParticleTracking.h"
#include "MainFrm.h"

#include <sstream>
#include <algorithm>

#include "DirichletTesselation.h"
#include "TrajectSelector.h"
#include "Primitives.h"

namespace EvaporatingParticle
{

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
CTrackDraw::CTrackDraw()
  : m_pTracker(NULL), m_hWnd(NULL), m_hDC(NULL), m_hRC(NULL), m_Color(clLtGray), m_pRegUnderCursor(NULL)
{
// User-defined flags:
  m_bDrawTracks = true;
  m_nDrawMode = dmFlatAndWire;
  m_bDrawNorm = false;
  m_bDrawSelFaces = true;   // enable/disable selected faces drawing.

  m_bWireframeReady = false;
  m_bFacesReady = false;
  m_bNormReady = false;
  m_bAuxReady = false;

// Material properties:
  m_fMatAmbient = 0.3f;
  m_fMatDiffuse = 0.2f;
  m_fMatSpec = 0.3f;
  m_fMatShininess = 0.3f;

// Run-time:
  m_bNewData = true;
  m_bSelRegFlag = false;
  m_bSelTrajectFlag = false;
  m_bSelFacesFlag = false;    // this flag shows whether the faces selection context is ON or OFF.

  m_nTrajUnderCursorId = -1;

  m_fScale = 1.;
  m_vShift = Vector2D(0, 0);
  m_nRegime = nRegimeNone;
  m_fOpacity = 1.0f;

  set_bkgr_color(clLtGray);
  m_bRotCenter = false;

  m_nOvrAxis = 0;
  m_bBusy = false;
}

CTrackDraw::~CTrackDraw()
{
  clear_clr_regions();

  HGLRC hRC = wglGetCurrentContext();
  if(m_hRC == hRC)
    wglMakeCurrent(NULL, NULL);

  if(m_hRC)
    wglDeleteContext(m_hRC);

  if(m_hDC)
    ReleaseDC(m_hWnd, m_hDC);
}

void CTrackDraw::clear()
{
  clear_clr_regions();

  m_vFacesSelRegionVert.clear();
  m_vSelFacesVert.clear();

  m_vWireFrame.clear();
  m_vAuxLines.clear();

  m_vCrossSectVert.clear();
  m_vSelRegVert.clear();
}

void CTrackDraw::clear_clr_regions()
{
  size_t nClrCount = m_vClrFaceVert.size();
  for(size_t i = 0; i < nClrCount; i++)
    m_vClrFaceVert.at(i).clear();

  m_vClrFaceVert.clear();
}

void CTrackDraw::set_window_handle(HWND hwnd)
{
  m_hWnd = hwnd;
  m_hDC = GetDC(hwnd);

  create_window_layer();
}

void CTrackDraw::set_tracker(CTracker* pTr)
{
  if(m_pTracker == NULL || !m_pTracker->is_ready())
  {
    m_pTracker = pTr;
    if(m_bNewData)
      set_view();

    set_hidden_reg_names();
  }
}

bool CTrackDraw::create_window_layer()
{
  PIXELFORMATDESCRIPTOR pfd ;
  memset(&pfd, 0, sizeof(PIXELFORMATDESCRIPTOR)) ;
  pfd.nSize = sizeof(PIXELFORMATDESCRIPTOR);  
  pfd.nVersion = 1;
  pfd.dwFlags = PFD_DOUBLEBUFFER | PFD_SUPPORT_OPENGL | PFD_DRAW_TO_WINDOW;
  pfd.iPixelType = PFD_TYPE_RGBA;

  pfd.cColorBits = 24;             // 24 of bits for color
  pfd.cDepthBits = 32;             // 32-bit depth buffer

// Choose the pixel format:
  int nPixelFormat = ChoosePixelFormat(m_hDC, &pfd);
  if(nPixelFormat == 0)
  {
    TRACE("ChoosePixelFormat Failed %d\r\n", GetLastError()) ;
    return false ;
  }
  TRACE("Pixel Format %d\r\n", nPixelFormat) ;

// Set the pixel format:
  BOOL bResult = SetPixelFormat(m_hDC, nPixelFormat, &pfd);
  if(!bResult)
  {
    TRACE("SetPixelFormat Failed %d\r\n", GetLastError()) ;
    return false ;
  }
  
// Create a rendering context:
  m_hRC = wglCreateContext(m_hDC);
  if(!m_hRC)
  {
    TRACE("wglCreateContext Failed %x\r\n", GetLastError()) ;
    return false;
  }

// Make the rendering context the calling thread's current rendering context:
  BOOL bOK = wglMakeCurrent(m_hDC, m_hRC);
  if(!bOK)
  {
    TRACE("wglMakeCurrent Failed %x\r\n", GetLastError()) ;
    return false;
  }

  return true;
}

void CTrackDraw::set_view(int nViewDir)
{
  if(m_pTracker == NULL)
    return;

  long nResX, nResY;
  get_resolution(nResX, nResY);
  if(nResX < 1)
    return;

  double fAspect = (double)nResY / nResX;

  const CBox& box = m_pTracker->get_box();
  Vector3D vC = box.get_center();
  m_vCenter = vC;

  double fSizeX, fSizeY, fRotX, fRotY, fRotZ;
  switch(nViewDir)
  {
    case dirNone:
    {
      fSizeX = box.vMax.x - box.vMin.x;
      fSizeY = box.vMax.y - box.vMin.y;
      fRotX = 20;
      fRotY = -150;
      fRotZ = 0;
      break;
    }
    case dirPlusX:
    {
      fSizeX = box.vMax.z - box.vMin.z;
      fSizeY = box.vMax.y - box.vMin.y;
      fRotX = 0;
      fRotY = -90;
      fRotZ = 0;
      break;
    }
    case dirMinusX:
    {
      fSizeX = box.vMax.z - box.vMin.z;
      fSizeY = box.vMax.y - box.vMin.y;
      fRotX = 0;
      fRotY = 90;
      fRotZ = 0;
      break;
    }
    case dirPlusY:
    {
      fSizeX = box.vMax.x - box.vMin.x;
      fSizeY = box.vMax.z - box.vMin.z;
      fRotX = 90;
      fRotY = 180;
      fRotZ = 0;
      break;
    }
    case dirMinusY:
    {
      fSizeX = box.vMax.x - box.vMin.x;
      fSizeY = box.vMax.z - box.vMin.z;
      fRotX = -90;
      fRotY = 0;
      fRotZ = 0;
      break;
    }
    case dirPlusZ:
    {
      fSizeX = box.vMax.x - box.vMin.x;
      fSizeY = box.vMax.y - box.vMin.y;
      fRotX = 0;
      fRotY = 0;
      fRotZ = 0;
      break;
    }
    case dirMinusZ:
    {
      fSizeX = box.vMax.x - box.vMin.x;
      fSizeY = box.vMax.y - box.vMin.y;
      fRotX = 0;
      fRotY = 180;
      fRotZ = 0;
      break;
    }
  }

  if((fSizeY < fSizeX) && ((fSizeY / fSizeX) < fAspect))
    m_fScrWidth = fSizeX;
  else
    m_fScrWidth = fSizeY / fAspect;

  m_nRegime = nRegimeNone;
  m_vShift = Vector2D(0, 0);
  m_fRotAngleX = 0;
  m_fRotAngleY = 0;
  m_fScale = 1;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glTranslated(vC.x, vC.y, vC.z);
  glRotated(fRotX, 1., 0., 0.);
  glRotated(fRotY, 0., 1., 0.);
  glRotated(fRotZ, 0., 0., 1.);
  glTranslated(-vC.x, -vC.y, -vC.z);
}

void CTrackDraw::get_resolution(long& nx, long& ny) const
{
  RECT rect;
  if(m_hWnd == NULL || !GetClientRect(m_hWnd, &rect))
    return;

  nx = rect.right - rect.left;
  ny = rect.bottom - rect.top;
}

void CTrackDraw::set_projection()
{
// Dimensions of the viewport.
  long nResX, nResY;
  get_resolution(nResX, nResY);
  if(nResX < 1)
    return;

  glViewport(0, 0, nResX, nResY);

  double fAspect = (double)nResY / nResX;

// Dimensions of the 3D object.
  Vector3D vC = m_pTracker->get_center();

  m_fScrWidth *= m_fScale;
  m_fPix2cm = m_fScrWidth / nResX;

  double fHalfW = 0.5 * m_fScrWidth;
  double fMinX = vC.x - fHalfW;
  double fMaxX = vC.x + fHalfW;

  double fHalfH = fAspect * fHalfW;
  double fMinY = vC.y - fHalfH;
  double fMaxY = vC.y + fHalfH;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(fMinX - m_vShift.x, fMaxX - m_vShift.x, fMinY - m_vShift.y, fMaxY - m_vShift.y, vC.z - 10000, vC.z + 10000);

  glMatrixMode(GL_MODELVIEW);

// Visualization of the Dirichlet cells support.
  vC = Vector3D(0, 0, 0);
  if(m_bRotCenter)  // && m_bDrawNorm)
  {
    vC = m_vCenter;
    glTranslated(vC.x, vC.y, vC.z);
  }

  glRotated(m_fRotAngleX, 1., 0., 0.);
  glRotated(m_fRotAngleY, 0., 1., 0.);
  glRotated(m_fRotAngleZ, 0., 0., 1.);

  glTranslated(-vC.x, -vC.y, -vC.z);

  m_fRotAngleX = 0;
  m_fRotAngleY = 0;
  m_fRotAngleZ = 0;
  m_fScale = 1;  
}

void CTrackDraw::build_arrays()
{
  if(!m_pTracker->is_ready())
    return; // nothing to draw.

  build_wireframe_array();

  build_norm_array();

  build_faces_array();
  build_aux_arrays();

  if(m_bNewData)
  {
    set_view();
    m_bNewData = false;
  }
}

void CTrackDraw::build_faces_array()
{
  if(m_bFacesReady)
    return;

  set_job_name("Building faces");
  CObject::set_progress(0);

  glFrontFace(GL_CCW);
  glDisable(GL_NORMALIZE);
  glDisable(GL_CULL_FACE);

  CFace* pFace = NULL;
  CRegion* pReg = NULL;

  clear_clr_regions();
  m_vFaceCrossSectVert.clear();

  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();

  CSelAreasColl* pSelAreasColl = CParticleTrackingApp::Get()->GetSelAreas();
  int nAreasCount = pSelAreasColl->size();
  for(int k = -1; k < nAreasCount; k++)
  {
    CSelectedAreas* pArea = (k >= 0) ? pSelAreasColl->at(k) : pSelAreasColl->get_default_area();
    size_t nCount = pArea->size();
    CColoredRegVertices vFaceVert;
    vFaceVert.set_reg_color(pArea->get_faces_color());
    for(size_t l = 0; l < nCount; l++)
    {
      std::string sRegName = pArea->at(l);
      int nRegId = CAnsysMesh::get_region_id(sRegName);
      if(nRegId < 0)
        continue;

      pReg = regions.at(nRegId);
      if(!pReg->bEnabled || pReg->bSelected || pReg->bCrossSection)
        continue;

      const CFacesCollection& faces = pReg->vFaces;
      size_t nFaceCount = faces.size();
      for(size_t j = 0; j < nFaceCount; j++)
      {
        pFace = faces.at(j);
        vFaceVert.push_back(CFaceVertex(pFace->p0->pos, pFace->norm));
        vFaceVert.push_back(CFaceVertex(pFace->p1->pos, pFace->norm));
        vFaceVert.push_back(CFaceVertex(pFace->p2->pos, pFace->norm));
      }
    }

    m_vClrFaceVert.push_back(vFaceVert);
  }

  build_cs_faces_array();

  m_bFacesReady = true;
}

void CTrackDraw::build_cs_faces_array()
{
  CFace* pFace = NULL;
  CRegion* pReg = NULL;

  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    pReg = regions.at(i);
    if(!pReg->bCrossSection || !pReg->bEnabled || pReg->bSelected)
      continue;

    const CFacesCollection& faces = pReg->vFaces;
    size_t nFaceCount = faces.size();

    for(size_t j = 0; j < nFaceCount; j++)
    {
      pFace = faces.at(j);
      m_vFaceCrossSectVert.push_back(CFaceVertex(pFace->p0->pos, pFace->norm));
      m_vFaceCrossSectVert.push_back(CFaceVertex(pFace->p1->pos, pFace->norm));
      m_vFaceCrossSectVert.push_back(CFaceVertex(pFace->p2->pos, pFace->norm));
    }
  }
}

void CTrackDraw::build_aux_arrays()
{
  if(m_bAuxReady)
    return;

  build_sel_regions_array();
  build_sel_faces_array();

  m_bAuxReady = true;
}

void CTrackDraw::build_sel_faces_array()
{
  CFace* pFace = NULL;
  CRegion* pReg = NULL;

  m_vSelFacesVert.clear();
  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();

  size_t nSelFacesCount = m_vSelFaces.size();
  for(size_t i = 0; i < nSelFacesCount; i++)
  {
    const CRegFacePair& face = m_vSelFaces.at(i);
    if(face.nReg >= nRegCount)
      continue;

    pReg = regions.at(face.nReg);
    if(face.nFace >= pReg->vFaces.size())
      continue;

    pFace = pReg->vFaces.at(face.nFace);
    m_vSelFacesVert.push_back(CFaceVertex(pFace->p0->pos, pFace->norm));
    m_vSelFacesVert.push_back(CFaceVertex(pFace->p1->pos, pFace->norm));
    m_vSelFacesVert.push_back(CFaceVertex(pFace->p2->pos, pFace->norm));
  }
}

void CTrackDraw::build_sel_regions_array()
{
  CFace* pFace = NULL;
  CRegion* pReg = NULL;

  m_vFacesSelRegionVert.clear();
  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    pReg = regions.at(i);
    if(!pReg->bEnabled || !pReg->bSelected)
      continue;

    const CFacesCollection& faces = pReg->vFaces;
    size_t nFaceCount = faces.size();
    for(size_t j = 0; j < nFaceCount; j++)
    {
      pFace = faces.at(j);
      m_vFacesSelRegionVert.push_back(CFaceVertex(pFace->p0->pos, pFace->norm));
      m_vFacesSelRegionVert.push_back(CFaceVertex(pFace->p1->pos, pFace->norm));
      m_vFacesSelRegionVert.push_back(CFaceVertex(pFace->p2->pos, pFace->norm));
    }
  }
}

void CTrackDraw::build_norm_array()
{
  if(m_bNormReady)
    return;

  m_vAuxLines.clear();
  if(!m_bDrawNorm)
    return;

// Temporarily interpret m_nDrawnCell as index of a region in the regions collection.
//  const CRegionsCollection& vRegs = m_pTracker->get_regions();
//  if(m_nDrawnCell >= vRegs.size())
//    return;

  const CNodesVector& vNodes = m_pTracker->get_nodes();
  size_t nNodeCount = vNodes.size();
  if(m_nDrawnCell >= nNodeCount)
    return;

//  CRegion* pReg = vRegs.at(m_nDrawnCell);
//  CNode3D* pTestNode = pReg->vFaces.at(pReg->vFaces.size() / 2)->p0;

  const CNode3D& TestNode = vNodes.at(m_nDrawnCell);
  if(TestNode.vNbrFaces.size() != 0)
    return;   // boundary cells are not built as 3D objects.

  m_vCenter = TestNode.pos;

  CDirichletTesselation test_obj(true);
  test_obj.set_mesh((CAnsysMesh*)m_pTracker);
  CDirichletCell* pCell = test_obj.build_cell_in_node(TestNode, vNodes);

// DEBUG - visualization of normals to the elliptical surface.
/*
  m_vAuxLines.clear();
  CFieldPtbCollection& vPtbColl = m_pTracker->get_field_ptb();
  size_t nPtbCount = vPtbColl.size();
  for(size_t i = 0; i < nPtbCount; i++)
  {
    CFieldPerturbation* pPtb = vPtbColl.at(i);
    if(pPtb->type() != CFieldPerturbation::ptbElliptSubstrRF)
      continue;

    Vector3F vPos, vNorm;
    CEllipticalSubstrateRF* pElliptPtb = (CEllipticalSubstrateRF*)pPtb;
    CEllipticalCylSector* pEllipse = (CEllipticalCylSector*)(pElliptPtb->get_bound_shape());
    CNodesCollection* pNodesColl = pEllipse->get_reg_nodes();
    size_t nNodesCount = pNodesColl->size();
    for(size_t j = 0; j < nNodesCount; j++)
    {
      vPos = pNodesColl->at(j)->pos;
      if(vPos.z > 0.01 || vPos.z < -0.01)
        continue;

      vNorm = pEllipse->get_loc_normal(vPos);

      m_vAuxLines.push_back(CEdgeVertex(vPos));
      m_vAuxLines.push_back(CEdgeVertex(vPos + vNorm * 0.1));
    }
  }
*/
// END DEBUG
  m_bNormReady = true;

/*
//  if(m_bNormReady)
//    return;

//  m_vAuxLines.clear();
//  if(!m_bDrawNorm)
//    return;

  const CRegionsCollection& vRegs = m_pTracker->get_regions();
//  const CNodesCollection& vNodes = m_pTracker->get_nodes();
//  size_t nNodeCount = vNodes.size();
  size_t nNbrCount, nReg, nFace;
  CNode3D* pNode = NULL;
  Vector3D vNorm;
  for(size_t i = 0; i < nNodeCount; i++)
  {
    pNode = vNodes.at(i);
    nNbrCount = pNode->vNbrFaces.size();
    if(nNbrCount == 0)
      continue;   // inner vertex.

    bool bTrivial = true;
    vNorm = Vector3D(0, 0, 0);
    for(size_t j = 0; j < nNbrCount; j++)
    {
      nReg = pNode->vNbrFaces.at(j).nReg;
      if(nReg >= vRegs.size() || !vRegs.at(nReg)->bEnabled) // do not show normals at disabled regions; (this, however, will distort the
        continue;                                           // normals at the bounary vertices, belonging to both enabled and disabled regions).

      nFace = pNode->vNbrFaces.at(j).nFace;
      if(nFace >= vRegs.at(nReg)->vFaces.size())
        continue;

      vNorm += vRegs.at(nReg)->vFaces.at(nFace)->norm;
      bTrivial = false;
    }

    if(bTrivial)
      continue;

    vNorm.normalize();
    m_vAuxLines.push_back(CEdgeVertex(pNode->pos));
    m_vAuxLines.push_back(CEdgeVertex(pNode->pos + 0.05 * vNorm));
  }

  m_bNormReady = true;
*/
}

void CTrackDraw::build_wireframe_array()
{
  if(m_bWireframeReady)
    return;

  m_vWireFrame.clear();

  set_job_name("Building regions from planes");
  CObject::set_progress(0);
  const CRegionsCollection& regions = m_pTracker->get_regions();

  CFace* pFace = NULL;
  CRegion* pReg = NULL;

  set_job_name("Building wireframe");
  CObject::set_progress(0);

  size_t nRegCount = regions.size();
  size_t nTotalCount = 0;
  for(size_t i = 0; i < nRegCount; i++)
  {
    pReg = regions.at(i);
    nTotalCount += pReg->bCrossSection ? 0 : pReg->vFaces.size(); // for set_progress() only.
  }

  size_t nProcessedCount = 0;
  for(size_t i = 0; i < nRegCount; i++)
  {
    pReg = regions.at(i);
    if(pReg->bCrossSection)
      continue;

    CEdgesVector vEdges;
    std::vector<CEdge>::iterator iter;

    const CFacesCollection& faces = pReg->vFaces;
    size_t nFaceCount = faces.size();
    for(size_t j = 0; j < nFaceCount; j++)
    {
      pFace = faces.at(j);
      if(!pFace->is_wireframe_face())
        continue;

      CEdge edges[] = { CEdge(pFace->p0, pFace->p1, pFace), CEdge(pFace->p1, pFace->p2, pFace), CEdge(pFace->p2, pFace->p0, pFace) };
      for(UINT k = 0; k < 3; k++)
      {
        iter = std::find(vEdges.begin(), vEdges.end(), edges[k]);
        if(iter != vEdges.end())
          (*iter).pFace1 = pFace;
        else
          vEdges.push_back(edges[k]);
      }

      nProcessedCount++;
      if(nProcessedCount % 100 == 0)
        CObject::set_progress(100 * nProcessedCount / nTotalCount);
    }

    size_t nEdgeCount = vEdges.size();
    for(size_t l = 0; l < nEdgeCount; l++)
    {
      const CEdge& edge = vEdges.at(l);
      if(edge.pFace1 != NULL)
        continue; // edges belonging to the wireframe belong to only one face in this region.

      if(edge.pNode0->is_wireframe_node() && edge.pNode1->is_wireframe_node())
      {
        m_vWireFrame.push_back(CEdgeVertex(edge.pNode0->pos));
        m_vWireFrame.push_back(CEdgeVertex(edge.pNode1->pos));
      }
    }
  }

  m_bWireframeReady = true;
}

void CTrackDraw::set_region_under_cursor(CRegion* pReg)
{
  m_pRegUnderCursor = pReg;
  if(m_pRegUnderCursor == NULL)
    return;

  CFace* pFace = NULL;
  m_vSelRegVert.clear();

  const CFacesCollection& faces = m_pRegUnderCursor->vFaces;
  size_t nFaceCount = faces.size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    pFace = faces.at(i);
    m_vSelRegVert.push_back(CEdgeVertex(pFace->p0->pos));
    m_vSelRegVert.push_back(CEdgeVertex(pFace->p1->pos));

    m_vSelRegVert.push_back(CEdgeVertex(pFace->p1->pos));
    m_vSelRegVert.push_back(CEdgeVertex(pFace->p2->pos));

    m_vSelRegVert.push_back(CEdgeVertex(pFace->p2->pos));
    m_vSelRegVert.push_back(CEdgeVertex(pFace->p0->pos));
  }
}

// Cross-sections of the calculators suppport:
void CTrackDraw::set_cross_sections_array(CExternalFaces* pFacesColl)
{
  m_vCrossSectVert.clear();
  if(pFacesColl == NULL)
    return;

  CFace* pFace = NULL;
  size_t nFaceCount = pFacesColl->size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    pFace = pFacesColl->at(i);
    m_vCrossSectVert.push_back(CEdgeVertex(pFace->p0->pos));
    m_vCrossSectVert.push_back(CEdgeVertex(pFace->p1->pos));

    m_vCrossSectVert.push_back(CEdgeVertex(pFace->p1->pos));
    m_vCrossSectVert.push_back(CEdgeVertex(pFace->p2->pos));

    m_vCrossSectVert.push_back(CEdgeVertex(pFace->p2->pos));
    m_vCrossSectVert.push_back(CEdgeVertex(pFace->p0->pos));
  }
}

void CTrackDraw::draw()
{
  if(m_pTracker == NULL)
    return;

  m_bBusy = true;
  build_arrays();

  set_global();
  set_projection();

  glClearColor(m_fBkgrRed, m_fBkgrGreen, m_fBkgrBlue, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glEnableClientState(GL_VERTEX_ARRAY);

  if(m_bWireframeReady && m_bFacesReady)  // execution could be terminated by the user, so this check is necessary.
  {
    draw_geometry();
    draw_tracks();

    draw_contours();
    draw_region_under_cursor();
    draw_cross_sections();
    draw_selected_faces();
  }

  draw_axes();

  glDisableClientState(GL_VERTEX_ARRAY);

  glFinish();
  glFlush();

  BOOL bOK = SwapBuffers(m_hDC);
  m_bBusy = false;

  set_progress("Ready", -1);
}

void CTrackDraw::draw_tracks()
{
  if(m_bDrawTracks)
    m_ColoredTracks.draw();
}

void CTrackDraw::draw_region_under_cursor()
{
  if(m_pRegUnderCursor == NULL || m_vSelRegVert.size() == 0)
    return;

  glDisable(GL_LIGHTING);
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);

  glColor3ub(120, 190, 230);
  UINT nStride = 3 * sizeof(GLdouble);
  glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_vSelRegVert[0].x));

  glDrawArrays(GL_LINES, 0, m_vSelRegVert.size());

  glEnable(GL_DEPTH_TEST);
}

void CTrackDraw::draw_cross_sections()
{
  size_t nSize = m_vCrossSectVert.size();
  if(nSize == 0)
    return;

  glDisable(GL_LIGHTING);
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);

  glColor3ub(255, 200, 0);
  UINT nStride = 3 * sizeof(GLdouble);
  glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_vCrossSectVert[0].x));

  glDrawArrays(GL_LINES, 0, nSize);

  glEnable(GL_DEPTH_TEST);
}

void CTrackDraw::draw_geometry()
{
  switch(m_nDrawMode)
  {
    case dmNone: return;
    case dmWire: draw_selected_regions(); draw_wire(); draw_cs_flat(); return;
    case dmFlatAndWire: draw_wire(); draw_flat(); draw_selected_regions(); draw_cs_flat(); draw_norm(); return;
    case dmFlatOnly: draw_flat(); draw_selected_regions(); draw_cs_flat(); draw_norm(); return;
  }
}

void CTrackDraw::draw_wire()
{
  size_t nSize = m_vWireFrame.size();
  if(nSize == 0)
    return;

  glDisable(GL_LIGHTING);
  glDisable(GL_ALPHA_TEST);
//  glDisable(GL_DEPTH_TEST);
//  glDisable(GL_BLEND);
//  glLineWidth(1.5);
//  glEnable(GL_LINE_SMOOTH);

  glColor3ub(128, 128, 128);
  UINT nStride = 3 * sizeof(GLdouble);
  glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_vWireFrame[0].x));

  glDrawArrays(GL_LINES, 0, nSize);

//  glEnable(GL_DEPTH_TEST);
//  glLineWidth(1);
//  glDisable(GL_LINE_SMOOTH);
}

void CTrackDraw::draw_norm()
{
  if(!m_bDrawNorm)
    return;

  size_t nSize = m_vAuxLines.size();
  if(nSize == 0)
    return;

  glDisable(GL_LIGHTING);
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);

  glColor3ub(0, 0, 255);
  UINT nStride = 3 * sizeof(GLdouble);
  glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_vAuxLines[0].x));

  glDrawArrays(GL_LINES, 0, nSize);

  glEnable(GL_DEPTH_TEST);
}

void CTrackDraw::set_lights()
{
  glEnable(GL_LIGHTING);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

// The first light in the scene:
  float pLight_0_Dir[4] = { 0.5f, -1.0f, 1.0f, 0.0f };  // this is direction of the parallel type of light.
  float pLight_0_Ambient[4] = { m_fMatAmbient, m_fMatAmbient, m_fMatAmbient, 1.0f };
  float pLight_0_Diffuse[4] = { m_fMatDiffuse, m_fMatDiffuse, m_fMatDiffuse, 1.0f };
  float pLight_0_Spec[4] = { 1.0f, 1.0f, 1.0f, 1.0f };

  glLightfv(GL_LIGHT0, GL_POSITION, pLight_0_Dir);
  glLightfv(GL_LIGHT0, GL_AMBIENT, pLight_0_Ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, pLight_0_Diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, pLight_0_Spec);

  glEnable(GL_LIGHT0);

// The second light is of opposite direction:
  float fDiff = 0.1 * m_fMatDiffuse;
  float pLight_1_Dir[4] = { -0.5f, 1.0f, -1.0f, 0.0f };     // this is direction of the parallel type of light.
  float pLight_1_Ambient[4] = { 0, 0, 0, 1.0f };
  float pLight_1_Diffuse[4] = { fDiff, fDiff, fDiff, 1.0f };
  float pLight_1_Spec[4] = { 0.1f, 0.1f, 0.1f, 1.0f };

  glLightfv(GL_LIGHT1, GL_POSITION, pLight_1_Dir);
  glLightfv(GL_LIGHT1, GL_AMBIENT, pLight_1_Ambient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, pLight_1_Diffuse);
  glLightfv(GL_LIGHT1, GL_SPECULAR, pLight_1_Spec);

  glEnable(GL_LIGHT1);
}

void CTrackDraw::set_materials()
{
  glEnable(GL_COLOR_MATERIAL);
    
  glDisable(GL_CULL_FACE);
  glDisable(GL_NORMALIZE);

  glEnable(GL_ALPHA_TEST);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  float specular[4] = { m_fMatSpec, m_fMatSpec, m_fMatSpec, 1.0f };
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);

  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 1.0f - m_fMatShininess);
}

void CTrackDraw::draw_flat()
{
  set_lights();
  set_materials();

  UINT nStride = 6 * sizeof(GLdouble);
  glEnableClientState(GL_NORMAL_ARRAY);

  size_t nAreasCount = m_vClrFaceVert.size();
  for(size_t i = 0; i < nAreasCount; i++)
  {
    const CColoredRegVertices& vFaceVert = m_vClrFaceVert.at(i);
    if(vFaceVert.size() == 0)
      continue;

    glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&vFaceVert[0].x));
    glNormalPointer(GL_DOUBLE, nStride, (const void*)(&vFaceVert[0].nx));

    COLORREF clr = vFaceVert.get_reg_color();
    glColor4ub(GetRValue(clr), GetGValue(clr), GetBValue(clr), (unsigned char)(255 * m_fOpacity));
    glDrawArrays(GL_TRIANGLES, 0, vFaceVert.size());
  }

  glDisableClientState(GL_NORMAL_ARRAY);
}

void CTrackDraw::draw_cs_flat()
{
  if(m_vFaceCrossSectVert.size() == 0)
    return;

  set_lights();
  set_materials();
  glEnableClientState(GL_NORMAL_ARRAY);

  UINT nStride = 6 * sizeof(GLdouble);
  glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_vFaceCrossSectVert[0].x));
  glNormalPointer(GL_DOUBLE, nStride, (const void*)(&m_vFaceCrossSectVert[0].nx));

  glColor4ub(200, 200, 0, (unsigned char)(128 * m_fOpacity));
  glDrawArrays(GL_TRIANGLES, 0, m_vFaceCrossSectVert.size());

  glDisableClientState(GL_NORMAL_ARRAY);
}

void CTrackDraw::draw_selected_regions()
{
  if(m_vFacesSelRegionVert.size() == 0)
    return;

  set_lights();
  set_materials();
  glEnableClientState(GL_NORMAL_ARRAY);

  UINT nStride = 6 * sizeof(GLdouble);
  glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_vFacesSelRegionVert[0].x));
  glNormalPointer(GL_DOUBLE, nStride, (const void*)(&m_vFacesSelRegionVert[0].nx));

//  glColor4ub(130, 200, 130, (unsigned char)(255 * m_fOpacity)); // pale green.
  glColor4ub(0, 200, 200, (unsigned char)(255 * m_fOpacity)); // bright cyan.
  glDrawArrays(GL_TRIANGLES, 0, m_vFacesSelRegionVert.size());

  glDisableClientState(GL_NORMAL_ARRAY);
}

void CTrackDraw::draw_selected_faces()
{
  if(!m_bDrawSelFaces)
    return;

  set_lights();
  set_materials();
  glEnableClientState(GL_NORMAL_ARRAY);

// Draw the selected faces...
  if(m_vSelFacesVert.size() != 0)
  {
    UINT nStride = 6 * sizeof(GLdouble);
    glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_vSelFacesVert[0].x));
    glNormalPointer(GL_DOUBLE, nStride, (const void*)(&m_vSelFacesVert[0].nx));

    glColor4ub(0, 0, 255, (unsigned char)(192 * m_fOpacity)); // semi-transparent bright blue.
    glDrawArrays(GL_TRIANGLES, 0, m_vSelFacesVert.size());
  }

  glDisableClientState(GL_NORMAL_ARRAY);

// ... and the face under the cursor (if any)
  const CRegionsCollection& regions = m_pTracker->get_regions();
  if(m_FaceUnderCursor.nReg < regions.size())
  {
    CRegion* pReg = regions.at(m_FaceUnderCursor.nReg);
    if(m_FaceUnderCursor.nFace < pReg->vFaces.size())
    {
      CFace* pFace = pReg->vFaces.at(m_FaceUnderCursor.nFace);

      CEdgeVertexColl vFaceVert(6);
      vFaceVert.push_back(CEdgeVertex(pFace->p0->pos));
      vFaceVert.push_back(CEdgeVertex(pFace->p1->pos));

      vFaceVert.push_back(CEdgeVertex(pFace->p1->pos));
      vFaceVert.push_back(CEdgeVertex(pFace->p2->pos));

      vFaceVert.push_back(CEdgeVertex(pFace->p2->pos));
      vFaceVert.push_back(CEdgeVertex(pFace->p0->pos));

      glDisable(GL_LIGHTING);
      glDisable(GL_ALPHA_TEST);
      glDisable(GL_DEPTH_TEST);
      glDisable(GL_BLEND);
      glLineWidth(2);

      glColor3ub(255, 0, 0);
      UINT nStride = 3 * sizeof(GLdouble);
      glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&vFaceVert[0].x));
      glDrawArrays(GL_LINES, 0, vFaceVert.size());

      glLineWidth(1);
      glEnable(GL_DEPTH_TEST);
    }
  }
}

void CTrackDraw::set_global()
{
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glShadeModel(GL_FLAT);
}

// Mouse events support
void CTrackDraw::on_mouse_move(const CPoint& point)
{
  if(!m_pTracker || !m_pTracker->is_ready())
    return;

  if(m_nRegime != nRegimeNone)
  {
    CPoint dp = point - m_StartPoint;
    if(m_nRegime == nRegimeMove)
    { // the shifts are in cm.
      m_vShift.x += m_fPix2cm * dp.x;
      m_vShift.y -= m_fPix2cm * dp.y;
    }
    else
    { // degrees, the coefficient is taken from experience.
      if(m_nContext == nContextRotX)
        m_fRotAngleX = -0.2 * dp.y;
      if(m_nContext == nContextRotY)
        m_fRotAngleY = 0.2 * dp.y;
      if(m_nContext == nContextRotZ)
        m_fRotAngleZ = 0.2 * dp.y;
    }

    m_StartPoint = point;
    draw();
  }
  else if(m_bSelRegFlag)
  {
    CRegFacePair face;
    CRay ray = get_view_dir(point);
    CRegion* pReg = intersect(ray, face);
    set_region_under_cursor(pReg);

    if(m_bCtrlPressed && (pReg != NULL))  // if the Control key is down, select all regions, over which the cursor is, on the fly.
    {
      pReg->bSelected = true;
      invalidate_faces();
      invalidate_aux();
    }

    draw();
  }
  else if(m_bSelTrajectFlag && m_bDrawTracks)
  {
    CTrajectSelector sel(m_pTracker->get_tracks(), point.x, point.y);
    m_nTrajUnderCursorId = sel.find_traject();
    draw();
  }
  else if(m_bSelFacesFlag)
  {
    CRay ray = get_view_dir(point);
    CRegion* pReg = intersect(ray, m_FaceUnderCursor);

    if(m_bCtrlPressed && (pReg != NULL) && (m_FaceUnderCursor.nFace != UINT_MAX))  // if the Control key is down, select all faces, over which the cursor is, on the fly.
    {
      CFaceIndices::iterator iter = std::find(m_vSelFaces.begin(), m_vSelFaces.end(), m_FaceUnderCursor);
      if(iter == m_vSelFaces.end())
      {
        m_vSelFaces.push_back(m_FaceUnderCursor);
        invalidate_faces();
        invalidate_aux();
      }
    }

    draw();
  }
  else
  {
    int nOldOvrAxis = m_nOvrAxis;
    get_over_axis(point);
    if(nOldOvrAxis != m_nOvrAxis)
      draw();
  }

  Vector3D w;
  screen_to_world(point, w);

  std::string cTextX = "X = " + dbl_to_str(10 * w.x);
  std::string cTextY = "Y = " + dbl_to_str(10 * w.y);
  std::string cTextZ = "Z = " + dbl_to_str(10 * w.z);

  CMainFrame* pMainWnd = (CMainFrame*)(CParticleTrackingApp::Get()->m_pMainWnd);
  if(pMainWnd == NULL)
    return;

  CMFCStatusBar* pStatusBar = pMainWnd->GetStatusBar();
  if(pStatusBar)
  {
    pStatusBar->SetPaneText(1, cTextX.c_str());
    pStatusBar->SetPaneText(2, cTextY.c_str());
    pStatusBar->SetPaneText(3, cTextZ.c_str());
  }
}

void CTrackDraw::on_mouse_wheel(short nDelta, const CPoint& point)
{
  if(m_pTracker == NULL)
    return;

  if(nDelta > 0)
  {
    m_fScale = 0.95;
  }
  else if(nDelta < 0)
  {
    m_fScale = 1.0526;
  }
}

void CTrackDraw::on_left_button_down(const CPoint& point)
{
  if(m_pTracker == NULL)
    return;

  if(m_bCtrlPressed)
  {
    m_nRegime = nRegimeNone;
  }
  else
  {
    if(m_nContext == nContextMove)
      m_nRegime = nRegimeMove;
    else
      m_nRegime = nRegimeRotate;

    m_StartPoint = point;
  }
}

void CTrackDraw::on_left_button_up(const CPoint& point)
{
  if(m_pTracker == NULL)
    return;

  m_nRegime = nRegimeNone;

  if(m_bSelRegFlag && (m_pRegUnderCursor != NULL) && (m_StartPoint == point))
  {
    m_pRegUnderCursor->bSelected = !m_pRegUnderCursor->bSelected;
    invalidate_faces();
    invalidate_aux();
    draw();
  }
  else if(m_bSelTrajectFlag && (m_nTrajUnderCursorId != -1) && (m_StartPoint == point))
  {
    std::vector<size_t>::iterator iter = std::find(m_vSelTrackIds.begin(), m_vSelTrackIds.end(), (size_t)m_nTrajUnderCursorId);
    if(iter == m_vSelTrackIds.end())
      m_vSelTrackIds.push_back(m_nTrajUnderCursorId);
    else
      m_vSelTrackIds.erase(iter);

    m_nTrajUnderCursorId = -1;
    draw();
  }
  else if(m_bSelFacesFlag && (m_FaceUnderCursor.nFace != UINT_MAX) && (m_StartPoint == point))
  {
    CFaceIndices::iterator iter = std::find(m_vSelFaces.begin(), m_vSelFaces.end(), m_FaceUnderCursor);
    if(iter == m_vSelFaces.end())
      m_vSelFaces.push_back(m_FaceUnderCursor);
    else
      m_vSelFaces.erase(iter);

    m_FaceUnderCursor.nFace = UINT_MAX;
    m_FaceUnderCursor.nReg = UINT_MAX;
    invalidate_faces();
    invalidate_aux();
    draw();
  }
}

void CTrackDraw::on_left_button_dblclk(const CPoint& point)
{
  if(m_pTracker == NULL || m_bDrawNorm || !m_bRotCenter) // m_bDrawNorm in fact means drawing the selected Dirichlet cell.
    return;

// If the Dirichlet cell drawing regime is NOT set, the double click will set the rotation center:
  CRay ray = get_view_dir(point);
  CBox box = m_pTracker->get_box();

  Vector3D v1, v2;
  for(int nface = 0; nface < 5; nface += 2)
  {
    if(box.intersect_plane(ray, nface, v1) && box.intersect_plane(ray, nface + 1, v2))
    {
      m_vCenter = 0.5 * (v1 + v2);
      return;
    }
  }

  m_vCenter = box.get_center();
}

void CTrackDraw::on_right_button_up(const CPoint& point)
{
  if(m_pTracker == NULL)
    return;

  set_view(m_nOvrAxis);
}

bool CTrackDraw::capture_image()
{
  return m_WndImage.CaptureWindow(m_hWnd);
}

bool CTrackDraw::save_image(const CString& cFileName)
{
  HRESULT hRes = m_WndImage.Save(cFileName);
  if(FAILED(hRes))
    return false;

  return true;
}

// Cursor coordinates in the status line:
void CTrackDraw::screen_to_world(const CPoint& p, Vector3D& w, bool bWorldDepth) const
{
  double pModelMtx[16], pProjMtx[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, pModelMtx);
  glGetDoublev(GL_PROJECTION_MATRIX, pProjMtx);

  int pViewPort[4];
  glGetIntegerv(GL_VIEWPORT, pViewPort);

  float x = float(p.x), y = (float)pViewPort[3] - float(p.y), z;
  glReadPixels(p.x, int(y), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z);
  if(bWorldDepth)
    gluUnProject(x, y, z, pModelMtx, pProjMtx, pViewPort, &w.x, &w.y, &w.z);
  else
    gluUnProject(x, y, 0.5f, pModelMtx, pProjMtx, pViewPort, &w.x, &w.y, &w.z);
}

bool CTrackDraw::world_to_screen(const Vector3D& pos, CPoint& scr) const
{
  double pModelMtx[16], pProjMtx[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, pModelMtx);
  glGetDoublev(GL_PROJECTION_MATRIX, pProjMtx);

  int pViewPort[4];
  glGetIntegerv(GL_VIEWPORT, pViewPort);

  double x, y, z;
  GLint nRes = gluProject(pos.x, pos.y, pos.z, pModelMtx, pProjMtx, pViewPort,&x, &y, &z);
  if(nRes == GLU_FALSE)
    return false;

  scr.x = int(x);
  scr.y = int((double)pViewPort[3] - y);

  return true;
}

CRay CTrackDraw::get_view_dir(const CPoint& point) const
{
  double pModelMtx[16], pProjMtx[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, pModelMtx);
  glGetDoublev(GL_PROJECTION_MATRIX, pProjMtx);

  int pViewPort[4];
  glGetIntegerv(GL_VIEWPORT, pViewPort);

  Vector3D vPos, vDir;
  float x = float(point.x), y = (float)pViewPort[3] - float(point.y);
  gluUnProject(x, y, 0, pModelMtx, pProjMtx, pViewPort, &vPos.x, &vPos.y, &vPos.z);
  gluUnProject(x, y, 1, pModelMtx, pProjMtx, pViewPort, &vDir.x, &vDir.y, &vDir.z);

  vDir -= vPos;
  vDir.normalize();
  return CRay(vPos, vDir);
}

std::string CTrackDraw::dbl_to_str(double val) const
{
  double trunc = val >= 0 ? 0.001 * int(0.5 + 1000 * val) : 0.001 * int(-0.5 + 1000 * val);
  std::ostringstream buffer;
  buffer << trunc;
  return buffer.str();
}

// Text output support:
void CTrackDraw::get_inv_rot(double* pInvMtx)
{
  double pRotMtx[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, pRotMtx);

  pInvMtx[15] = 1;
  pInvMtx[3] = pInvMtx[7] = pInvMtx[11] = pInvMtx[12] = pInvMtx[13] = pInvMtx[14] = 0;
  for(UINT i = 0; i < 3; i++)
    for(UINT j = 0; j < 3; j++)
      pInvMtx[j + 4 * i] = pRotMtx[i + 4 * j];
}

void CTrackDraw::draw_axes()
{
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glLineWidth(2.5f);

// Screen position of the axes triad:
  long nResX, nResY;
  get_resolution(nResX, nResY);

// World position of the axes:
  Vector3D vWorldPos(0, 0, 0);
  CPoint cScrPos(int(0.92 * nResX), int(0.92 * nResY));
  screen_to_world(cScrPos, vWorldPos, false);

  m_fAxisLen = 0.05 * m_fScrWidth;

  m_pAxesTriad[0] = m_pAxesTriad[6] = m_pAxesTriad[9] = m_pAxesTriad[12] = m_pAxesTriad[15] = vWorldPos.x;
  m_pAxesTriad[1] = m_pAxesTriad[4] = m_pAxesTriad[7] = m_pAxesTriad[13] = m_pAxesTriad[16] = vWorldPos.y;
  m_pAxesTriad[2] = m_pAxesTriad[5] = m_pAxesTriad[8] = m_pAxesTriad[11] = m_pAxesTriad[14] = vWorldPos.z;

  m_pAxesTriad[3] = vWorldPos.x + m_fAxisLen;
  m_pAxesTriad[10] = vWorldPos.y + m_fAxisLen;
  m_pAxesTriad[17] = vWorldPos.z + m_fAxisLen;

  UINT nStride = 3 * sizeof(GLdouble);
  glColor3ub(255, 0, 0);
  glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_pAxesTriad[0]));
  glDrawArrays(GL_LINES, 0, 2);

  glColor3ub(0, 128, 0);
  glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_pAxesTriad[6]));
  glDrawArrays(GL_LINES, 0, 2);

  glColor3ub(0, 0, 255);
  glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_pAxesTriad[12]));
  glDrawArrays(GL_LINES, 0, 2);

  double fHalfH = 0.0035 * m_fScrWidth;
  double fHalfW = 0.67 * fHalfH;
  for(UINT i = 0; i < 48; i++)
    m_pXYZ[i] = 0.;

// X:
  m_pXYZ[0] = -fHalfW; m_pXYZ[1] = +fHalfH;
  m_pXYZ[3] = +fHalfW; m_pXYZ[4] = -fHalfH;

  m_pXYZ[6] = +fHalfW; m_pXYZ[7] = +fHalfH;
  m_pXYZ[9] = -fHalfW; m_pXYZ[10] = -fHalfH;

// Y:
  m_pXYZ[12] = -fHalfW; m_pXYZ[13] = +fHalfH;
  m_pXYZ[15] = 0; m_pXYZ[16] = 0;

  m_pXYZ[18] = +fHalfW; m_pXYZ[19] = +fHalfH;
  m_pXYZ[21] = 0; m_pXYZ[22] = 0;

  m_pXYZ[24] = 0; m_pXYZ[25] = 0;
  m_pXYZ[27] = 0; m_pXYZ[28] = -fHalfH;

// Z:
  fHalfH *= 0.9;
  m_pXYZ[30] = -fHalfW; m_pXYZ[31] = +fHalfH;
  m_pXYZ[33] = +fHalfW; m_pXYZ[34] = +fHalfH;

  m_pXYZ[36] = +fHalfW; m_pXYZ[37] = +fHalfH;
  m_pXYZ[39] = -fHalfW; m_pXYZ[40] = -fHalfH;

  m_pXYZ[42] = -fHalfW; m_pXYZ[43] = -fHalfH;
  m_pXYZ[45] = +1.05 * fHalfW; m_pXYZ[46] = -fHalfH;

// The minus sign:
  m_pXYZ[48] = -3 * fHalfW;  m_pXYZ[51] = -1.5 * fHalfW;

  double pInvRot[16];
  get_inv_rot(pInvRot); // the inversion of the rotational part of the model view matrix.

  glMatrixMode(GL_MODELVIEW);

  m_fLetDist = 1.13 * m_fAxisLen;
  draw_letter(pInvRot, vWorldPos, dirPlusX);
  draw_letter(pInvRot, vWorldPos, dirPlusY);
  draw_letter(pInvRot, vWorldPos, dirPlusZ);
  if(m_nOvrAxis == dirMinusX)
    draw_letter(pInvRot, vWorldPos, dirMinusX);
  else if(m_nOvrAxis == dirMinusY)
    draw_letter(pInvRot, vWorldPos, dirMinusY);
  else if(m_nOvrAxis == dirMinusZ)
    draw_letter(pInvRot, vWorldPos, dirMinusZ);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glLineWidth(1.0f);
}

void CTrackDraw::draw_letter(double* pInvRot, const Vector3D& vWorldPos, int nLetterType)
{
  UINT nStart, nLength;
  Vector3D vLetterPos = vWorldPos;
  switch(nLetterType)
  {
    case dirPlusX:  vLetterPos.x += m_fLetDist; nStart =  0; nLength = 4; glColor3ub(255, 0, 0); break;
    case dirMinusX: vLetterPos.x -= m_fAxisLen; nStart =  0; nLength = 4; glColor3ub(255, 0, 0); break;
    case dirPlusY:  vLetterPos.y += m_fLetDist; nStart = 12; nLength = 6; glColor3ub(0, 128, 0); break;
    case dirMinusY: vLetterPos.y -= m_fAxisLen; nStart = 12; nLength = 6; glColor3ub(0, 128, 0); break;
    case dirPlusZ:  vLetterPos.z += m_fLetDist; nStart = 30; nLength = 6; glColor3ub(0, 0, 255); break;
    case dirMinusZ: vLetterPos.z -= m_fAxisLen; nStart = 30; nLength = 6; glColor3ub(0, 0, 255); break;
  }

  float fLineWidth = 2.5f;
  if(m_nOvrAxis == nLetterType)
  {
    fLineWidth = 3.0f;
    if(nLetterType == dirMinusX || nLetterType == dirMinusY || nLetterType == dirMinusZ)
      fLineWidth = 2.5f;
  }

  glLineWidth(fLineWidth);

  UINT nStride = 3 * sizeof(GLdouble);

  glPushMatrix();
  glTranslated(vLetterPos.x, vLetterPos.y, vLetterPos.z);
  glMultMatrixd(pInvRot);   // rotate the text making it look always into the camera.
  glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_pXYZ[nStart]));
  glDrawArrays(GL_LINES, 0, nLength);
  if(nLetterType == dirMinusX || nLetterType == dirMinusY || nLetterType == dirMinusZ)
  {
    glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_pXYZ[48]));
    glDrawArrays(GL_LINES, 0, 2);
  }
  glPopMatrix();
}

void CTrackDraw::get_over_axis(const CPoint& point)
{
  m_nOvrAxis = 0;

  long nResX, nResY, nLimX, nLimY;
  get_resolution(nResX, nResY);

  nLimX = 3 * nResX / 4;
  if(point.x < nLimX)
    return;

  nLimY = 3 * nResY / 4;
  if(point.y < nLimY)
    return;

// World position of the axes:
  Vector3D vWorldPos(0, 0, 0);
  CPoint cScrPos(int(0.92 * nResX), int(0.92 * nResY));

  double pModelMtx[16], pProjMtx[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, pModelMtx);
  glGetDoublev(GL_PROJECTION_MATRIX, pProjMtx);

  int pViewPort[4];
  glGetIntegerv(GL_VIEWPORT, pViewPort);

  float x = float(cScrPos.x), y = (float)pViewPort[3] - float(cScrPos.y);
  gluUnProject(x, y, 0.5f, pModelMtx, pProjMtx, pViewPort, &vWorldPos.x, &vWorldPos.y, &vWorldPos.z);

  double fScrX, fScrY, fScrZ;
  const int nTol = 15;

// Try X:
  Vector3D vAxis = Vector3D(vWorldPos.x + m_fLetDist, vWorldPos.y, vWorldPos.z);
  gluProject(vAxis.x, vAxis.y, vAxis.z, pModelMtx, pProjMtx, pViewPort, &fScrX, &fScrY, &fScrZ);
  double fScrRevY = (double)pViewPort[3] - fScrY;
  if((point.x >= fScrX - nTol) && (point.x <= fScrX + nTol) && (point.y >= fScrRevY - nTol) && (point.y <= fScrRevY + nTol))
  {
    m_nOvrAxis = dirPlusX;
    return;
  }
// Try -X:
  vAxis = Vector3D(vWorldPos.x - m_fAxisLen, vWorldPos.y, vWorldPos.z);
  gluProject(vAxis.x, vAxis.y, vAxis.z, pModelMtx, pProjMtx, pViewPort, &fScrX, &fScrY, &fScrZ);
  fScrRevY = (double)pViewPort[3] - fScrY;
  if((point.x >= fScrX - nTol) && (point.x <= fScrX + nTol) && (point.y >= fScrRevY - nTol) && (point.y <= fScrRevY + nTol))
  {
    m_nOvrAxis = dirMinusX;
    return;
  }
// Try Y:
  vAxis = Vector3D(vWorldPos.x, vWorldPos.y + m_fLetDist, vWorldPos.z);
  gluProject(vAxis.x, vAxis.y, vAxis.z, pModelMtx, pProjMtx, pViewPort, &fScrX, &fScrY, &fScrZ);
  fScrRevY = (double)pViewPort[3] - fScrY;
  if((point.x >= fScrX - nTol) && (point.x <= fScrX + nTol) && (point.y >= fScrRevY - nTol) && (point.y <= fScrRevY + nTol))
  {
    m_nOvrAxis = dirPlusY;
    return;
  }
// Try -Y:
  vAxis = Vector3D(vWorldPos.x, vWorldPos.y - m_fAxisLen, vWorldPos.z);
  gluProject(vAxis.x, vAxis.y, vAxis.z, pModelMtx, pProjMtx, pViewPort, &fScrX, &fScrY, &fScrZ);
  fScrRevY = (double)pViewPort[3] - fScrY;
  if((point.x >= fScrX - nTol) && (point.x <= fScrX + nTol) && (point.y >= fScrRevY - nTol) && (point.y <= fScrRevY + nTol))
  {
    m_nOvrAxis = dirMinusY;
    return;
  }
// Try Z:
  vAxis = Vector3D(vWorldPos.x, vWorldPos.y, vWorldPos.z + m_fLetDist);
  gluProject(vAxis.x, vAxis.y, vAxis.z, pModelMtx, pProjMtx, pViewPort, &fScrX, &fScrY, &fScrZ);
  fScrRevY = (double)pViewPort[3] - fScrY;
  if((point.x >= fScrX - nTol) && (point.x <= fScrX + nTol) && (point.y >= fScrRevY - nTol) && (point.y <= fScrRevY + nTol))
  {
    m_nOvrAxis = dirPlusZ;
    return;
  }
// Try -Z:
  vAxis = Vector3D(vWorldPos.x, vWorldPos.y, vWorldPos.z - m_fAxisLen);
  gluProject(vAxis.x, vAxis.y, vAxis.z, pModelMtx, pProjMtx, pViewPort, &fScrX, &fScrY, &fScrZ);
  fScrRevY = (double)pViewPort[3] - fScrY;
  if((point.x >= fScrX - nTol) && (point.x <= fScrX + nTol) && (point.y >= fScrRevY - nTol) && (point.y <= fScrRevY + nTol))
  {
    m_nOvrAxis = dirMinusZ;
    return;
  }
}

void CTrackDraw::set_progress(const char* cJobName, int nPercent) const
{
  CMainFrame* pMainWnd = (CMainFrame*)(CParticleTrackingApp::Get()->m_pMainWnd);
  if(pMainWnd == NULL)
    return;

  CMFCStatusBar* pStatusBar = pMainWnd->GetStatusBar();
  if(pStatusBar == NULL)
    return;

  char buff[4];
  std::string cPercentage("");
  if((nPercent >= 0) && (nPercent <= 100) && (_itoa_s(nPercent, buff, 4, 10) == 0))
    cPercentage = std::string("  ") + std::string(buff) + std::string(" %");

  std::string cStatusLine = std::string(cJobName) + cPercentage;

  pStatusBar->SetPaneText(0, cStatusLine.c_str());
}

// Streaming:
void CTrackDraw::save(CArchive& ar)
{
  const UINT nVersion = 9;  // 9 - m_fMatShininess; 8 - m_fMatAmbient, m_fMatDiffuse and m_fMatSpec; 7 - m_vSelFaces; 6 - m_vSelTrackIds; 5 - m_nDrawMode instead of m_bDrawFaces; 4 - m_ColoredTracks; 3 - m_bDrawTracks and m_vContours.
  ar << nVersion;

  ar << m_nDrawMode; // ar << m_bDrawFaces;
  ar << m_fOpacity;
  ar << m_Color;

// Projection parameters:
  ar << m_fScrWidth;
  ar << m_vShift.x;
  ar << m_vShift.y;

  double pMtx[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, pMtx);
  for(UINT i = 0; i < 16; i++)
    ar << pMtx[i];

// Since version 1:
  ar << m_BkgrColor;

  size_t nHiddenRegCount = m_vHiddenRegNames.size();
  ar << nHiddenRegCount;
  for(size_t j = 0; j < nHiddenRegCount; j++)
  {
    CString cName(m_vHiddenRegNames.at(j).c_str());
    ar << cName;
  }

  ar << m_bDrawTracks;
  size_t nContCount = m_vContours.size();
  ar << nContCount;
  for(size_t k = 0; k < nContCount; k++)
  {
    CColorContour* pObj = m_vContours.at(k);
    pObj->save(ar);
  }

  m_ColoredTracks.save(ar); // since version 4.

  size_t nSelTracksCount = m_vSelTrackIds.size();
  ar << nSelTracksCount;
  for(size_t m = 0; m < nSelTracksCount; m++)
    ar << m_vSelTrackIds.at(m);

  size_t nSelFacesCount = m_vSelFaces.size();
  ar << nSelFacesCount;
  for(size_t n = 0; n < nSelFacesCount; n++)
  {
    const CRegFacePair& face = m_vSelFaces.at(n);
    ar << face.nReg;
    ar << face.nFace;
  }

  ar << m_fMatAmbient;
  ar << m_fMatDiffuse;
  ar << m_fMatSpec;
  ar << m_fMatShininess;
}

void CTrackDraw::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  if(nVersion >= 5)
  {
    ar >> m_nDrawMode;
  }
  else
  {
    bool bDrawFaces;
    ar >> bDrawFaces;
    m_nDrawMode = bDrawFaces ? dmFlatAndWire : dmNone;
  }

  ar >> m_fOpacity;
  ar >> m_Color;

// Projection parameters:
  ar >> m_fScrWidth;
  ar >> m_vShift.x;
  ar >> m_vShift.y;

  double pMtx[16];
  for(UINT i = 0; i < 16; i++)
    ar >> pMtx[i];

  if(nVersion == 1)
  {
    bool bHideSymm;
    ar >> bHideSymm;
  }
  if(nVersion >= 1)
  {
    ar >> m_BkgrColor;
    set_bkgr_color(m_BkgrColor);
  }

  if(nVersion >= 2)
  {
    size_t nHiddenRegCount;
    ar >> nHiddenRegCount;
    m_vHiddenRegNames.clear();
    for(size_t j = 0; j < nHiddenRegCount; j++)
    {
      CString cName;
      ar >> cName;
      std::string sName((const char*)cName);
      m_vHiddenRegNames.push_back(sName);
    }
  }

  if(nVersion >= 3)
  {
    ar >> m_bDrawTracks;
    size_t nContCount;
    ar >> nContCount;
    m_vContours.clear();
    for(size_t k = 0; k < nContCount; k++)
    {
      CColorContour* pObj = new CColorContour();
      pObj->load(ar);
      m_vContours.push_back(pObj);
    }
  }

  if(nVersion >= 4)
    m_ColoredTracks.load(ar);

  if(nVersion >= 6)
  {
    m_vSelTrackIds.clear();
    size_t nSelTracksCount, nId;
    ar >> nSelTracksCount;
    for(size_t m = 0; m < nSelTracksCount; m++)
    {
      ar >> nId;
      m_vSelTrackIds.push_back(nId);
    }
  }

  if(nVersion >= 7)
  {
    m_vSelFaces.clear();
    size_t nSelFacesCount;
    ar >> nSelFacesCount;
    UINT nReg, nFace;
    for(size_t n = 0; n < nSelFacesCount; n++)
    {
      ar >> nReg;
      ar >> nFace;
      m_vSelFaces.push_back(CRegFacePair(nReg, nFace));
    }
  }

  if(nVersion >= 8)
  {
    ar >> m_fMatAmbient;
    ar >> m_fMatDiffuse;
    ar >> m_fMatSpec;
  }

  if(nVersion >= 9)
    ar >> m_fMatShininess;

  glMatrixMode(GL_MODELVIEW);
  glLoadMatrixd(pMtx);

  m_bNewData = false;
}

void CTrackDraw::set_hidden_reg_names()
{
  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();
  for(size_t j = 0; j < nRegCount; j++)
  {
    CRegion* pReg = regions.at(j);
    if(pReg->bCrossSection)
      continue;   // visibility of cross-sections is controlled ONLY by the "Enable Plane" check box on the "Drawing" tab.

    pReg->bEnabled = std::find(m_vHiddenRegNames.begin(), m_vHiddenRegNames.end(), pReg->sName) == m_vHiddenRegNames.end();
  }

  invalidate_faces();
  invalidate_aux();
}

CString CTrackDraw::get_hidden_names_str() const
{
  return CObject::compile_string(m_vHiddenRegNames);
}

void CTrackDraw::enter_sel_context(CNamesVector* pRegNames, bool bAllowSelect)
{
  if(m_pTracker == NULL)
    return;

  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();
  for(size_t j = 0; j < nRegCount; j++)
  {
    CRegion* pReg = regions.at(j);
    pReg->bSelected = std::find(pRegNames->begin(), pRegNames->end(), pReg->sName) != pRegNames->end();
  }

  m_bSelRegFlag = bAllowSelect;
  invalidate_faces();
  invalidate_aux();
}

void CTrackDraw::exit_sel_context(CNamesVector* pRegNames)
{
  if(m_pTracker == NULL)
    return;

  pRegNames->clear();
  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();
  for(size_t j = 0; j < nRegCount; j++)
  {
    CRegion* pReg = regions.at(j);
    if(pReg->bSelected)
      pRegNames->push_back(pReg->sName);
  }

// We do not know whose collection pRegNames is, but if its owner is "this", the CRegion::bEnabled flag must be updated.
  if(pRegNames == &m_vHiddenRegNames)
    set_hidden_reg_names();
  else
    CParticleTrackingApp::Get()->SelectedRegionChanged(pRegNames);

  m_bSelRegFlag = false;
  set_region_under_cursor(NULL);
  invalidate_faces();
  invalidate_aux();
}

void CTrackDraw::hide_selected()
{
  if(m_pTracker == NULL)
    return;

  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();
  for(size_t j = 0; j < nRegCount; j++)
  {
    CRegion* pReg = regions.at(j);
    pReg->bSelected = false;
  }

  invalidate_faces();
  invalidate_aux();
}

void CTrackDraw::show_all_regions()
{
  if(m_pTracker == NULL)
    return;

  m_vHiddenRegNames.clear();
  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();
  for(size_t j = 0; j < nRegCount; j++)
  {
    CRegion* pReg = regions.at(j);
    if(pReg->bCrossSection)
      continue; // visibility of cross-sections is controlled ONLY by the "Enable Plane" check box on the "Drawing" tab.

    pReg->bEnabled = true;
  }

// Set visibility flags for all areas in the scene:
  CParticleTrackingApp::Get()->GetSelAreas()->make_all_visible();

  m_vFacesSelRegionVert.clear();
  invalidate_faces();
  invalidate_aux();
}

void CTrackDraw::set_visibility_status(CNamesVector* pRegNames, bool bVisible)
{
  if(m_pTracker == NULL || pRegNames == NULL || pRegNames->size() == 0)
    return;

  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();
  for(size_t j = 0; j < nRegCount; j++)
  {
    CRegion* pReg = regions.at(j);
    if(pReg->bCrossSection)
      continue; // visibility of cross-sections is controlled ONLY by the "Enable Plane" check box on the "Drawing" tab.

    if(std::find(pRegNames->begin(), pRegNames->end(), pReg->sName) != pRegNames->end())
      pReg->bEnabled = bVisible;
  }

  invalidate_faces();
  invalidate_aux();
}

CRegion* CTrackDraw::intersect(const CRay& ray, CRegFacePair& face) const
{
  if(m_pTracker == NULL)
    return NULL;

  UINT nFaceID;
  CRegion* pSelReg = NULL;
  double fMinDist = FLT_MAX, fDist;
  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();
  for(size_t j = 0; j < nRegCount; j++)
  {
    CRegion* pReg = regions.at(j);
    if(pReg->bEnabled && pReg->intersect(ray, fDist, nFaceID) && (fDist < fMinDist))
    {
      face.nReg = j;
      face.nFace = nFaceID;
      fMinDist = fDist;
      pSelReg = pReg;
    }
  }

  return pSelReg;
}

const char* CTrackDraw::get_draw_mode_name(int nMode) const
{
  switch(nMode)
  {
    case dmWire: return _T("Wireframe");
    case dmFlatAndWire: return _T("Flat and Wireframe");
    case dmFlatOnly: return _T("Flat Only");
  }

  return _T("None");
}

CString CTrackDraw::get_sel_faces_square_str() const
{
  double fS = get_sel_faces_square();
  CString sVal = dbl_to_str(fS).c_str();
  return sVal;
}

double CTrackDraw::get_sel_faces_square(bool bSquaredMillimeters) const
{
  if(m_pTracker == NULL)
    return 0;

  double fS = 0;
  size_t nSelCount = m_vSelFaces.size();
  const CRegionsCollection& regions = m_pTracker->get_regions();
  size_t nRegCount = regions.size();
  CRegion* pReg = NULL;

  for(size_t i = 0; i < nSelCount; i++)
  {
    CRegFacePair face = m_vSelFaces.at(i);
    if(face.nReg >= nRegCount)
      continue;

    pReg = regions.at(face.nReg);
    if(face.nFace >= pReg->vFaces.size())
      continue;

    fS += pReg->vFaces.at(face.nFace)->square();
  }

  return bSquaredMillimeters ? 100 * fS : fS;
}

//-------------------------------------------------------------------------------------------------
//  Contours and vector plots:
//-------------------------------------------------------------------------------------------------
void CTrackDraw::add_contour()
{
  m_vContours.push_back(new CColorContour());
}

void CTrackDraw::remove_contour(CColorContour* pItem)
{
  size_t nCount = m_vContours.size();
  for(size_t i = 0; i < nCount; i++)
  {
    if(pItem == m_vContours.at(i))
    {
      delete pItem;
      m_vContours.erase(m_vContours.begin() + i);
      break;
    }
  }
}

void CTrackDraw::invalidate_contour(DWORD_PTR pRegNames)
{
  size_t nCount = m_vContours.size();
  for(size_t i = 0; i < nCount; i++)
  {
    CColorContour* pObj = m_vContours.at(i);
    if(pObj->get_drawn_reg_names_ptr() == pRegNames)
    {
      pObj->invalidate();
      break;
    }
  }
}

void CTrackDraw::draw_contours()
{
  size_t nContCount = get_contours_count();
  for(size_t i = 0; i < nContCount; i++)
  {
    CColorContour* pObj = get_contour(i);
    if(pObj->get_enable_image())
      pObj->draw();
  }
}

void CTrackDraw::set_phi_to_nodes() const
{
  CFieldDataColl* pAllFields = CParticleTrackingApp::Get()->GetFields();
  size_t nFieldCount = pAllFields->size();

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CFieldPtbCollection& vPtbs = pObj->get_field_ptb();
  size_t nPtbCount = vPtbs.size();

  if((nFieldCount == 0) && (nPtbCount == 0))
    return;

  CElectricFieldData* pField = NULL;
  CFieldPerturbation* pPtb = NULL;

  CNodesVector& vNodes = pObj->get_nodes();
  size_t nNodeCount = vNodes.size();

  for(size_t i = 0; i < nNodeCount; i++)
  {
    CNode3D& node = vNodes.at(i);
    node.phi = 0;
    for(size_t j = 0; j < nFieldCount; j++)
    {
      pField = pAllFields->at(j);
      if(pField->get_enable_vis())
      {  
        if(pField->get_type() != CElectricFieldData::typeMirror)
          node.phi += pField->get_phi(i) * pField->get_scale(); // Note that if the j-th field is not ready, its get_phi(i) returns 0.
        else
          node.phi += pField->get_clmb_phi(i);
      }
    }

    for(size_t k = 0; k < nPtbCount; k++)
    {
      pPtb = vPtbs.at(k);
      if(pPtb->get_enable())
        node.phi += pPtb->get_phi(i);
    }
  }
}

bool CTrackDraw::save_sel_faces(const char* pFile)
{
  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, pFile, (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return false;

  size_t nPairsCount = m_vSelFaces.size();
  fprintf(pStream, "%zd\n", nPairsCount);
  for(size_t i = 0; i < nPairsCount; i++)
  {
    const CRegFacePair& face = m_vSelFaces.at(i);
    fprintf(pStream, "%d, %d\n", face.nReg, face.nFace);
  }

  fclose(pStream);
  return true;
}

bool CTrackDraw::load_sel_faces(const char* pFile)
{
  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, pFile, (const char*)("r"));
  if(nErr != 0 || pStream == 0)
    return false;

  UINT nReg, nFace;
  size_t nPairsCount = 0;
  int nRes = fscanf_s(pStream, "%zd\n", &nPairsCount);
  if(nRes == EOF || nRes == 0 || nPairsCount == 0)
  {
    fclose(pStream);
    return false;
  }

  m_vSelFaces.clear();
  m_vSelFaces.resize(nPairsCount, CRegFacePair());
  for(size_t i = 0; i < nPairsCount; i++)
  {
    nRes = fscanf_s(pStream, "%d, %d", &nReg, &nFace);
    if(nRes == EOF || nRes == 0)
    {
      fclose(pStream);
      return false;
    }

    m_vSelFaces.at(i).nReg = nReg;
    m_vSelFaces.at(i).nFace = nFace;
  }

  fclose(pStream);
  invalidate_aux();
  return true;
}

};  // namespace EvaporatingParticle