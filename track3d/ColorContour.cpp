
#include "stdafx.h"
#include "ColorContour.h"
#include "ParticleTracking.h"
#include "TrackItem.h"
#include <float.h>

#include <GL/gl.h>
#include <gl/glu.h>

namespace EvaporatingParticle
{

CColorContour::CColorContour()
{
  set_default();
}

CColorContour::~CColorContour()
{
}

void CColorContour::set_default()
{
  CColorImage::set_default();

  m_nVar = varPress;
  m_bDrawLines = false;

  for(UINT i = 0; i < varCount; i++)
  {
    m_vUserDefMin[i] = 0;
    m_vUserDefMax[i] = 0;
  }
}

bool CColorContour::build_draw_array()
{
  m_vIsolines.clear();
  m_vFaceColors.clear();
  m_vFaceVert.clear();
  get_phi();  // the electric potential visualization requires summation over all enabled fields.

  if(m_nLevelsCount < 2)
    return false;

  get_glob_min_max();  // minimum and maximum of the contoured parameter values over all simulation domain.
  
  get_values_array();
  if(m_vValues.size() == 0)
    return false;

  get_colors_array();

  RGB_Color clr;
  double fCurrLevel, fKsi02, fKsi01, fKsi12;
  double pFaceValues[3];  // face values of the parameter to be contoured.
// Face vertices arranged in such a way that pFaceVert[0] corresponds to the minimal value, pFaceVert[2] corresponds to the maximal
// value and pFaceVert[1] corresponds to the intermediate value of the parameter.
  Vector3D pFaceVert[3];
  Vector3D* p02 = new Vector3D[m_nLevelsCount];
  Vector3D* p01 = new Vector3D[m_nLevelsCount];
  Vector3D* p12 = new Vector3D[m_nLevelsCount];

  CRegionsCollection vRegions = get_regions();
  size_t nRegCount = vRegions.size();                     /*                                                                      */
  for(size_t i = 0; i < nRegCount; i++)                   /*                P1, Umin <= U1 <= Umax                                */
  {                                                       /*                   /\                                                 */
    CRegion* pReg = vRegions.at(i);                       /*                  /  \  / Current Level                               */
    size_t nFaceCount = pReg->vFaces.size();              /*                 /    \/                                              */
    for(size_t j = 0; j < nFaceCount; j++)                /*                /     /\ P12, U12 = U1 + fKsi12 * (U2 - U1)           */
    {                                                     /*               /     /  \                                             */
      CFace* pFace = pReg->vFaces.at(j);                  /*              /     /    \                                            */
      reorder_vertices(pFace, pFaceVert, pFaceValues);    /*             /     /      \                                           */
                                                          /*            /_____/________\ P2, Umax                                 */
      if(pFaceValues[2] > pFaceValues[0])                 /*   P0, Umin      / P02                                                */
      {                                                   /*                   U02 = U0 + fKsi02 * (U2 - U0)                      */
        UINT n01 = 0, n12 = 0, n02 = 0, k0 = 0, k1;       /*                                                                      */
        for(UINT k = 0; k < m_nLevelsCount; k++)
        {
          fCurrLevel = m_vValues.at(k);
          if(fCurrLevel < pFaceValues[0] || fCurrLevel > pFaceValues[2])
            continue;

// If there are any isolines intersecting this face, they all must intersect the edge 02.
          fKsi02 = (fCurrLevel - pFaceValues[0]) / (pFaceValues[2] - pFaceValues[0]); // there must be 0 <= fKsi02 <= 1
          p02[k] = pFaceVert[0] + fKsi02 * (pFaceVert[2] - pFaceVert[0]);

          fKsi01 = (fCurrLevel - pFaceValues[0]) / (pFaceValues[1] - pFaceValues[0]);
          fKsi12 = (fCurrLevel - pFaceValues[1]) / (pFaceValues[2] - pFaceValues[1]);

          if((fKsi01 > 0) && (fKsi01 <= 1))
          {
            p01[k] = pFaceVert[0] + fKsi01 * (pFaceVert[1] - pFaceVert[0]);
            n01++;
            if(n02 == 0)
              k0 = k;   // remember the index of first intersection.
          }
          else if((fKsi12 > 0) && (fKsi12 < 1))
          {
            p12[k] = pFaceVert[1] + fKsi12 * (pFaceVert[2] - pFaceVert[1]);
            n12++;
            if(n02 == 0)
              k0 = k;   // remember the index of first intersection.
          }

          n02 = n01 + n12;  // total count of intersections.
        }

        if(n02 == 0) // the whole face {P0, P1, P2} must be colored by a single color.
        {
          add_face(pFaceVert[0], pFaceVert[1], pFaceVert[2], get_color(pFaceValues[0]));
        }
        else if(n01 == 0)
        {
          clr = m_vColors.at(k0);
          add_face(pFaceVert[0], pFaceVert[1], p12[k0], clr);
          add_face(pFaceVert[0], p12[k0], p02[k0], clr);
          add_edge(p02[k0], p12[k0]);

          k1 = k0 + n12;
          for(UINT k = k0 + 1; k < k1; k++)
          {
            clr = m_vColors.at(k);
            add_face(p02[k - 1], p12[k - 1], p12[k], clr);
            add_face(p02[k - 1], p12[k], p02[k], clr);
            add_edge(p02[k], p12[k]);
          }

          clr = m_vColors.at(k1);
          add_face(p02[k1 - 1], p12[k1 - 1], pFaceVert[2], clr);
        }
        else if(n12 == 0)
        {
          clr = m_vColors.at(k0);
          add_face(pFaceVert[0], p01[k0], p02[k0], clr);
          add_edge(p02[k0], p01[k0]);

          k1 = k0 + n01;
          for(UINT k = k0 + 1; k < k1; k++)
          {
            clr = m_vColors.at(k);
            add_face(p02[k - 1], p01[k - 1], p01[k], clr);
            add_face(p02[k - 1], p01[k], p02[k], clr);
            add_edge(p02[k], p01[k]);
          }

          clr = m_vColors.at(k1);
          add_face(p02[k1 - 1], p01[k1 - 1], pFaceVert[1], clr);
          add_face(p02[k1 - 1], pFaceVert[1], pFaceVert[2], clr);
        }
        else  // both n01 and n12 are non-zero.
        {
          clr = m_vColors.at(k0);
          add_face(pFaceVert[0], p01[k0], p02[k0], clr);  // there is at least one intersection with the edge 01.
          add_edge(p02[k0], p01[k0]);

          k1 = k0 + n01;
          for(UINT k = k0 + 1; k < k1; k++)
          {
            clr = m_vColors.at(k);
            add_face(p02[k - 1], p01[k - 1], p01[k], clr);
            add_face(p02[k - 1], p01[k], p02[k], clr);
            add_edge(p02[k], p01[k]);
          }

          clr = m_vColors.at(k1);
          add_face(p02[k1 - 1], p01[k1 - 1], pFaceVert[1], clr);

          k1 = k0 + n02;
          p12[k0 + n01 - 1] = pFaceVert[1];
          for(UINT k = k0 + n01; k < k1; k++)
          {
            clr = m_vColors.at(k);
            add_face(p02[k - 1], p12[k - 1], p12[k], clr);
            add_face(p02[k - 1], p12[k], p02[k], clr);
            add_edge(p02[k], p12[k]);
          }

          clr = m_vColors.at(k1);
          add_face(p02[k1 - 1], p12[k1 - 1], pFaceVert[2], clr);
        }
      }
      else  // the whole face {P0, P1, P2} must be colored by a single color.
      {
        add_face(pFaceVert[0], pFaceVert[1], pFaceVert[2], get_color(pFaceValues[0]));
      }
    }
  }

  delete[] p02;
  delete[] p01;
  delete[] p12;

  m_bReady = true;
  return true;
}

void CColorContour::add_face(const Vector3D& v0, const Vector3D& v1, const Vector3D& v2, const RGB_Color& clr)
{
  m_vFaceVert.push_back(v0);
  m_vFaceColors.push_back(clr);

  m_vFaceVert.push_back(v1);
  m_vFaceColors.push_back(clr);

  m_vFaceVert.push_back(v2);
  m_vFaceColors.push_back(clr);
}

void CColorContour::add_edge(const Vector3D& v0, const Vector3D& v1)
{
  m_vIsolines.push_back(v0);
  m_vIsolines.push_back(v1);
}

void CColorContour::draw()
{
  if(!m_bReady && !build_draw_array())
    return;

// Color contours:
  if(m_vFaceVert.size() != 0)
  {
    glDisable(GL_LIGHTING);

    glEnableClientState(GL_COLOR_ARRAY);

    glVertexPointer(3, GL_DOUBLE, 0, (const void*)(&m_vFaceVert[0].x));
    glColorPointer(3, GL_UNSIGNED_BYTE, 0, (const void*)(&m_vFaceColors[0].red));

    glDrawArrays(GL_TRIANGLES, 0, m_vFaceVert.size());

    glDisableClientState(GL_COLOR_ARRAY);
  }

// Isolines:
  if(m_bDrawLines && (m_vIsolines.size() != 0))
  {
    glDisable(GL_LIGHTING);
    glDisable(GL_ALPHA_TEST);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);

    glColor3ub(255, 255, 255);
    UINT nStride = 3 * sizeof(GLdouble);
    glVertexPointer(3, GL_DOUBLE, nStride, (const void*)(&m_vIsolines[0].x));

    glDrawArrays(GL_LINES, 0, m_vIsolines.size());

    glEnable(GL_DEPTH_TEST);
  }
}

void CColorContour::get_min_max()
{
  if(m_bUserDefRange)
    return;

  const CRegionsCollection& vRegions = get_regions();

  m_fMinVal = FLT_MAX;
  m_fMaxVal = -FLT_MAX;

  double pFaceVal[3];
  size_t nRegCount = vRegions.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    CRegion* pReg = vRegions.at(i);
    size_t nFaceCount = pReg->vFaces.size();
    for(size_t j = 0; j < nFaceCount; j++)
    {
      CFace* pFace = pReg->vFaces.at(j);
      get_face_values_array(pFace, pFaceVal); // pFaceVal is filled by values of the contoured parameter taken from the face vertices.

      for(UINT k = 0; k < 3; k++)
      {
        if(m_fMinVal > pFaceVal[k])
          m_fMinVal = pFaceVal[k];
        if(m_fMaxVal < pFaceVal[k])
          m_fMaxVal = pFaceVal[k];
      }
    }
  }
}

void CColorContour::get_glob_min_max()
{
  if(m_bUserDefRange)
    return;

  m_fMinVal = FLT_MAX;
  m_fMaxVal = -FLT_MAX;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CNodesCollection& vNodes = pObj->get_nodes();
  CNode3D* pNode = NULL;
  double fVal = 0;

  size_t nNodeCount = vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    pNode = vNodes.at(i);
    fVal = get_node_value(pNode);
    if(m_fMinVal > fVal)
      m_fMinVal = fVal;
    if(m_fMaxVal < fVal)
      m_fMaxVal = fVal;
  }
}

void CColorContour::reorder_vertices(CFace* pFace, Vector3D* pFaceVert, double* pVal) const
{
  double pValTmp[3];
  get_face_values_array(pFace, pValTmp);

  UINT jMin = 0, jMax = 0;
  double fMin = pValTmp[0], fMax = pValTmp[0];
  for(UINT j = 1; j < 3; j++)
  {
    if(pValTmp[j] < fMin)
    {
      fMin = pValTmp[j];
      jMin = j;
    }
    if(pValTmp[j] > fMax)
    {
      fMax = pValTmp[j];
      jMax = j;
    }
  }

  UINT jInt;
  CNode3D* pFaceNodes[3] = { pFace->p0, pFace->p1, pFace->p2 };
  if(jMax != jMin)
  {
    UINT jSum = jMin + jMax;
    jInt = jSum == 1 ? 2 : (jSum == 2 ? 1 : 0);

    pFaceVert[0] = pFaceNodes[jMin]->pos;
    pVal[0] = pValTmp[jMin];

    pFaceVert[1] = pFaceNodes[jInt]->pos;
    pVal[1] = pValTmp[jInt];

    pFaceVert[2] = pFaceNodes[jMax]->pos;
    pVal[2] = pValTmp[jMax];
  }
  else
  {
    for(UINT i = 0; i < 3; i++)
    {
      pFaceVert[i] = pFaceNodes[i]->pos;
      pVal[i] = pValTmp[i];
    }
  }
}

void CColorContour::get_face_values_array(CFace* pFace, double* pVal) const
{
  CNode3D* pNodes[3] = { pFace->p0, pFace->p1, pFace->p2 };
  for(UINT i = 0; i < 3; i++)
    pVal[i] = get_node_value(pNodes[i]);
}

double CColorContour::get_node_value(CNode3D* pNode) const
{
  switch(m_nVar)
  {
    case varPress:   return pNode->press;
    case varDens:    return pNode->dens;
    case varTemp:    return pNode->temp;
    case varAbsVel:  return pNode->vel.length();
    case varVelX:    return pNode->vel.x;
    case varVelY:    return pNode->vel.y;
    case varVelZ:    return pNode->vel.z;
    case varAbsClmb: return pNode->clmb.length();
    case varClmbX:   return pNode->clmb.x;
    case varClmbY:   return pNode->clmb.y;
    case varClmbZ:   return pNode->clmb.z;
    case varPhi:     return pNode->phi;
  }

  return 0;
}

bool CColorContour::get_phi() const
{
  CFieldDataColl* pAllFields = CParticleTrackingApp::Get()->GetFields();
  size_t nFieldCount = pAllFields->size();
  if(nFieldCount == 0)
    return false;

  CNode3D* pNode = NULL;
  CElectricFieldData* pField = NULL;
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CNodesCollection& vNodes = pObj->get_nodes();
  size_t nNodeCount = vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    pNode = vNodes.at(i);
    pNode->phi = 0;
    for(size_t j = 0; j < nFieldCount; j++)
    {
      pField = pAllFields->at(j);
      if(pField->get_enable_vis())
      {  
        if(pField->get_type() != CElectricFieldData::typeMirror)
          pNode->phi += pField->get_phi(i) * pField->get_scale(); // Note that if the j-th field is not ready, its get_phi(i) returns 0.
        else
          pNode->phi += pField->get_clmb_phi(i);
      }
    }
  }

  return true;
}

void CColorContour::save(CArchive& ar)
{
  UINT nVersion = 1;
  ar << nVersion;

  CColorImage::save(ar);

  ar << m_nVar;
  ar << m_bDrawLines;

// Since version 1:
  ar << varCount;
  for(UINT i = 0; i < varCount; i++)
  {
    ar << m_vUserDefMin[i];
    ar << m_vUserDefMax[i];
  }
}

void CColorContour::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CColorImage::load(ar);

  ar >> m_nVar;
  ar >> m_bDrawLines;

  if(nVersion >= 1)
  {
    int nCount;
    ar >> nCount;
    for(UINT i = 0; i < nCount; i++)
    {
      if(i >= varCount)
        continue;

      ar >> m_vUserDefMin[i];
      ar >> m_vUserDefMax[i];
    }

    if(m_bUserDefRange && (m_nVar < varCount))
      restore_user_range();
  }
}

double CColorContour::get_min_val() const
{
  return m_fMinVal * get_multiplier(false);
}

void CColorContour::set_min_val(double fVal)
{
  double fValCGS = fVal * get_multiplier(true);
  if(fabs(m_fMinVal - fValCGS) > Const_Almost_Zero)
  {
    m_fMinVal = fValCGS;
    m_bReady = false;
    if(m_bUserDefRange)
      m_vUserDefMin[m_nVar] = m_fMinVal;
  }
}

double CColorContour::get_max_val() const
{
  return m_fMaxVal * get_multiplier(false);
}

void CColorContour::set_max_val(double fVal)
{
  double fValCGS = fVal * get_multiplier(true);
  if(fabs(m_fMaxVal - fValCGS) > Const_Almost_Zero)
  {
    m_fMaxVal = fValCGS;
    m_bReady = false;
    if(m_bUserDefRange)
      m_vUserDefMax[m_nVar] = m_fMaxVal;
  }
}

void CColorContour::restore_user_range()
{
  m_fMinVal = m_vUserDefMin[m_nVar];
  m_fMaxVal = m_vUserDefMax[m_nVar];
}

double CColorContour::get_multiplier(bool bSI_to_CGS) const
{
  switch(m_nVar)
  {
    case varPress: 
    {
      return bSI_to_CGS ? SI_to_CGS_Press : CGS_to_SI_Press;
    }
    case varDens:
    {
      return bSI_to_CGS ? SI_to_CGS_Dens : CGS_to_SI_Dens;
    }
    case varTemp:
    {
      return 1.0;
    }
    case varAbsVel:
    case varVelX:
    case varVelY:
    case varVelZ:
    {
      return bSI_to_CGS ? SI_to_CGS_Vel : CGS_to_SI_Vel;
    }
    case varAbsClmb:
    case varClmbX:
    case varClmbY:
    case varClmbZ:
    {
      return bSI_to_CGS ? SI_to_CGS_Voltage : CGS_to_SI_Voltage;
    }
  }

  return 1.0;
}

const char* CColorContour::get_var_name(int nVar) const
{
  switch(nVar)
  {
    case varPress:    return _T("Pressure");
    case varDens:     return _T("Density");
    case varTemp:     return _T("Temperature");
    case varAbsVel:   return _T("Absolute Velocity");
    case varVelX:     return _T("Velocity X");
    case varVelY:     return _T("Velocity Y");
    case varVelZ:     return _T("Velocity Z");
    case varAbsClmb:  return _T("Absolute Coulomb Field");
    case varClmbX:    return _T("Coulomb Field X");
    case varClmbY:    return _T("Coulomb Field Y");
    case varClmbZ:    return _T("Coulomb Field Z");
    case varPhi:      return _T("Electric Potential");
  }

  return _T("");
}


//---------------------------------------------------------------------------------------
// CColoredTracks.
//---------------------------------------------------------------------------------------
static const int snClrCount = 30;
static const COLORREF sClrArray[snClrCount] =
{
  RGB(255, 0, 0), RGB(255, 160, 0), RGB(255, 255, 0), RGB(180, 255, 0), RGB(0, 255, 0), RGB(0, 255, 255), RGB(0, 160, 255),
  RGB(0, 0, 255), RGB(160, 0, 255), RGB(255, 0, 255), RGB(255, 0, 160), RGB(0, 0, 10), RGB(190, 0, 0), RGB(190, 100, 0), RGB(190, 190, 0),
  RGB(80, 190, 0), RGB(0, 190, 190), RGB(0, 80, 190), RGB(0, 0, 190), RGB(190, 0, 190), RGB(255, 128, 128), RGB(255, 200, 128), RGB(128, 128, 128),
  RGB(255, 255, 190), RGB(190, 255, 190), RGB(190, 255, 255), RGB(190, 190, 255), RGB(100, 230, 190), RGB(255, 160, 200), RGB(255, 255, 255)
};

CColoredTracks::CColoredTracks()
{
  set_default();
}

CColoredTracks::~CColoredTracks()
{
  clear();
}

void CColoredTracks::clear()
{
  size_t nTracksCount = m_vTracksVert.size();
  for(size_t i = 0; i < nTracksCount; i++)
  {
    CVert3DColl& track_vert = m_vTracksVert.at(i);
    track_vert.clear();

    CColorVector& track_colors = m_vTracksColors.at(i);
    track_colors.clear();
  }

  m_vTracksVert.clear();
  m_vTracksColors.clear();
}

void CColoredTracks::set_default()
{
  CColorImage::set_default();

  m_nVarColorMap = tcmDefault;
}

void CColoredTracks::draw()
{
  if(!m_bReady && !build_draw_array())
    return;

  size_t nTracksCount = m_vTracksVert.size();
  if(nTracksCount != 0)
  {
// Draw ordinary tracks:
    glDisable(GL_LIGHTING);
    glEnableClientState(GL_COLOR_ARRAY);

    UINT nStrideDbl = 3 * sizeof(GLdouble);
    UINT nStrideUINT = 3 * sizeof(GLubyte);
    for(size_t i = 0; i < nTracksCount; i++)
    {
      CVert3DColl& vVertArray = m_vTracksVert.at(i);
      CColorVector& vColorArray = m_vTracksColors.at(i);
      if(vVertArray.size() != 0)
      {
        glVertexPointer(3, GL_DOUBLE, nStrideDbl, (const void*)(&vVertArray[0].x));
        glColorPointer(3, GL_UNSIGNED_BYTE, nStrideUINT, (const void*)(&vColorArray[0].red));

        glDrawArrays(GL_LINE_STRIP, 0, vVertArray.size());
      }
    }

    glDisableClientState(GL_COLOR_ARRAY);

    glLineWidth(2);
    glDisable(GL_DEPTH_TEST);
// Draw selected tracks ...
    glColor3ub(255, 0, 0);
    CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
    const CIdsVector& vIds = pDrawObj->get_sel_traject_ids();
    size_t nSelCount = vIds.size();
    for(size_t j = 0; j < nSelCount; j++)
    {
      CVert3DColl& vVertArray = m_vTracksVert.at(vIds.at(j));
      if(vVertArray.size() != 0)
      {
        glVertexPointer(3, GL_DOUBLE, nStrideDbl, (const void*)(&vVertArray[0].x));
        glDrawArrays(GL_LINE_STRIP, 0, vVertArray.size());
      }
    }
// ... and the track under the cursor (if any):
    glColor3ub(255, 0, 255);
    int nTrackUnderCursor = pDrawObj->get_traject_under_cursor_id();
    if((nTrackUnderCursor >= 0) && (nTrackUnderCursor < nTracksCount))
    {
      CVert3DColl& vVertArray = m_vTracksVert.at(nTrackUnderCursor);
      glVertexPointer(3, GL_DOUBLE, nStrideDbl, (const void*)(&vVertArray[0].x));
      glDrawArrays(GL_LINE_STRIP, 0, vVertArray.size());
    }

    glEnable(GL_DEPTH_TEST);
    glLineWidth(1);
  }
}

void CColoredTracks::restore_user_range()
{
  m_fMinVal = m_vUserDefMin[m_nVarColorMap];
  m_fMaxVal = m_vUserDefMax[m_nVarColorMap];
}

bool CColoredTracks::build_draw_array()
{
  clear();
  get_min_max();
  get_values_array();
  if((m_vValues.size() == 0) && (m_nVarColorMap != tcmDefault))
    return false;

  get_colors_array();

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CTrackVector& vTracks = pObj->get_tracks();

  double fVal;
  CTrackItem item;
  size_t nTracksCount = vTracks.size();
  for(size_t i = 0; i < nTracksCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    size_t nItemsCount = track.size();

    CVert3DColl track_vert;
    track_vert.reserve(nItemsCount);

    CColorVector track_colors;
    track_colors.reserve(nItemsCount);

    for(size_t j = 0; j < nItemsCount; j++)
    {
      track.get_track_item(j, item);
      fVal = get_vert_value(item, i);

      track_vert.push_back(item.pos);
      if(m_nVarColorMap == tcmDefault)
        track_colors.push_back(RGB_Color(sClrArray[i % snClrCount]));
      else
        track_colors.push_back(get_color(fVal));
    }

    m_vTracksVert.push_back(track_vert);
    m_vTracksColors.push_back(track_colors);
  }

  m_bReady = true;
  return true;
}

double CColoredTracks::get_vert_value(const CTrackItem& item, size_t nTrackIndex) const
{
  switch(m_nVarColorMap)
  {
    case tcmIonTemp:     return item.temp;
    case tcmIonEqTemp:   return item.tempinf;
    case tcmIonConc:     return item.unfragm;
    case tcmStartRadius: return get_start_radius(nTrackIndex);
  }

  return 0;
}

double CColoredTracks::get_start_radius(size_t nTrackIndex)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CTrackVector& vTracks = pObj->get_tracks();
  if(nTrackIndex >= vTracks.size())
    return 0;

  const CTrack& track = vTracks.at(nTrackIndex);
  if(track.size() < 1)
    return 0;

  CBaseTrackItem* pItem = track.at(0);
  return (pItem->pos - pObj->get_src()->get_src_pos()).length();
}

void CColoredTracks::get_min_max()
{
  if(m_bUserDefRange)
    return;

  double fVal;
  m_fMinVal = FLT_MAX;
  m_fMaxVal = -FLT_MAX;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CTrackVector& vTracks = pObj->get_tracks();

  CTrackItem item;
  size_t nTracksCount = vTracks.size();
  for(size_t i = 0; i < nTracksCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    size_t nItemsCount = track.size();
    for(size_t j = 0; j < nItemsCount; j++)
    {
      track.get_track_item(j, item);
      fVal = get_vert_value(item, i);
      if(fVal < m_fMinVal)
        m_fMinVal = fVal;
      if(fVal > m_fMaxVal)
        m_fMaxVal = fVal;
    }
  }
}

const char* CColoredTracks::get_var_name(int nVar) const
{
  switch(nVar)
  {
    case tcmDefault:     return _T("None");
    case tcmIonTemp:     return _T("Ion Temperature");
    case tcmIonEqTemp:   return _T("Equilibrium Ion Temperature");
    case tcmIonConc:     return _T("Ion Concentration");
    case tcmStartRadius: return _T("Ion Start Radius");
  }

  return _T("None");
}

void CColoredTracks::save(CArchive& ar)
{
  UINT nVersion = 1;
  ar << nVersion;

  CColorImage::save(ar);

  ar << m_nVarColorMap;

// Since version 1:
  ar << tcmCount;
  for(UINT i = 0; i < tcmCount; i++)
  {
    ar << m_vUserDefMin[i];
    ar << m_vUserDefMax[i];
  }
}

void CColoredTracks::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CColorImage::load(ar);

  ar >> m_nVarColorMap;

  if(nVersion >= 1)
  {
    int nCount;
    ar >> nCount;
    for(UINT i = 0; i < nCount; i++)
    {
      if(i >= tcmCount)
        continue;

      ar >> m_vUserDefMin[i];
      ar >> m_vUserDefMax[i];
    }

    if(m_bUserDefRange && (m_nVarColorMap < tcmCount))
      restore_user_range();
  }
} 

//---------------------------------------------------------------------------------------
// CColoredCrossSection - a class for output of the ion beam cross-section colored by 
//                        initial radius. The data could be taken and plotted by Origin.
//---------------------------------------------------------------------------------------
CColoredCrossSection::CColoredCrossSection()
{
  set_default();
}

CColoredCrossSection::~CColoredCrossSection()
{
  clear();
}


void CColoredCrossSection::draw()
{
  if(!m_bReady)
    build_draw_array();
}

void CColoredCrossSection::restore_user_range()
{
}

const char* CColoredCrossSection::get_var_name(int nVar) const
{
  switch(nVar)
  {
    case varIonTemp: return _T("Ion Temperature");
    case varStartRadius: return _T("Ion Start Radius");
  }

  return _T("None");
}

void CColoredCrossSection::get_min_max()
{
// Temporarly:
  if(m_nVar == varIonTemp)
  {
    m_fMinVal = 190;
    m_fMaxVal = 1600;
    return;
  }

  double fVal;
  m_fMinVal = FLT_MAX;
  m_fMaxVal = -FLT_MAX;

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CTrackVector& vTracks = pObj->get_tracks();

  CTrackItem item;
  size_t nTracksCount = vTracks.size();
  for(size_t i = 0; i < nTracksCount; i++)
  {
    const CTrack& track = vTracks.at(i);
    size_t nItemsCount = track.size();
    for(size_t j = 0; j < nItemsCount; j++)
    {
      track.get_track_item(j, item);
      fVal = get_vert_value(item, i);
      if(fVal < m_fMinVal)
        m_fMinVal = fVal;
      if(fVal > m_fMaxVal)
        m_fMaxVal = fVal;
    }
  }
}

void CColoredCrossSection::set_default()
{
  CColorImage::set_default();

  m_fCrossSectX = 0;
  m_nVar = varIonTemp;
}

bool CColoredCrossSection::build_draw_array()
{
  clear();
  get_min_max();
  get_values_array();
  if(m_vValues.size() == 0)
    return false;

  set_job_name("Calculating Cross-Section...");

  get_colors_array();
  collect_intersections();
  return true;
}

void CColoredCrossSection::collect_intersections()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CTrackVector& vTracks = pObj->get_tracks();

  double fKsi, fdX, fVal;
  Vector3D vPoint(m_fCrossSectX, 0, 0);

  CTrackItem item;
  size_t nTrackCount = vTracks.size();
  for(size_t i = 0; i < nTrackCount; i++)
  {
    set_progress(100 * i / nTrackCount);

    const CTrack& track = vTracks.at(i);
    size_t nItemCount = track.size();

    Vector3D vPrev = track.at(0)->pos, vCurr;
    for(size_t j = 1; j < nItemCount; j++)
    {
      vCurr = track.at(j)->pos;
      if((vCurr.x >= m_fCrossSectX) && (vPrev.x < m_fCrossSectX) || (vCurr.x < m_fCrossSectX) && (vPrev.x >= m_fCrossSectX))
      {
        fdX = vCurr.x - vPrev.x;
        if(fabs(fdX) < Const_Almost_Zero)
          continue;

        fKsi = (m_fCrossSectX - vPrev.x) / fdX;
        vPoint.y = vPrev.y + fKsi * (vCurr.y - vPrev.y);
        vPoint.z = vPrev.z + fKsi * (vCurr.z - vPrev.z);
        m_vPoints.push_back(vPoint);

        track.get_track_item(j, item);
        fVal = get_vert_value(item, i);
        m_vCrossSectColors.push_back(get_color(fVal));
        m_vCrossSectValues.push_back(fVal);
        break;
      }

      vPrev = vCurr;
    }
  }
}

double CColoredCrossSection::get_vert_value(const CTrackItem& item, size_t nTrackIndex) const
{
  switch(m_nVar)
  {
    case varIonTemp: return item.temp;
    case varStartRadius: return CColoredTracks::get_start_radius(nTrackIndex);
  }

  return 0;
}

void CColoredCrossSection::get_colored_point(size_t nIndex, Vector3D& vPoint, double& fVal, RGB_Color& clr) const
{
  if(nIndex >= m_vPoints.size())
    return;

  vPoint = m_vPoints.at(nIndex);
  clr = m_vCrossSectColors.at(nIndex);
  fVal = m_vCrossSectValues.at(nIndex);
}

void CColoredCrossSection::clear()
{
  m_vPoints.clear();
  m_vCrossSectValues.clear();
  m_vCrossSectColors.clear();
}

void CColoredCrossSection::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CColorImage::save(ar);

  ar << m_fCrossSectX;
  ar << m_nVar;

  for(int i = 0; i < varCount; i++)
  {
    ar << m_vUserDefMin[i];
    ar << m_vUserDefMax[i];
  }
}

void CColoredCrossSection::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CColorImage::load(ar);

  ar >> m_fCrossSectX;
  ar >> m_nVar;

  for(int i = 0; i < varCount; i++)
  {
    ar >> m_vUserDefMin[i];
    ar >> m_vUserDefMax[i];
  }
}

};  // namespace EvaporatingParticle
