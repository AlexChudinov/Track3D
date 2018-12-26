
#include "stdafx.h"

#include "TrajectSelector.h"

#include "float.h"
#include "constant.hpp"

#include <GL/gl.h>
#include <gl/glu.h>

#include "../utilities/ParallelFor.h"


namespace EvaporatingParticle
{

CTrajectSelector::CTrajectSelector(const CTrackVector& vTracks, int x0, int y0)
  : m_vTracks(vTracks), m_vCur(x0, y0)
{
  init(x0, y0);
}

void CTrajectSelector::init(int x0, int y0)
{
  m_nTrajectId = -1;

  m_nR = 2;
  m_nXmin = x0 - m_nR;
  m_nXmax = x0 + m_nR;
  m_nYmin = y0 - m_nR;
  m_nYmax = y0 + m_nR;
  m_nR2 = m_nR * m_nR;

  glGetDoublev(GL_MODELVIEW_MATRIX, m_pModelMtx);
  glGetDoublev(GL_PROJECTION_MATRIX, m_pProjMtx);
  glGetIntegerv(GL_VIEWPORT, m_pViewPort);
}

int CTrajectSelector::find_traject()
{
  size_t nTracksCount = m_vTracks.size();
  ThreadPool::splitInPar(nTracksCount,
	[&](size_t i) 
  {
    const CTrack& track = m_vTracks.at(i);
    if(track_under_cursor(track))
    {
      m_nTrajectId = (int)i;
      return;
    }
  });

  return m_nTrajectId;
}

bool CTrajectSelector::track_under_cursor(const CTrack& track) const
{
  GLint nRes = 0;
  size_t nItemsCount = track.size();
  Vector3D vA, vB;  // a single track edge AB.
  double xA, yA, zA, xB, yB, zB, alpha, sqlenAB, r2;
  Vector2D v2dAO, v2dAB;
  for(size_t i = 1; i < nItemsCount; i++)
  {
    vA = track.at(i - 1)->pos;
    nRes = gluProject(vA.x, vA.y, vA.z, m_pModelMtx, m_pProjMtx, m_pViewPort, &xA, &yA, &zA);
    if(nRes == GLU_FALSE)
      continue;

    vB = track.at(i)->pos;
    nRes = gluProject(vB.x, vB.y, vB.z, m_pModelMtx, m_pProjMtx, m_pViewPort, &xB, &yB, &zB);
    if(nRes == GLU_FALSE)
      continue;

    if(((xA < m_nXmin) && (xB < m_nXmin)) || ((xA > m_nXmax) && (xB > m_nXmax)))
      continue;

    yA = (double)(m_pViewPort[3]) - yA;
    yB = (double)(m_pViewPort[3]) - yB;
    if(((yA < m_nYmin) && (yB < m_nYmin)) || ((yA > m_nYmax) && (yB > m_nYmax)))
      continue;

    if((xA >= m_nXmin) && (xA <= m_nXmax) && (xB >= m_nXmin) && (xB <= m_nXmax))
      return true;

    if((yA >= m_nYmin) && (yA <= m_nYmax) && (yB >= m_nYmin) && (yB <= m_nYmax))
      return true;

    v2dAO = Vector2D(m_vCur.x - xA, m_vCur.y - yA);
    v2dAB = Vector2D(xB - xA, yB - yA);
    sqlenAB = v2dAB.sqlength();
    if(sqlenAB < Const_Almost_Zero)
      continue;

    alpha = (v2dAO & v2dAB) / sqlenAB;
    r2 = (alpha * v2dAB - v2dAO).sqlength();
    if(r2 < m_nR2)
      return true;
  }

  return false;
}


}; // namespace EvaporatingParticle
