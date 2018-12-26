#pragma once

#ifndef _TRAJECT_SELECTOR_
#define _TRAJECT_SELECTOR_

#include "TrackItem.h"
#include "vector2d.hpp"

namespace EvaporatingParticle
{

class CTrajectSelector
{
public:
  CTrajectSelector(const CTrackVector& vTracks, int x0, int y0);

  int                 find_traject();

protected:
  void                init(int x0, int y0);
  bool                track_under_cursor(const CTrack& track) const;

private:
  int                 m_nTrajectId;

  const CTrackVector& m_vTracks;

  Vector2D            m_vCur;

  int                 m_nR,
                      m_nR2;

  int                 m_nXmin,
                      m_nXmax,
                      m_nYmin,
                      m_nYmax;

  double              m_pModelMtx[16];
  double              m_pProjMtx[16];

  int                 m_pViewPort[4];
};

};  // namespace EvaporatingParticle

#endif
