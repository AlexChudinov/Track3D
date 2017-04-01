
#include "stdafx.h"
#include "TrackItem.h"
//#include "Integrators.hpp"
//#include "BarnesHut.h"
//#include "StochasticGasDynamic.hpp"
//#include "Tracker.hpp"
//#include "ParticleTracking.h"
//#include <algorithm>


namespace EvaporatingParticle
{

CBaseTrackItem* CBaseTrackItem::create(int nType)
{
  switch(nType)
  {
    case CTrack::ptDroplet:
    {
      CDropletTrackItem* pDropletItem = new CDropletTrackItem();
      return (CBaseTrackItem*)pDropletItem;
    }
    case CTrack::ptIon:
    {
      CIonTrackItem* pIonItem = new CIonTrackItem();
      return (CBaseTrackItem*)pIonItem;
    }
  }

  return NULL;
}

//---------------------------------------------------------------------------------------
//
// CDropletTrackItem:
// pState is an array of size DROPLET_STATE_SIZE = 8.
// pItem->pos(pState[0], pState[1], pState[2]);
// pItem->vel(pState[3], pState[4], pState[5]);
// pItem->temp = pState[6];
// pItem->mass = pState[7];
//
// CIonTrackItem:
// pState is an array of size ION_STATE_SIZE = 8.
// pItem->pos(pState[0], pState[1], pState[2]);
// pItem->vel(pState[3], pState[4], pState[5]);
// pItem->temp = pState[6];
// pItem->unfragm = pState[7];
//
//---------------------------------------------------------------------------------------
CBaseTrackItem* CBaseTrackItem::create(int nType, int nElemId, double fTime, double* pState)
{
  if(pState == NULL)
    return NULL;

  Vector3D vPos(pState[0], pState[1], pState[2]);
  Vector3D vVel(pState[3], pState[4], pState[5]);
  switch(nType)
  {
    case CTrack::ptDroplet:
    {
      double fTemp = pState[6];
      double fMass = pState[7];
      CDropletTrackItem* pDropletItem = new CDropletTrackItem(nElemId, vPos, vVel, fTemp, fMass, fTime);

      return (CBaseTrackItem*)pDropletItem;
    }
    case CTrack::ptIon:
    {
      double fIonTemp = pState[6];
      double fUnfragm = pState[7];
      CIonTrackItem* pIonItem = new CIonTrackItem(nElemId, vPos, vVel, fIonTemp, fUnfragm, fTime);

      return (CBaseTrackItem*)pIonItem;
    }
  }
  
  return NULL;  
}

void CBaseTrackItem::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  ar << nElemId;
  ar << pos.x;
  ar << pos.y;
  ar << pos.z;
  ar << vel.x;
  ar << vel.y;
  ar << vel.z;
  ar << time;
}

void CBaseTrackItem::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> nElemId;
  ar >> pos.x;
  ar >> pos.y;
  ar >> pos.z;
  ar >> vel.x;
  ar >> vel.y;
  ar >> vel.z;
  ar >> time;
}

//[AC 31/03/2017] Memory manager
void CBaseTrackItem::operator delete(void * ptr, size_t n)
{
	((CBaseTrackItem*)ptr)->deleteObj();
}
void CIonTrackItem::deleteObj()
{
	BlockPool<CIonTrackItem>::getInstance().freeBlock(this);
}
void CDropletTrackItem::deleteObj()
{
	BlockPool<CDropletTrackItem>::getInstance().freeBlock(this);
}
//[/AC]

//---------------------------------------------------------------------------------------
//  CIonTrackItem
//---------------------------------------------------------------------------------------
void CIonTrackItem::state(double* pState) const
{
  pState[0] = pos.x;
  pState[1] = pos.y;
  pState[2] = pos.z;
  pState[3] = vel.x;
  pState[4] = vel.y;
  pState[5] = vel.z;
  pState[6] = temp;
  pState[7] = unfragm;
}

CBaseTrackItem* CIonTrackItem::copy() const
{
  CIonTrackItem* pItem = new CIonTrackItem(nElemId, pos, vel, temp, unfragm, time);
  pItem->tempinf = tempinf;
  return (CBaseTrackItem*)pItem;
}

double CIonTrackItem::get_mass() const
{
  return 0.;
}

double CIonTrackItem::get_temp() const
{
  return temp;
}

void CIonTrackItem::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CBaseTrackItem::save(ar);

  ar << temp;
  ar << tempinf;
  ar << unfragm;
}

void CIonTrackItem::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CBaseTrackItem::load(ar);

  ar >> temp;
  ar >> tempinf;
  ar >> unfragm;
}

//---------------------------------------------------------------------------------------
//  CDropletTrackItem
//---------------------------------------------------------------------------------------
void CDropletTrackItem::state(double* pState) const
{
  pState[0] = pos.x;
  pState[1] = pos.y;
  pState[2] = pos.z;
  pState[3] = vel.x;
  pState[4] = vel.y;
  pState[5] = vel.z;
  pState[6] = temp;
  pState[7] = mass;
}

CBaseTrackItem* CDropletTrackItem::copy() const
{
  CDropletTrackItem* pItem = new CDropletTrackItem(nElemId, pos, vel, temp, mass, time);
  return (CBaseTrackItem*)pItem;
}

double CDropletTrackItem::get_mass() const
{
  return mass;
}

double CDropletTrackItem::get_temp() const
{
  return temp;
}

void CDropletTrackItem::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CBaseTrackItem::save(ar);

  ar << temp;
  ar << mass;
}

void CDropletTrackItem::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CBaseTrackItem::load(ar);

  ar >> temp;
  ar >> mass;
}

//---------------------------------------------------------------------------------------
// CTrack
//---------------------------------------------------------------------------------------
void CTrack::get_track_item(size_t nInd, CTrackItem& item) const
{
  CBaseTrackItem* pItem = at(nInd);
  item.nElemId = pItem->nElemId;
  item.pos = pItem->pos;
  item.vel = pItem->vel;
  item.time = pItem->time;

  switch(m_nType)
  {
    case ptDroplet:
    {
      CDropletTrackItem* pDropletItem = (CDropletTrackItem*)pItem;
      item.temp = pDropletItem->temp;
      item.mass = pDropletItem->mass;
      item.tempinf = 0;
      item.unfragm = 1;
    }
    case ptIon:
    {
      CIonTrackItem* pIonItem = (CIonTrackItem*)pItem;
      item.temp = pIonItem->temp;
      item.tempinf = pIonItem->tempinf;
      item.unfragm = pIonItem->unfragm;
      item.mass = 0;
    }
  }
}

};  // namespace EvaporatingParticle