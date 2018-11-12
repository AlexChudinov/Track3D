
// ParticleTracking.h : main header file for the ParticleTracking application
//
#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"       // main symbols

#include "ExportOpenFOAM.h"
#include "Calculator.h"
#include "ElectricField.h"
#include "DirichletTesselation.h"
#include "DomainCrossSection.h"
#include "DrawTrack.h"

// CParticleTrackingApp:
// See ParticleTracking.cpp for the implementation of this class
//

class CParticleTrackingApp : public CWinAppEx
{
public:
	CParticleTrackingApp();

  static CParticleTrackingApp*                Get();

  EvaporatingParticle::CTracker*              GetTracker();
  EvaporatingParticle::CDirichletTesselation* GetDirichletTess();
  EvaporatingParticle::CTrackDraw*            GetDrawObj();
  EvaporatingParticle::CExportOpenFOAM*       GetExporter();
  EvaporatingParticle::CCalcCollection*       GetCalcs();
  EvaporatingParticle::CFieldDataColl*        GetFields();
  EvaporatingParticle::CCrossSectColl*        GetPlanes();

  bool                                        is_terminated() const { return m_bTerminate; }
  void                                        terminate(bool bTerm) { m_bTerminate = bTerm; }

  void                                        SelectedRegionChanged(EvaporatingParticle::CNamesVector* pRegNames);

// Overrides
public:
	virtual                                     BOOL InitInstance();
  virtual BOOL                                OnIdle(LONG lCount); // return TRUE if more idle processing.

// Implementation
	UINT                                        m_nAppLook;
	BOOL                                        m_bHiColorIcons;

	virtual void                                PreLoadState();
	virtual void                                LoadCustomState();
	virtual void                                SaveCustomState();

  virtual BOOL                                PreTranslateMessage(MSG* pMsg); // for key press events.

	afx_msg void                                OnAppAbout();

private:
  EvaporatingParticle::CTracker               m_Tracker;
  EvaporatingParticle::CDirichletTesselation  m_DirichletTess;
  EvaporatingParticle::CTrackDraw             m_Drawer;
  EvaporatingParticle::CExportOpenFOAM        m_Exporter;
  EvaporatingParticle::CCalcCollection        m_vCalcs;
  EvaporatingParticle::CFieldDataColl         m_vFields;
  EvaporatingParticle::CCrossSectColl         m_vPlanes;

  bool                                        m_bTerminate;

	DECLARE_MESSAGE_MAP()
};

extern CParticleTrackingApp theApp;
