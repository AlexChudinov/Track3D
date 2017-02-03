
// ParticleTrackingDoc.cpp : implementation of the CParticleTrackingDoc class
//

#include "stdafx.h"
#include "ParticleTracking.h"

#include "ParticleTrackingDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CParticleTrackingDoc

IMPLEMENT_DYNCREATE(CParticleTrackingDoc, CDocument)

BEGIN_MESSAGE_MAP(CParticleTrackingDoc, CDocument)
END_MESSAGE_MAP()


// CParticleTrackingDoc construction/destruction

CParticleTrackingDoc::CParticleTrackingDoc()
{
	// TODO: add one-time construction code here

}

CParticleTrackingDoc::~CParticleTrackingDoc()
{
}

BOOL CParticleTrackingDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}




// CParticleTrackingDoc serialization

void CParticleTrackingDoc::Serialize(CArchive& ar)
{
  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();

	if(ar.IsStoring())
	{
    pDrawObj->save(ar);
    pObj->save(ar);
	}
	else
	{
    const char* pFName = (const char*)ar.m_strFileName;
    std::string sFileName(pFName);
    pDrawObj->load(ar);
    pObj->load(ar); // CTracker::load() must be the last as CTracker::update_interface() is called from it.
	}
}


// CParticleTrackingDoc diagnostics

#ifdef _DEBUG
void CParticleTrackingDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CParticleTrackingDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// CParticleTrackingDoc commands
