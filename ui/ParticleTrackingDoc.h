
// ParticleTrackingDoc.h : interface of the CParticleTrackingDoc class
//

#pragma once

class CParticleTrackingDoc : public CDocument
{
protected: // create from serialization only
	CParticleTrackingDoc();
	DECLARE_DYNCREATE(CParticleTrackingDoc)

// Attributes
public:

// Operations
public:

// Overrides
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);

// Implementation
public:
	virtual ~CParticleTrackingDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
};


