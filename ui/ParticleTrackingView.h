
// ParticleTrackingView.h : interface of the CParticleTrackingView class
//


#pragma once


class CParticleTrackingView : public CView
{
protected: // create from serialization only
	CParticleTrackingView();
	DECLARE_DYNCREATE(CParticleTrackingView)

// Attributes
public:
	CParticleTrackingDoc* GetDocument() const;

// Operations
public:

// Overrides
public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:

// Implementation
public:
	virtual ~CParticleTrackingView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
  afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
  afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
  afx_msg void OnMouseMove(UINT nFlags, CPoint point);

  afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint point);

//	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
  afx_msg void OnLButtonDblClk(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in ParticleTrackingView.cpp
inline CParticleTrackingDoc* CParticleTrackingView::GetDocument() const
   { return reinterpret_cast<CParticleTrackingDoc*>(m_pDocument); }
#endif

