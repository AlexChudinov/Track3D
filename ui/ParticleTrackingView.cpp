
// ParticleTrackingView.cpp : implementation of the CParticleTrackingView class
//

#include "stdafx.h"
#include "ParticleTracking.h"

#include "ParticleTrackingDoc.h"
#include "ParticleTrackingView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CParticleTrackingView

IMPLEMENT_DYNCREATE(CParticleTrackingView, CView)

BEGIN_MESSAGE_MAP(CParticleTrackingView, CView)
  ON_WM_LBUTTONDOWN()
  ON_WM_LBUTTONUP()
  ON_WM_MOUSEMOVE()
  ON_WM_RBUTTONUP()
  ON_WM_MOUSEWHEEL()
END_MESSAGE_MAP()

// CParticleTrackingView construction/destruction

CParticleTrackingView::CParticleTrackingView()
{
	// TODO: add construction code here

}

CParticleTrackingView::~CParticleTrackingView()
{
}

BOOL CParticleTrackingView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

// CParticleTrackingView drawing

void CParticleTrackingView::OnDraw(CDC* pDC)
{
	CParticleTrackingDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if(!pDoc)
		return;

  EvaporatingParticle::CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  if(pDC != NULL)
    pObj->update_interface();

  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  
  if(pDrawObj->get_window_handle() == NULL)
    pDrawObj->set_window_handle(GetSafeHwnd());

  if(pDrawObj->get_tracker() == NULL)
    pDrawObj->set_tracker(pObj);

  pDrawObj->draw();
}

// Mouse events

void CParticleTrackingView::OnLButtonDown(UINT nFlags, CPoint point)
{
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->on_left_button_down(point);
}

void CParticleTrackingView::OnLButtonUp(UINT nFlags, CPoint point)
{
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->on_left_button_up(point);
}

void CParticleTrackingView::OnMouseMove(UINT nFlags, CPoint point)
{
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->on_mouse_move(point);
//  OnDraw(NULL);
}

BOOL CParticleTrackingView::OnMouseWheel(UINT nFlags, short zDelta, CPoint point)
{
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->on_mouse_wheel(zDelta, point);
  OnDraw(NULL);

  return CView::OnMouseWheel(nFlags, zDelta, point);
}

void CParticleTrackingView::OnRButtonUp(UINT nFlags, CPoint point)
{
//	ClientToScreen(&point);
//	OnContextMenu(this, point);
  EvaporatingParticle::CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->on_right_button_up(point);
  OnDraw(NULL);
}

void CParticleTrackingView::OnContextMenu(CWnd* pWnd, CPoint point)
{
	theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
}


// CParticleTrackingView diagnostics

#ifdef _DEBUG
void CParticleTrackingView::AssertValid() const
{
	CView::AssertValid();
}

void CParticleTrackingView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CParticleTrackingDoc* CParticleTrackingView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CParticleTrackingDoc)));
	return (CParticleTrackingDoc*)m_pDocument;
}
#endif //_DEBUG


// CParticleTrackingView message handlers
