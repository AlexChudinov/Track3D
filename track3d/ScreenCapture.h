
#pragma once

namespace EvaporatingParticle
{

class CScreenImage : public CImage
{
public:
  BOOL CaptureRect(const CRect& rect) throw();
  BOOL CaptureWindow(HWND hWnd) throw();
  BOOL CaptureScreen() throw();
};

}; // namespace EvaporatingParticle