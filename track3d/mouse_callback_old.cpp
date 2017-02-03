/*
void CTrackDraw::mouse_callback(int nEvent, int nXpos, int nYpos, int nFlags, void* pData)
{
  CTrackDraw* pObj = (CTrackDraw*)pData;
  switch(nEvent)
  {
    case CV_EVENT_MOUSEMOVE:
    {
      Vector2D vPos = pObj->pix_2_world(cvPoint(nXpos, nYpos));
      std::string strX = "X = " + pObj->dbl_to_str(10 * vPos.x);
      std::string strY = "Y = " + pObj->dbl_to_str(10 * vPos.y);

      cvRectangle(pObj->m_pImg, cvPoint(0, 0), cvPoint(250, 30), cvScalar(0, 0, 0, 0), CV_FILLED);
      cvPutText(pObj->m_pImg, strX.c_str(), cvPoint(10, 20), &(pObj->m_Font), sClrArray[10]);
      cvPutText(pObj->m_pImg, strY.c_str(), cvPoint(130, 20), &(pObj->m_Font), sClrArray[10]);
      cvShowImage(pObj->window_caption(), pObj->m_pImg);
      break;
    }
    case CV_EVENT_LBUTTONDOWN:
    {
      pObj->m_nRegime = nFlags & CV_EVENT_FLAG_SHIFTKEY ? nRegimeScale : nRegimeMove;
      pObj->m_vStart = cvPoint(nXpos, nYpos);
      break;
    }
    case CV_EVENT_LBUTTONUP:
    {
      CvPoint shift = cvPoint(nXpos - pObj->m_vStart.x, nYpos - pObj->m_vStart.y);
      if(pObj->m_nRegime == nRegimeMove)
      {
        double dx = -pObj->m_fScale * shift.x;
        double dy = pObj->m_fScale * shift.y;

        pObj->m_fXMin += dx;
        pObj->m_fXMax += dx;
        pObj->m_fYMin += dy;
        pObj->m_fYMax += dy;
      }
      else if(pObj->m_nRegime == nRegimeScale)
      {
        Vector2D vCenter = pObj->pix_2_world(cvPoint(pObj->m_nXRes / 2, pObj->m_nYRes / 2));
        double fScale = 1. + 4. * fabs((double)shift.y) / pObj->m_nYRes;
        if(shift.y >= 0)
        {
          pObj->m_fXRange *= fScale;
          pObj->m_fYRange *= fScale;
          pObj->m_fScale *= fScale;
          pObj->m_fInvScale = 1. / pObj->m_fScale;
        }
        else
        {
          pObj->m_fXRange /= fScale;
          pObj->m_fYRange /= fScale;
          pObj->m_fScale /= fScale;
          pObj->m_fInvScale = 1. / pObj->m_fScale;
        }

        pObj->m_fXMin = vCenter.x - 0.5 * pObj->m_fXRange;
        pObj->m_fXMax = vCenter.x + 0.5 * pObj->m_fXRange;
        pObj->m_fYMin = vCenter.y - 0.5 * pObj->m_fYRange;
        pObj->m_fYMax = vCenter.y + 0.5 * pObj->m_fYRange;
      }

      pObj->m_nRegime = nRegimeNone;
      pObj->draw();
      break;
    }
    case CV_EVENT_RBUTTONUP:
    {
      pObj->default_ranges();
      pObj->draw();
      break;
    }
  }
}
*/