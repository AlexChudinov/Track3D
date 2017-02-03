#ifndef _InterpolationTest_
#define _InterpolationTest_

#include "Elements.h"

namespace EvaporatingParticle
{

class CInterpolationTest
{
public:
  CInterpolationTest();
  virtual ~CInterpolationTest();

  void            do_test();
  void            create_objects();
  void            set_default();

  void            test_func();
  void            test_param_at_nodes();
  void            test_inverse();

  void            test_inside_outside();

  void            apply_rotation(CElem3D* pElem, const Matrix3D& m);

private:
  CTetra*         m_pTetra;
  CPyramid*       m_pPyramid;
  CWedge*         m_pWedge;
  CHexa*          m_pHexa;

  Vector3D        m_vTetraVert[4];
  Vector3D        m_vPyramidVert[5];
  Vector3D        m_vWedgeVert[6];
  Vector3D        m_vHexaVert[8];
};

};  // end of namespace EvaporatingParticle

#endif // #ifndef _InterpolationTest_