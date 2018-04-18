#pragma once

#include "AnsysMesh.h"
#include "DirichletTesselation.h"


namespace EvaporatingParticle
{

class CFiniteVolumesSolver : public CObject
{
public:
  CFiniteVolumesSolver(CAnsysMesh* pMesh, CDirichletTesselation* pTess);

protected:

private:
  CAnsysMesh*             m_pMesh;
  CDirichletTesselation*  m_pTess;
  std::vector<float>      m_vBoundCond;
};

};  // namespace EvaporatingParticle
