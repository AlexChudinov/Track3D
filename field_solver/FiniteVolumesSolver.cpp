
#include "stdafx.h"


#include "../field_solver/FiniteVolumesSolver.h"
#include "../utilities/ParallelFor.h"

// DEBUG
#include "ParticleTracking.h"
// END DEBUG


namespace EvaporatingParticle
{

CFiniteVolumesSolver::CFiniteVolumesSolver(CAnsysMesh* pMesh, CDirichletTesselation* pTess)
  : m_pMesh(pMesh), m_pTess(pTess)
{
  init();
}

CFiniteVolumesSolver::~CFiniteVolumesSolver()
{
  clear();
}

void CFiniteVolumesSolver::init()
{
  m_vRes.clear();
  m_vRes.assign(m_pMesh->get_nodes().size(), CNodeInfo());
}

void CFiniteVolumesSolver::clear()
{
  m_vRes.clear();
  m_vU0.clear();
}

void CFiniteVolumesSolver::set_boundary_conditions(const CIndexVector& vNodeIds, const CFloatArray& vVals)
{
  size_t nCount = vNodeIds.size();
  if(nCount == 0 || nCount != vVals.size())
    return;

  float fVal;
  UINT nNodeId;
  for(size_t i = 0; i < nCount; i++)
  {
    nNodeId = vNodeIds.at(i);
    CNodeInfo& info = m_vRes.at(nNodeId);
    info.fValue = vVals.at(i);
    info.nType = 1;
  }
}

void CFiniteVolumesSolver::set_boundary_conditions(const CIndexVector& vNodeIds, float fVal)
{
  size_t nCount = vNodeIds.size();
  if(nCount == 0)
    return;

  UINT nNodeId;
  for(size_t i = 0; i < nCount; i++)
  {
    nNodeId = vNodeIds.at(i);
    CNodeInfo& info = m_vRes.at(nNodeId);
    info.fValue = fVal;
    info.nType = 1;
  }
}

void CFiniteVolumesSolver::set_boundary_conditions(const CIndexVector& vNodeIds)
{
  size_t nCount = vNodeIds.size();
  if(nCount == 0)
    return;

  UINT nNodeId;
  for(size_t i = 0; i < nCount; i++)
  {
    nNodeId = vNodeIds.at(i);
    CNodeInfo& info = m_vRes.at(nNodeId);
    if(info.nType == 1)
      continue; // do NOT change the type from 1 to 2.

    info.nType = 2;
  }
}

CSolutionInfo CFiniteVolumesSolver::solve(float fTol, UINT nIterCount, CFloatArray& vRes, bool bMultiThread)
{
  set_job_name("Solving...");
  set_progress(0);

  prepare();

  CSolutionInfo info(nIterCount);

  int nType;
  float fMaxRelErr = FLT_MAX;
  const CNodesCollection& vNodes = m_pMesh->get_nodes();
  size_t nNodesCount = vNodes.size();
  vRes.resize(nNodesCount);
  vRes = m_vU0;

  if(bMultiThread)
  {
    for(UINT k = 0; k < nIterCount; k++)
    {
		set_progress(100 * (k + 1) / nIterCount);

      ThreadPool::splitInPar(vNodes.size(),
	      [&](size_t i) 
        {
		  CNode3D* pNode = vNodes[i];
		  CDirichletCell* pCell = m_pTess->get_cell(i);
          float fDiagCoeff = m_vDiagCoeff.at(i);
          int nType = m_vRes.at(i).nType;
          if(nType != 1)
	          m_vRes[i].fValue = single_node_iter(pNode, pCell, m_vU0, vNodes, fDiagCoeff, nType);
        });

      fMaxRelErr = get_max_relative_error();

      info.vMaxErrHist.at(k) = fMaxRelErr;
      if(fMaxRelErr <= fTol)
      {
        info.bSuccess = true;
        break;
      }
    }
  }
  else
  {
    for(UINT i = 0; i < nIterCount; i++)
    {
		set_progress(100 * (i + 1) / nIterCount);
		for (size_t j = 0; j < vNodes.size(); ++j)
		{
			CNode3D* pNode = vNodes[j];
			CDirichletCell * pCell = m_pTess->get_cell(j);
			double fDiagCoeff = m_vDiagCoeff.at(j);
			nType = m_vRes.at(j).nType;
			if (nType != 1)
				m_vRes[j].fValue = single_node_iter(pNode, pCell, m_vU0, vNodes, fDiagCoeff, nType);
		}

		fMaxRelErr = get_max_relative_error();

		info.vMaxErrHist.at(i) = fMaxRelErr;
		if (fMaxRelErr <= fTol)
		{
			info.bSuccess = true;
			break;
		}
    }
  }

  for(size_t j = 0; j < nNodesCount; j++)
    vRes[j] = m_vRes.at(j).fValue;

  return info;
}

static const Vector3D scvNull(0, 0, 0);

float CFiniteVolumesSolver::one_iter()
{
  CNode3D* pNode = NULL;
  CDirichletCell* pCell = NULL;
  const CNodesCollection& vNodes = m_pMesh->get_nodes();
  Vector3D vNorm, vNbr, vF, vB;
  Matrix3D mMtx;
  float fSum;

  size_t nNodesCount = vNodes.size(), nNbrCount;
  for(size_t i = 0; i < nNodesCount; i++)
  {
    pNode = vNodes.at(i);
    nNbrCount = pNode->vNbrNodes.size();
    pCell = m_pTess->get_cell(i);

    const CNodeInfo& info = m_vRes.at(i);
    switch(info.nType)
    {
      case 0: // inner node.
      {
        fSum = 0;
        for(size_t j = 0; j < nNbrCount; j++)
          fSum += m_vU0[pNode->vNbrNodes.at(j)] * pCell->pFaceSquare[j] / pCell->pNbrDist[j];

        m_vRes.at(i).fValue = fSum / m_vDiagCoeff.at(i);
        break;
      }
      case 2: // boundary node, zero gradient condition.
      {
        vB = scvNull;
        vNorm = Vector3D(pCell->pFaceSquare[7], pCell->pFaceSquare[8], pCell->pFaceSquare[9]);  // normal at the boundary node.
        mMtx = m_pTess->get_bound_cell_mtx(pCell);
        vF = mMtx * vNorm;
        for(size_t j = 0; j < nNbrCount; j++)
        {
          vNbr = m_pTess->get_neighbor_vector(i, j);  // note: in this function i is the global node index, but j is the local neighbor index of this node.
          vNbr /= (pCell->pNbrDist[j] * pCell->pNbrDist[j]);
          vB += vNbr * m_vU0[pNode->vNbrNodes.at(j)];
        }

        m_vRes.at(i).fValue = (vF & vB) / m_vDiagCoeff.at(i);
        break;
      }
    }
  }

  return get_max_relative_error();
}

float CFiniteVolumesSolver::single_node_iter(CNode3D* pNode, CDirichletCell* pCell, const CFloatArray& vU0, const CNodesCollection& vNodes, float fDiagCoeff, int nType)
{
  size_t nNbr, nNbrCount = pNode->vNbrNodes.size();
  switch(nType)
  {
    case 0:
    {
      float fSum = 0;
      for(size_t j = 0; j < nNbrCount; j++)
        fSum += vU0.at(pNode->vNbrNodes.at(j)) * pCell->pFaceSquare[j] / pCell->pNbrDist[j];

      return fSum / fDiagCoeff;
    }
    case 2:
    {
      Vector3D vB = scvNull, vNbr;
      Vector3D vNorm(pCell->pFaceSquare[7], pCell->pFaceSquare[8], pCell->pFaceSquare[9]);  // normal at the boundary node.

      Matrix3D mMtx(pCell->pFaceSquare[0], pCell->pFaceSquare[1], pCell->pFaceSquare[2],
                    pCell->pFaceSquare[1], pCell->pFaceSquare[3], pCell->pFaceSquare[4],
                    pCell->pFaceSquare[2], pCell->pFaceSquare[4], pCell->pFaceSquare[5]);

      Vector3D vF = mMtx * vNorm;
      for(size_t j = 0; j < nNbrCount; j++)
      {
        nNbr = pNode->vNbrNodes.at(j);
        vNbr = vNodes.at(nNbr)->pos - pNode->pos;
        vNbr /= (pCell->pNbrDist[j] * pCell->pNbrDist[j]);
        vB += vNbr * vU0.at(nNbr);
      }

      return float(vF & vB) / fDiagCoeff;
    }
  }

  return vU0.at(pNode->nInd);
}

void CFiniteVolumesSolver::prepare()
{
  size_t nNodesCount = m_pMesh->get_nodes().size();
  m_vDiagCoeff.assign(nNodesCount, 0);

  calc_boundary_normals();
  calc_diag();

// Initialize run-time arrays:
  m_vU0.assign(nNodesCount, 0);
  for(size_t i = 0; i < nNodesCount; i++)
  {
    CNodeInfo& info = m_vRes.at(i);
    m_vU0[i] = info.fValue;
  }
}

void CFiniteVolumesSolver::calc_diag()
{
  CDirichletCell* pCell = NULL;
  const CNodesCollection& vNodes = m_pMesh->get_nodes();
  Vector3D vNorm, vNbr, vF, vB;
  Matrix3D mMtx;
  float fSum;

  size_t nNodesCount = vNodes.size();
  for(size_t i = 0; i < nNodesCount; i++)
  {
    pCell = m_pTess->get_cell(i);
    const CNodeInfo& info = m_vRes.at(i);
    switch(info.nType)
    {
      case 0:   // inner node.
      {
        fSum = 0;
        for(size_t j = 0; j < pCell->nFaceCount; j++)
          fSum += (pCell->pFaceSquare[j] / pCell->pNbrDist[j]);

        m_vDiagCoeff[i] = fSum;
        break;
      }
      case 1:   // boundary node with 1-st type boundary conditions.
      {
        m_vDiagCoeff[i] = 1;
        break;
      }
      case 2:   // boundary node with 2-nd type boundary conditions.
      {
        vNorm = Vector3D(pCell->pFaceSquare[7], pCell->pFaceSquare[8], pCell->pFaceSquare[9]);  // normal at the boundary node.
        mMtx = m_pTess->get_bound_cell_mtx(pCell);
        vF = mMtx * vNorm;

        vB = scvNull;
        for(size_t j = 0; j < pCell->nFaceCount; j++)
        {
          vNbr = m_pTess->get_neighbor_vector(i, j);
          vNbr /= (pCell->pNbrDist[j] * pCell->pNbrDist[j]);
          vB += vNbr;
        }

        m_vDiagCoeff[i] = vF & vB;
        break;
      }
    }
  }
// DEBUG
  int nIndMin = 0, nIndMax = 0, nIndAbsMin = 0;
  float fDiagCoeffMin = FLT_MAX, fDiagCoeffMax = -FLT_MAX, fAbsCoeffMin = FLT_MAX, fCoeff;
  for(size_t i = 0; i < nNodesCount; i++)
  {
    fCoeff = m_vDiagCoeff[i];
    if(fDiagCoeffMin > fCoeff)
    {
      fDiagCoeffMin = fCoeff;
      nIndMin = i;
    }
    if(fDiagCoeffMax < fCoeff)
    {
      fDiagCoeffMax = fCoeff;
      nIndMax = i;
    }

    fCoeff = fabsf(fCoeff);
    if(fAbsCoeffMin > fCoeff)
    {
      fAbsCoeffMin = fCoeff;
      nIndAbsMin = i;
    }
  }

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  COutputEngine& out_engine = pObj->get_output_engine();
  std::string cPath = out_engine.get_full_path(pObj->get_filename());
  std::string cName("Diag_Coeff_Min_Max");
  std::string cExt(".csv");
  std::string cFileName = cPath + cName + cExt;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if((nErr == 0) && (pStream != NULL))
  {
    fputs("Diag Coeff Min, Diag Coeff Max, Abs Coeff Min\n", pStream);
    fprintf(pStream, "%f,         %f,         %f\n", fDiagCoeffMin, fDiagCoeffMax, fAbsCoeffMin);
    fprintf(pStream, "%d,         %d,         %d\n", nIndMin, nIndMax, nIndAbsMin);

    fclose(pStream);
  }

  cName = "Bound_Nodes";
  cFileName = cPath + cName + cExt;
  nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if((nErr == 0) && (pStream != NULL))
  {
    fputs("Bound Node#, Diag Coeff\n", pStream);
    for(size_t i = 0; i < nNodesCount; i++)
    {
      CNodeInfo& info = m_vRes.at(i);
      if(info.nType == 0)
        continue;

      fprintf(pStream, "  %zd,       %f\n", i, m_vDiagCoeff[i]);
    }

    fclose(pStream);
  }
// END DEBUG
}

void CFiniteVolumesSolver::calc_boundary_normals()
{
  Vector3D vNorm;
  CDirichletCell* pCell = NULL;
  size_t nNodesCount = m_pMesh->get_nodes().size();
  for(size_t i = 0; i < nNodesCount; i++)
  {
    const CNodeInfo& info = m_vRes.at(i);
    if(info.nType != 2)
      continue;

    vNorm = get_bound_norm(i);
    pCell = m_pTess->get_cell(i);
    pCell->pFaceSquare[7] = vNorm.x;
    pCell->pFaceSquare[8] = vNorm.y;
    pCell->pFaceSquare[9] = vNorm.z;
  }
}

Vector3D CFiniteVolumesSolver::get_bound_norm(size_t nNodeId) const
{
  CNode3D* pNode = m_pMesh->get_nodes().at(nNodeId);
  size_t nNbrCount = pNode->vNbrFaces.size();
  if(nNbrCount == 0)
    return scvNull;

  const CRegionsCollection& vRegs = m_pMesh->get_regions();
  size_t nRegCount = vRegs.size();

  UINT nReg, nFace;
  CFace* pFace = NULL;
  Vector3D vNorm(0, 0, 0);
  for(size_t i = 0; i < nNbrCount; i++)
  {
    nReg = pNode->vNbrFaces.at(i).nReg;
    if(nReg >= nRegCount)
      continue;

    nFace = pNode->vNbrFaces.at(i).nFace;
    if(nFace >= vRegs.at(nReg)->vFaces.size())
      continue;

    pFace = vRegs.at(nReg)->vFaces.at(nFace);
// Take into account only those faces (except the first one), in all vertices of which the zero gradient condition is set.
    if((i > 0) && (m_vRes.at(pFace->p0->nInd).nType != 2 || m_vRes.at(pFace->p1->nInd).nType != 2 || m_vRes.at(pFace->p2->nInd).nType != 2))
      continue;

    vNorm += pFace->norm;
  }

  vNorm.normalize();
  return vNorm;
}

float CFiniteVolumesSolver::get_max_relative_error()
{
  float fMaxRelErr = 0, fRelErr, fMaxVal;
  size_t nNodesCount = m_pMesh->get_nodes().size();
  for(size_t i = 0; i < nNodesCount; i++)
  {
    fMaxVal = max(fabsf(m_vRes.at(i).fValue), fabsf(m_vU0[i]));
    if(fMaxVal > Const_Almost_Zero)
    {
      fRelErr = fabsf(m_vRes.at(i).fValue - m_vU0[i]) / fMaxVal;
      if(fMaxRelErr < fRelErr)
        fMaxRelErr = fRelErr;
    }

    m_vU0[i] = m_vRes.at(i).fValue;
  }

  return fMaxRelErr;
}

};  // namespace EvaporatingParticle