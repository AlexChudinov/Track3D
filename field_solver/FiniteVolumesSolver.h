#pragma once

#include "AnsysMesh.h"
#include "DirichletTesselation.h"


namespace EvaporatingParticle
{

typedef std::vector<float> CFloatArray;

struct CNodeInfo
{
  CNodeInfo(unsigned char t = 0, float v = 0)
    : nType(t), fValue(v)
  {
  }

  unsigned char nType;    // 0 - inner node, 1 - BC of the first type, 2 - zero gradient.
  float         fValue;   // value in the node.
};

struct CSolutionInfo
{
  CSolutionInfo(UINT nCount, bool b = false)
    : bSuccess(b)
  {
    vMaxErrHist.assign(nCount, 0);
  }

  bool          bSuccess;
  CFloatArray   vMaxErrHist;
};

typedef std::vector<CNodeInfo> CNodeInfoArray;
//-------------------------------------------------------------------------------------------------
// CFiniteVolumesSolver - implementation of the Laplace equation solver by the finite volumes method.
//-------------------------------------------------------------------------------------------------
class CFiniteVolumesSolver : public CObject
{
public:
  CFiniteVolumesSolver(CAnsysMesh* pMesh, CDirichletTesselation* pTess);
  ~CFiniteVolumesSolver();

// Use this function to set non-trivial boundary conditions of the 1-st type, like Coulomb potential:
  void                    set_boundary_conditions(const CIndexVector& vNodeIds, const CFloatArray& vVals);
// Boundary conditions of the 1-st type with a constant value on the whole region:
  void                    set_boundary_conditions(const CIndexVector& vNodeIds, float fVal);
// Boundary conditions of the 2-nd type:
  void                    set_boundary_conditions(const CIndexVector& vNodeIds);

  CSolutionInfo           solve(float fTol, UINT nTimeStepsCount, CFloatArray& vRes, bool bMultiThread = true);

protected:
  void                    init();
  void                    clear();

  float                   one_iter();   // one iteration in solving the system of linear equations.
// Multithreading support:
  float                   single_node_iter(CNode3D* pNode, CDirichletCell* pCell, const CFloatArray& vU0, const CNodesCollection& vNodes, float fDiagCoeff, int nType);

  void                    prepare();
  void                    calc_diag();  // calculation of diagonal components of the system matrix.

// Returns maximal relative discrepancy between U(r) at current (m_vRes) and previous (m_vU0) iterations.
  float                   get_max_relative_error();

private:
  CAnsysMesh*             m_pMesh;
  CDirichletTesselation*  m_pTess;

// Values of U(r) for every iteration. Contains also initial and boundary conditions.
  CNodeInfoArray          m_vRes;

// Run-time:
  CFloatArray             m_vU0;        // values of U(r) at the previous iteration, updated at the end of every iteration in get_max_relative_error().
  CFloatArray             m_vDiagCoeff; // diagonal coefficients of the system matrix;
};

};  // namespace EvaporatingParticle
