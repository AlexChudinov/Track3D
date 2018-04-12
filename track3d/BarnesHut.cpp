
#include "stdafx.h"

#include "BarnesHut.h"
#include "constant.hpp"
#include "math.h"

namespace EvaporatingParticle
{

//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
OctoTreeCell::OctoTreeCell(const Vector3D& bc, double s)
  : box_center(bc), size(s), ndep(0)
{
  for(int i = 0; i < 8; i++)
    pChild[i] = NULL;

  Qxx = Qyy = Qzz = Qxy = Qyz = Qxz = 0;
}

OctoTreeCell::~OctoTreeCell()
{
  for(int i = 0; i < 8; i++)
  {
    if(pChild[i] != NULL)
      delete pChild[i];
  }
}

//---------------------------------------------------------------------------
//  CBarnesHut - this class approximately computes Coulomb interactions 
//  between charged particles in an ion cloud. 
//---------------------------------------------------------------------------
CBarnesHut::CBarnesHut()
  : m_pTree(NULL),
    m_fDistCoeff(1.5),
    m_bEnableQuadTerms(false),
    m_bReady(false)
{
  set_crit_radius(0.001);  // default value of m_fCritRadius is 0.001 cm.
  set_max_rec_depth(12);
  set_sym_type(symXYonly);
}

CBarnesHut::~CBarnesHut()
{
  if(m_pTree)
  {
    for(int i = m_pTree->charges.size() - 1; i >= 0; i--)
      delete m_pTree->charges.at(i);

    delete m_pTree;
  }

  m_vCells.clear();
}

void CBarnesHut::add_particle(const Vector3D& pos, double q)
{
  OctoTreeCell* pCell = m_pTree;
  ChargeData* pData = new ChargeData(pos, q);
  pCell->charges.push_back(pData);
  while((pCell->charges.size() > 1) && (pCell->ndep < m_nMaxRecDepth))
  {
// There can be two cases. 
    if(pCell->charges.size() > 2)  // 1. There are more than two particles in the cell.
    {
      int i = child_cell_index(pos - pCell->box_center);
      if(pCell->pChild[i] != NULL)
      {
// Add the particle to the collection of the daughter cell.
        pCell->pChild[i]->charges.push_back(pData);
// Re-assign pCell and repeat the body of "while" again.
        pCell = pCell->pChild[i];
      }
      else
      {
// Create a new cell and settle the only particle there.
        pCell->pChild[i] = create_child_cell(i, pCell);
        pCell->pChild[i]->charges.push_back(pData);
        break;
      }
    }
    else  // 2. There are only two particles in the cell. In this case the children are not initialized yet.
    {
// There are only two particles and child cells are not initialized. The cell must be sub-divided
// into 8 sub-cells and BOTH particles must be settled in proper cells. Again, there can be two cases.
      int i = child_cell_index(pCell->charges.at(0)->pos - pCell->box_center);
      int j = child_cell_index(pCell->charges.at(1)->pos - pCell->box_center);
      if(i != j)
      {
// 1. The daughter cells for the two particles are different. After creating cells and settling
//    particles (one in each sub-cell) break the "while" cycle.
        pCell->pChild[i] = create_child_cell(i, pCell);
        pCell->pChild[i]->charges.push_back(pCell->charges.at(0));
        pCell->pChild[j] = create_child_cell(j, pCell);
        pCell->pChild[j]->charges.push_back(pCell->charges.at(1));
        break;
      }
      else
      {
// 2. Both particles come upon the same daughter cell. Create this cell, put both particles there.
//    Re-assign pCell and repeat the body of "while" cycle again.
        pCell->pChild[i] = create_child_cell(i, pCell);
        pCell->pChild[i]->charges.push_back(pCell->charges.at(0));
        pCell->pChild[i]->charges.push_back(pCell->charges.at(1));
        pCell = pCell->pChild[i];
      }
    }
  }
}

int CBarnesHut::child_cell_index(const Vector3D& d)
{
  int i = 0;
  if(d.x > 0) i++;
  if(d.y > 0) i += 2;
  if(d.z > 0) i += 4;
  return i;
}

OctoTreeCell* CBarnesHut::create_child_cell(int child_index, OctoTreeCell* pParentCell)
{
  Vector3D c(pParentCell->box_center);
  double half_size = 0.5 * pParentCell->size;  // size of the daughter's cell is twice as small as the parent's cell.
  switch(child_index)
  {
    case 0: c.x -= half_size; c.y -= half_size; c.z -= half_size; break;
    case 1: c.x += half_size; c.y -= half_size; c.z -= half_size; break;
    case 2: c.x -= half_size; c.y += half_size; c.z -= half_size; break;
    case 3: c.x += half_size; c.y += half_size; c.z -= half_size; break;

    case 4: c.x -= half_size; c.y -= half_size; c.z += half_size; break;
    case 5: c.x += half_size; c.y -= half_size; c.z += half_size; break;
    case 6: c.x -= half_size; c.y += half_size; c.z += half_size; break;
    case 7: c.x += half_size; c.y += half_size; c.z += half_size; break;
  }

  OctoTreeCell* pCell = new OctoTreeCell(c, half_size);
  pCell->ndep = pParentCell->ndep + 1;
  m_vCells.push_back(pCell);
  return pCell;
}

void CBarnesHut::create_main_cell(const Vector3D& cube_center, double cube_edge)
{
  m_pTree = new OctoTreeCell(cube_center, cube_edge);
  m_vCells.push_back(m_pTree);
}

void CBarnesHut::prepare(CalcThreadVector& vThreads, CNodesCollection& vNodes, UINT nIter)
{
  compute_moments(vThreads);

  CNode3D* pNode = NULL;
  size_t nNodesCount = vNodes.size();
  for(size_t i = 0; i < nNodesCount; i++)
  {
    pNode = vNodes.at(i);
    pNode->clmb = nIter == 1 ? coulomb_force(pNode->pos) : (double(nIter - 1) * pNode->clmb + coulomb_force(pNode->pos)) / (double)nIter;
// DEBUG: Visualization of the Coulomb potential.
    pNode->phi = coulomb_phi(pNode->pos);
// END DEBUG
  }

  m_bReady = true;
}

UINT CBarnesHut::build_moments(LPVOID pCalcThread)
{
  CalcThread* pThread = (CalcThread*)pCalcThread;
  CBarnesHut* pObj = (CBarnesHut*)pThread->m_pData;
  const CellsColl& vCells = pObj->get_cells();

  for(UINT nJob = pThread->get_first_job(); nJob <= pThread->get_last_job(); nJob++)
	{
    OctoTreeCell* pCell = vCells.at(nJob);
   
// Computation of the "center of charge" of the cell. All the moments will be computed with respect to this point.
    pCell->qSum = 0;
    Vector3D vPosChargeSum(0, 0, 0);
    UINT nChargesCount = pCell->charges.size();
    for(UINT i = 0; i < nChargesCount; i++)
    {
      ChargeData* pData = pCell->charges.at(i);
      vPosChargeSum += pData->pos * pData->charge;
      pCell->qSum += pData->charge;
    }

    if(pCell->qSum < Const_Almost_Zero)
      continue;
    
    pCell->charge_center = vPosChargeSum / pCell->qSum;

    if(pObj->m_bEnableQuadTerms)
    {
      double dx, dy, dz, dx2, dy2, dz2, dr2sum = 0;
      for(UINT i = 0; i < nChargesCount; i++)
      {
        ChargeData* pData = pCell->charges.at(i);

        dx = pData->pos.x - pCell->charge_center.x;
        dy = pData->pos.y - pCell->charge_center.y;
        dz = pData->pos.z - pCell->charge_center.z;

        dx2 = dx * dx;
        dy2 = dy * dy;
        dz2 = dz * dz;

        dr2sum += dx2 + dy2 + dz2;

        pCell->Qxx += dx2;
        pCell->Qyy += dy2;
        pCell->Qzz += dz2;
        pCell->Qxy += dx * dy;
        pCell->Qyz += dy * dz;
        pCell->Qxz += dz * dx;
      }

      dr2sum *= 0.3333333;
      pCell->Qxx -= dr2sum;
      pCell->Qyy -= dr2sum;
      pCell->Qzz -= dr2sum;

      pCell->Qxx *= 1.5;
      pCell->Qyy *= 1.5;
      pCell->Qzz *= 1.5;
      pCell->Qxy *= 1.5;
      pCell->Qyz *= 1.5;
		  pCell->Qxz *= 1.5;
    }

		pThread->done_job();
  }

  return 0;
}

void CBarnesHut::compute_moments(CalcThreadVector& vThreads)
{
	int nCellCount = m_vCells.size();
	if(nCellCount == 0)
		return;

  vThreads.distribute_jobs(0, nCellCount - 1, build_moments, (void*)this);
	vThreads.start_execution();
	vThreads.wait();
}

double CBarnesHut::coulomb_phi(const Vector3D& vPos)
{
  return coulomb_phi_cell(vPos, m_pTree); // potential in CGS.
}

double CBarnesHut::coulomb_phi_cell(const Vector3D& vPos, OctoTreeCell* pCell)
{
  Vector3D c = vPos - pCell->charge_center;
  double fSqrDist = c.sqlength();
  if(fSqrDist < Const_Almost_Zero)
    return 0;

  double fDist = sqrt(fSqrDist), phi = 0;
  UINT nChargesCount = pCell->charges.size();
  if(fDist > m_fDistCoeff * pCell->size || nChargesCount == 1 || pCell->ndep >= m_nMaxRecDepth)
  {
    phi = pCell->qSum / fDist;
    if(m_bEnableQuadTerms && (nChargesCount > 1))
    {
      double r5 = fSqrDist * fSqrDist * fDist;
      phi += quad_phi_cell(c, r5, pCell);
    }

    phi += correct_phi_for_symm(vPos, pCell);
	  return phi;
  }

  for(UINT i = 0; i < 8; i++)
  {
    OctoTreeCell* pChildCell = pCell->pChild[i];
    if(pChildCell != NULL)
      phi += coulomb_phi_cell(vPos, pChildCell);
  }

  return phi;
}

double CBarnesHut::correct_phi_for_symm(const Vector3D& vPos, OctoTreeCell* pCell)
{
  double fPhi = 0;
  switch(m_nSymmType)
  {
    case symXYonly:
    {
      fPhi += coulomb_phi_cell_symm(vPos, pCell, symXYonly);
      break;
    }
    case symXZonly:
    {
      fPhi += coulomb_phi_cell_symm(vPos, pCell, symXZonly);
      break;
    }
    case symBoth:
    {
      fPhi += coulomb_phi_cell_symm(vPos, pCell, symXYonly);
      fPhi += coulomb_phi_cell_symm(vPos, pCell, symXZonly);
      fPhi += coulomb_phi_cell_symm(vPos, pCell, symBoth);
      break;
    }
  }

  return fPhi;
}

double CBarnesHut::coulomb_phi_cell_symm(const Vector3D& vPos, OctoTreeCell* pCell, int nSymType)
{
  double fInv = -1;
  Vector3D vChargeCenter = pCell->charge_center;
  switch(nSymType)
  {
    case symXYonly: vChargeCenter.z *= fInv; break;
    case symXZonly: vChargeCenter.y *= fInv; break;
    case symBoth: vChargeCenter.y *= fInv; vChargeCenter.z *= fInv; break;
  }

  Vector3D c = vPos - vChargeCenter;
  double fSqrDist = c.sqlength();
  if(fSqrDist < Const_Almost_Zero)
    return 0;

  double fDist = sqrt(fSqrDist);
  double fPhi = pCell->qSum / fDist;
  if(m_bEnableQuadTerms && (pCell->charges.size() > 1))
  {
    double r5 = fSqrDist * fSqrDist * fDist;
    fPhi += quad_phi_cell(c, r5, pCell);
  }

  return fPhi;
}

double CBarnesHut::quad_phi_cell(const Vector3D c, double r5, OctoTreeCell* pCell)
{
  return pCell->qSum *
            (c.x * c.x * pCell->Qxx + c.y * c.y * pCell->Qyy + c.z * c.z * pCell->Qzz +
		         c.x * c.y * pCell->Qxy + c.y * c.z * pCell->Qyz + c.x * c.z * pCell->Qxz) / r5;
}

Vector3D CBarnesHut::coulomb_force(const Vector3D& vPos)
{
  return coulomb_force_cell(vPos, m_pTree); // E, field strength in CGS (to get force this must be multiplied by a probe charge).
}

Vector3D CBarnesHut::coulomb_force_cell(const Vector3D& vPos, OctoTreeCell* pCell)
{
  Vector3D c = vPos - pCell->charge_center;
  double fSqrDist = c.sqlength();
  if(fSqrDist < Const_Almost_Zero)
    return vNull;

  Vector3D f = vNull;
  double fDist = sqrt(fSqrDist);
  UINT nChargesCount = pCell->charges.size();
  if(fDist > m_fDistCoeff * pCell->size || nChargesCount == 1 || pCell->ndep >= m_nMaxRecDepth)
  {
    if((nChargesCount == 1) && (fDist < m_fCritRadius))
      return c * (pCell->qSum * m_fOne_ovr_R3);

    double rr2 = 1. / fSqrDist;
    double rr3 = rr2 / fDist;
    f = c * (pCell->qSum * rr3);  // repulsion force is always directed from cloud of ions, see definition of c above.

    if(m_bEnableQuadTerms && (nChargesCount > 1)) // quadripole terms.
      f += quad_field_cell(c, rr2, rr3, pCell);

    f += correct_field_for_symm(vPos, pCell);
    return f;
  }

  for(UINT i = 0; i < 8; i++)
  {
    OctoTreeCell* pChildCell = pCell->pChild[i];
    if(pChildCell != NULL)
      f += coulomb_force_cell(vPos, pChildCell);
  }

  return f;
}

Vector3D CBarnesHut::correct_field_for_symm(const Vector3D& vPos, OctoTreeCell* pCell)
{
  Vector3D f = vNull;
  switch(m_nSymmType)
  {
    case symXYonly:
    {
      f += coulomb_field_cell_symm(vPos, pCell, symXYonly);
      break;
    }
    case symXZonly:
    {
      f += coulomb_field_cell_symm(vPos, pCell, symXZonly);
      break;
    }
    case symBoth:
    {
      f += coulomb_field_cell_symm(vPos, pCell, symXYonly);
      f += coulomb_field_cell_symm(vPos, pCell, symXZonly);
      f += coulomb_field_cell_symm(vPos, pCell, symBoth);
      break;
    }
  }

  return f;
}

Vector3D CBarnesHut::coulomb_field_cell_symm(const Vector3D& vPos, OctoTreeCell* pCell, int nSymType)
{
  double fInv = -1;
  Vector3D vChargeCenter = pCell->charge_center;
  switch(nSymType)
  {
    case symXYonly: vChargeCenter.z *= fInv; break;
    case symXZonly: vChargeCenter.y *= fInv; break;
    case symBoth: vChargeCenter.y *= fInv; vChargeCenter.z *= fInv; break;
  }

  Vector3D c = vPos - vChargeCenter;
  double fSqrDist = c.sqlength();
  if(fSqrDist < Const_Almost_Zero)
    return vNull;

  Vector3D f = vNull;
  double fDist = sqrt(fSqrDist);
  UINT nChargesCount = pCell->charges.size();
  if((nChargesCount == 1) && (fDist < m_fCritRadius))
    return c * (pCell->qSum * m_fOne_ovr_R3);

  double rr2 = 1. / fSqrDist;
  double rr3 = rr2 / fDist;
  f = c * (pCell->qSum * rr3);  // repulsion force is always directed from cloud of ions, see definition of c above.

  if(m_bEnableQuadTerms && (nChargesCount > 1)) // quadripole terms.
    f += quad_field_cell(c, rr2, rr3, pCell);

  return f;
}

Vector3D CBarnesHut::quad_field_cell(const Vector3D c, double rr2, double rr3, OctoTreeCell* pCell)
{
  double rr5 = rr3 * rr2;
  double rr7 = rr5 * rr2;

  double xrr5 = c.x * rr5;
  double yrr5 = c.y * rr5;
  double zrr5 = c.z * rr5;

  double x2rr7 = 5 * c.x * c.x * rr7;
  double y2rr7 = 5 * c.y * c.y * rr7;
  double z2rr7 = 5 * c.z * c.z * rr7;
  double xyzrr7 = 5 * c.x * c.y * c.z * rr7;

  Vector3D f;

  f.x = pCell->Qxx * (2 * xrr5 - c.x * x2rr7) - c.x * (pCell->Qyy * y2rr7 + pCell->Qzz * z2rr7)
      + pCell->Qxy * (yrr5 - c.y * x2rr7) - pCell->Qyz * xyzrr7 + pCell->Qxz * (zrr5 - c.z * x2rr7);

  f.y = pCell->Qyy * (2 * yrr5 - c.y * y2rr7) - c.y * (pCell->Qxx * x2rr7 + pCell->Qzz * z2rr7)
      + pCell->Qxy * (xrr5 - c.x * y2rr7) + pCell->Qyz * (zrr5 - c.z * y2rr7) - pCell->Qxz * xyzrr7;

  f.z = pCell->Qzz * (2 * zrr5 - c.z * z2rr7) - c.z * (pCell->Qxx * x2rr7 + pCell->Qyy * y2rr7)
      - pCell->Qxy * xyzrr7 + pCell->Qyz * (yrr5 - c.y * z2rr7) + pCell->Qxz * (xrr5 - c.x * z2rr7);

  return pCell->qSum * f;
}

void CBarnesHut::set_crit_radius(double fR)
{
  m_fCritRadius = fR;
  m_fOne_ovr_R3 = fR * fR * fR;
  if(m_fOne_ovr_R3 > Const_Almost_Zero)
    m_fOne_ovr_R3 = 1. / m_fOne_ovr_R3;
  else
    m_fOne_ovr_R3 = 0.;
}

void CBarnesHut::scale_all_charges(double fCoeff)
{
  if(m_pTree == NULL)
    return;

  size_t nCount = m_pTree->charges.size();
  for(size_t i = 0; i < nCount; i++)
  {
    ChargeData* pData = m_pTree->charges.at(i);
    pData->charge *= fCoeff;
  }
}

};  // namespace EvaporatingParticle
