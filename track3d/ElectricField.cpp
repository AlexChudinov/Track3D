
#include "stdafx.h"

#include <algorithm>
#include <exception>

#include "ElectricField.h"
#include "ParticleTracking.h"
#include "BarnesHut.h"
#include "SelectedAreas.h"

#include "../field_solver/DCMeshAdapter.h"
#include "../field_solver/FiniteVolumesSolver.h"


namespace EvaporatingParticle
{

//---------------------------------------------------------------------------------------
// CPotentialBoundCond
//---------------------------------------------------------------------------------------
CPotentialBoundCond::CPotentialBoundCond(BoundaryMesh::BoundaryType type, int val)
  : m_nType(type), m_nFixedValType(val), m_bVisible(true)
{
  m_fStartX = 6.0;
  m_fStepX = 0.1;
  m_fEndX = 12.0;
// Linear potential function:
  m_fEndPhi = 0;
  m_nStepsCount = 10;
  m_fCenterFirstElectr = m_fStartX;
// Parabolic potential function:
  m_fStartPhi = 130;    // V.
  m_fFirstStepPhi = 2;  // V.
  m_fLastStepPhi = 4.7; // V.
// Merging options (if a Named Area is selected in addition to existing, manually selected regions).
  m_nMergeOpt = CSelectedAreas::optAdd;
  m_sLastMerged = _T("None");
}

const char* CPotentialBoundCond::get_bc_type_name(int nType)
{
  switch(nType)
  {
    case BoundaryMesh::FIXED_VAL: return _T("Fixed Value");
    case BoundaryMesh::ZERO_GRAD: return _T("Zero Normal Derivative");
  }

  return _T("None");
}

const char* CPotentialBoundCond::get_fixed_value_name(int nType)
{
  switch(nType)
  {
    case fvPlusUnity:     return _T("+1");
    case fvMinusUnity:    return _T("-1");
    case fvLinearStepsX:  return _T("Linear X Stepwise Potential");
    case fvLinearStepsY:  return _T("Linear Y Stepwise Potential");
    case fvQuadricStepsX: return _T("Quadric X Stepwise Potential");
    case fvQuadricStepsY: return _T("Quadric Y Stepwise Potential");
    case fvCoulomb:       return _T("Coulomb Potential");
  }

  return _T("0");
}

const char* CPotentialBoundCond::get_control_title(int nFixedValType, int nCtrlType)
{
  bool bX = nFixedValType == fvLinearStepsX || nFixedValType == fvQuadricStepsX;
  switch(nCtrlType)
  {
    case uitStartX:       return bX ? _T("Start X, mm") : _T("Start Y, mm");
    case uitStepX:        return bX ? _T("Period X, mm") : _T("Period Y, mm");
    case uitEndX:         return bX ? _T("End X, mm") : _T("End Y, mm");
    case uitStartPhi:     return _T("Start Potential, V");
    case uitEndPhi:       return _T("Dimensionless End Potential");
    case uitFirstStepPhi: return _T("Potential Diff First, V");
    case uitLastStepPhi:  return _T("Potential Diff Last, V");
    case uitStepsCount:  return _T("Electrodes Count");
    case uitCenterFirst:  return _T("First Electrode Center, mm");
  }

  return _T("");
}

const char* CPotentialBoundCond::get_hint(int nFixedValType, int nCtrlType)
{
  bool bX = nFixedValType == fvLinearStepsX || nFixedValType == fvQuadricStepsX;
  CString csXY = bX ? _T("X") : _T("Y");
  switch(nCtrlType)
  {
    case uitStartX:       return CString(_T("Specify the ")) + csXY + CString(_T("-coordinate (in mm) of the first electrode center in the step-wise potentials set."));
    case uitStepX:        return CString(_T("Specify the ")) + csXY + CString(_T("-period, i.e. the space between centers of the equidistant electrodes, mm."));
    case uitEndX:         return CString(_T("Specify the ")) + csXY + CString(_T("-coordinate (in mm) of the last electrode center in the step-wise potentials set."));
    case uitStartPhi:     return _T("Specify the start potential value in V. Important: Make sure the <Voltage Scale> edit control reads unity!");
    case uitEndPhi:       return _T("Dimensionless End Potential");
    case uitFirstStepPhi: return _T("Set the potential difference between the second and first electrodes in V. Important: Make sure the <Voltage Scale> edit control reads unity!");
    case uitLastStepPhi:  return _T("Set the potential difference between the last and last but one electrodes in V. Important: Make sure the <Voltage Scale> edit control reads unity!");
    case uitStepsCount:  return _T("Count of electrodes, on which linear step-like potential is to be set.");
    case uitCenterFirst:  return _T("The coordinate for the center of the first electrode, on which linear step-like potential is to be set, in mm");
  }

  return _T("");
}

double CPotentialBoundCond::linear_potential(double x)
{
  double dx = m_fEndX - m_fStartX;
  if(fabs(dx) < Const_Almost_Zero)
    return 0;

  return 1 + (m_fEndPhi - 1) * (x - m_fStartX) / dx;
}

void CPotentialBoundCond::save(CArchive& ar)
{
  UINT nVersion = 5;  // 5 - Merging option; 4 - nStepsCount and fCenterFirstElectr; 3 - linear and quadric stepwise potentials; 2 - hide/show regions; 1 - step-like potential is supported.
  ar << nVersion;

  ar << m_nType;
  ar << m_nFixedValType;

  CString cStr;
  size_t nCount = m_vRegNames.size();
  ar << nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    cStr = CString(m_vRegNames.at(i).c_str());
    ar << cStr;
  }

  ar << m_fStartX;
  ar << m_fStepX;
  ar << m_fEndX;

  ar << m_fStartPhi;
  ar << m_fEndPhi;
  ar << m_fFirstStepPhi;
  ar << m_fLastStepPhi;

  ar << m_bVisible;

  ar << m_fCenterFirstElectr;
  ar << m_nStepsCount;

  ar << m_nMergeOpt;
  ar << m_sLastMerged;
}

void CPotentialBoundCond::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_nType;
  ar >> m_nFixedValType;

  CString cStr;
  size_t nCount;
  ar >> nCount;
  m_vRegNames.reserve(nCount);

  for(size_t i = 0; i < nCount; i++)
  {
    ar >> cStr;
    std::string sRegName((const char*)cStr);
    m_vRegNames.push_back(sRegName);
  }

  if(nVersion >= 1)
  {
    ar >> m_fStartX;
    ar >> m_fStepX;
    ar >> m_fEndX;
  }

  if(nVersion >= 3)
  {
    ar >> m_fStartPhi;
    ar >> m_fEndPhi;
    ar >> m_fFirstStepPhi;
    ar >> m_fLastStepPhi;
  }

  if(nVersion >= 2)
    ar >> m_bVisible;

  if(nVersion >= 4)
  {
    ar >> m_fCenterFirstElectr;
    ar >> m_nStepsCount;
  }

  if(nVersion >= 5)
  {
    ar >> m_nMergeOpt;
    ar >> m_sLastMerged;
  }
}

//-------------------------------------------------------------------------------------------------
// CFieldDataCollection.
//-------------------------------------------------------------------------------------------------
CFieldDataColl::CFieldDataColl()
{
  m_nCurrField = -1;
}

CFieldDataColl::~CFieldDataColl()
{
  clear_fields();
}

void CFieldDataColl::clear_fields()
{
  for(size_t i = 0; i < size(); i++)
  {
    CElectricFieldData* pObj = at(i);
    delete pObj;
  }

  clear();
}

bool CFieldDataColl::calc_fields(bool bMirrorClmb)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  HWND hProgress, hJobName, hDlgWnd;
  pObj->get_handlers(hJobName, hProgress, hDlgWnd);

// bMirrorClmb == false:
//   The call is from CPropertiesWnd::OnDoTracking(); only the steady-state (not Mirror) fields must be calculated if. In fact, only those fields
//   will be re-calculated, which have their m_bNeedRecalc flags raised.
// bMirrorClmb == true:
//   The call is from CTracker::create_BH_object(); only the Mirror field must be re-calculated on every iteration.

  if(!bMirrorClmb)
    clear_fields_in_nodes();

  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = at(i);
    bool bDoCalc = bMirrorClmb ? pData->get_type() == CElectricFieldData::typeMirror : pData->get_type() != CElectricFieldData::typeMirror;
    if(!bDoCalc)
      continue;

    if(pObj->get_terminate_flag())
      return false;

    pData->set_handlers(hJobName, hProgress, hDlgWnd);
    if(!pData->calc_field())
    {
      pData->set_handlers(NULL, NULL, NULL);
      return false;
    }

    pData->set_handlers(NULL, NULL, NULL);
  }

  if(!bMirrorClmb)
    pObj->apply_perturbations();  // [MS] 20-01-2019, the perturbations are applied at every node before the tracking starts instead of being called during tracking.

  return true;
}

static const Vector3D scvNull(0, 0, 0);

void CFieldDataColl::clear_fields_in_nodes()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  pObj->set_job_name("Clearing field in nodes...");
  CNodesCollection& vNodes = pObj->get_nodes();

  CNode3D* pNode = NULL;
  size_t nCount = vNodes.size();
  for(size_t i = 0; i < nCount; i++)
  {
    if((i % 100 == 0) && !pObj->get_terminate_flag())
      pObj->set_progress(100 * i / nCount);

    pNode = vNodes.at(i);
    pNode->field = scvNull;
    pNode->clmb = scvNull;
    pNode->rf = scvNull;
  }
}

bool CFieldDataColl::sel_region_changed(CStringVector* pRegNames)
{
  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = at(i);
    size_t nBoundCondCount = pData->get_bc_count();
    for(size_t j = 0; j < nBoundCondCount; j++)
    {
      CPotentialBoundCond* pBC = pData->get_bc(j);
      if((DWORD_PTR)pRegNames == pBC->get_region_names_ptr())
      {
        pData->invalidate();
        return true;
      }
    }
  }

  return false;
}

bool CFieldDataColl::remove_bound_cond(CPotentialBoundCond* pBC)
{
  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();

  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = at(i);
    size_t nBoundCondCount = pData->get_bc_count();
    for(size_t j = 0; j < nBoundCondCount; j++)
    {
      if(pBC == pData->get_bc(j))
      {
// Set "true" for the visibility flag of all the regions belonging to this BC before deleting the BC.
        CStringVector* pRegNames = (CStringVector*)(pBC->get_region_names_ptr());
        pDrawObj->set_visibility_status(pRegNames, true);
        pData->remove_bc(j);
        return true;
      }
    }
  }

  return false;
}

void CFieldDataColl::update_visibility_status()
{
  CPotentialBoundCond* pBC = NULL;
  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();

  size_t nCount = size();
  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = at(i);
    size_t nBoundCondCount = pData->get_bc_count();
    for(size_t j = 0; j < nBoundCondCount; j++)
    {
      pBC = pData->get_bc(j);
      CStringVector* pRegNames = (CStringVector*)(pBC->get_region_names_ptr());
      pDrawObj->set_visibility_status(pRegNames, pBC->get_visibility_flag());
    }
  }
}

void CFieldDataColl::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  size_t nCount = size();
  ar << nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = at(i);
    pData->save(ar);
  }
}

void CFieldDataColl::load(CArchive& ar)
{
  clear_fields();

  UINT nVersion;
  ar >> nVersion;

  size_t nCount;
  ar >> nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    CElectricFieldData* pData = new CElectricFieldData();
    pData->load(ar);
    push_back(pData);
  }

  m_nCurrField = size() > 0 ? 0 : -1;
}

//---------------------------------------------------------------------------------------
// CElectricFieldData - an auxiliary class for simulating scalable electric fields.
//---------------------------------------------------------------------------------------
CElectricFieldData::CElectricFieldData(int nType)
  : m_nType(nType)
{
  set_default();
}

CElectricFieldData::~CElectricFieldData()
{
  clear_bc();
}

void CElectricFieldData::set_default()
{
  m_bEnable = true;
  m_bMultiThread = true;  // for tests only, it is never saved to the stream.
  m_bEnableVis = false;   // visualization flag.

  m_nCalcMethod = cmLaplacian3;   // the oldest and most tested method so far.

  m_fScale = 10 * SI_to_CGS_Voltage;  // 10 V in CGS.
  set_freq(1.0e+6);                   // 1 MHz by default.

  m_bNeedRecalc = true;
  m_nIterCount = 3000;
  m_fTol = 1e-5;

// An attempt to get analytic field in the flatapole. Alpha version.
  m_bAnalytField = false;
  m_fRadius = 0.21;   // cm, inscribed radius of the flatapole electrodes.
  m_fLowLimX = 14.7;  // cm, an analytic formula will be used if m_fLowLimX < x < m_fHighLimX;
  m_fHighLimX = 1e+5;

  m_sName = default_name();
}

bool CElectricFieldData::calc_field()
{
  if(!m_bEnable)
    return true;

  bool bOK = m_vField.size() == CParticleTrackingApp::Get()->GetTracker()->get_nodes().size();

  if(!m_bNeedRecalc && bOK)
    return get_result();

  switch(m_nCalcMethod)
  {
    case cmLaplacian3: return calc_lap3();
    case cmDirTessLap3: return calc_dirichlet_lap3();
    case cmFinVolJacobi: return calc_finite_vol_jacobi();
	case cmEigenLibSolver: return calc_eigen_lib_lap();
  }

  return false;
}

bool CElectricFieldData::calc_lap3()
{
  try
  {
    terminate(false);
    CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
    const CNodesCollection& vNodes = pObj->get_nodes();
    CMeshAdapter mesh(pObj->get_elems(), vNodes);
    mesh.progressBar()->set_handlers(m_hJobNameHandle, m_hProgressBarHandle, m_hDlgWndHandle);

    check_regions();
    if(!set_default_boundary_conditions(mesh) || !set_boundary_conditions(mesh))
      return false;

    size_t nNodesCount = vNodes.size();
    std::vector<double> field(nNodesCount, 0.0);

// Laplacian Solver:
    mesh.boundaryMesh()->applyBoundaryVals(field);
    CMeshAdapter::PScalFieldOp pOp = mesh.createOperator(CMeshAdapter::LaplacianSolver3);
    if(pOp == NULL)
      return false;

    set_job_name(job_name(jobCalcField));
    for(int i = 0; i < m_nIterCount; ++i)
    {
      if(get_terminate_flag())
      {
        BlockPoolInterface::cleanUpEveryPool();
        return false;
      }

      pOp->applyToField(field);
      set_progress(100 * i / m_nIterCount);
    }

// Gradient calculation:
    CMeshAdapter::PScalFieldOp pGrad = mesh.createOperator(CMeshAdapter::GradX);
    std::vector<double> dPhiDx(field);
    pGrad->applyToField(dPhiDx);

    pGrad = mesh.createOperator(CMeshAdapter::GradY);
    std::vector<double> dPhiDy(field);
    pGrad->applyToField(dPhiDy);

    pGrad = mesh.createOperator(CMeshAdapter::GradZ);
    std::vector<double> dPhiDz(field);
    pGrad->applyToField(dPhiDz);

    CNode3D* pNode = NULL;
    m_vPhi.resize(nNodesCount, 0);
    m_vField.resize(nNodesCount);
    for(size_t i = 0; i < nNodesCount; i++)
    {
      m_vPhi[i] = (float)field[i];
      m_vField[i] = Vector3F(-(float)dPhiDx[i], -(float)dPhiDy[i], -(float)dPhiDz[i]);

// An attempt to get analytic field in the flatapole. Alpha version.
      if(m_bAnalytField && (m_nType == typeFieldRF))
        apply_analytic_field(vNodes.at(i)->pos, m_vField[i], m_vPhi[i]);
    }

    m_bNeedRecalc = m_nType == typeMirror;  // false;
	}
  catch (const std::exception& ex)
  {
    AfxMessageBox(ex.what());
  }

  bool bRes = get_result();
  notify_scene();

  /*[AC 03.05.2017] Memory clean up */
  BlockPoolInterface::cleanUpEveryPool();
  /*[/AC]*/

  return bRes;
}

bool CElectricFieldData::calc_dirichlet_lap3()
{
  CDirichletTesselation* pTessObj = CParticleTrackingApp::Get()->GetDirichletTess();
  try
  {
    terminate(false);
    CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
    const CNodesCollection& vNodes = pObj->get_nodes();
    DCMeshAdapter mesh(pObj->get_elems(), vNodes, *pTessObj);
    mesh.progressBar()->set_handlers(m_hJobNameHandle, m_hProgressBarHandle, m_hDlgWndHandle);

    check_regions();
    if(!set_default_boundary_conditions(mesh) || !set_boundary_conditions(mesh))
      return false;

    size_t nNodesCount = vNodes.size();
    std::vector<double> field(nNodesCount, 0.0);

// Laplacian Solver:
    mesh.boundaryMesh()->applyBoundaryVals(field);
    DCMeshAdapter::PScalFieldOp pOp = mesh.createOperator(DCMeshAdapter::LaplacianSolver, NULL);
    if(pOp == NULL)
      return false;

    set_job_name(job_name(jobCalcField));
    CFloatArray vMaxRelErr;
    vMaxRelErr.assign(m_nIterCount, 0);
    double fErr;
    for(int i = 0; i < m_nIterCount; ++i)
    {
      if(get_terminate_flag())
      {
        BlockPoolInterface::cleanUpEveryPool();
        return false;
      }

      fErr = 0;
      pOp->applyToField(field, &fErr);
      vMaxRelErr[i] = (float)fErr;
      set_progress(100 * i / m_nIterCount);
      if(fErr < m_fTol)
        break;
    }

// Gradient calculation:
    DCMeshAdapter::PScalFieldOp pGrad = mesh.createOperator(DCMeshAdapter::GradX, NULL);
    std::vector<double> dPhiDx(field);
    pGrad->applyToField(dPhiDx);

    pGrad = mesh.createOperator(DCMeshAdapter::GradY, NULL);
    std::vector<double> dPhiDy(field);
    pGrad->applyToField(dPhiDy);

    pGrad = mesh.createOperator(DCMeshAdapter::GradZ, NULL);
    std::vector<double> dPhiDz(field);
    pGrad->applyToField(dPhiDz);

    CNode3D* pNode = NULL;
    m_vPhi.resize(nNodesCount, 0);
    m_vField.resize(nNodesCount);
    for(size_t i = 0; i < nNodesCount; i++)
    {
      m_vPhi[i] = (float)field[i];
      m_vField[i] = Vector3F(-(float)dPhiDx[i], -(float)dPhiDy[i], -(float)dPhiDz[i]);

// An attempt to get analytic field in the flatapole. Alpha version.
      if(m_bAnalytField && (m_nType == typeFieldRF))
        apply_analytic_field(vNodes.at(i)->pos, m_vField[i], m_vPhi[i]);
    }

    m_bNeedRecalc = m_nType == typeMirror;

    output_convergence_history(vMaxRelErr);
	}
  catch (const std::exception& ex)
  {
    AfxMessageBox(ex.what());
  }

  bool bRes = get_result();
  notify_scene();

  /*[AC 03.05.2017] Memory clean up */
  BlockPoolInterface::cleanUpEveryPool();
  /*[/AC]*/

  return bRes;
}

bool CElectricFieldData::calc_finite_vol_jacobi()
{
  CAnsysMesh* pMesh = (CAnsysMesh*)CParticleTrackingApp::Get()->GetTracker();
  CDirichletTesselation* pTess = CParticleTrackingApp::Get()->GetDirichletTess();

  try
  {
    terminate(false);
    CFiniteVolumesSolver solver(pMesh, pTess);

    solver.set_handlers(m_hJobNameHandle, m_hProgressBarHandle, m_hDlgWndHandle);

    check_regions();
    set_boundary_conditions(solver);

    CFloatArray vPhi;
    float fTol = (float)m_fTol;
    set_job_name(job_name(jobCalcField));
    CSolutionInfo info = solver.solve(fTol, m_nIterCount, vPhi, m_bMultiThread);

    solver.set_handlers(NULL, NULL, NULL);
    if(!info.bSuccess)
      return false;

    CNodesCollection& vNodes = pMesh->get_nodes();
    size_t nNodeCount = vNodes.size();
    m_vPhi.resize(nNodeCount);
    m_vField.resize(nNodeCount);
    Vector3D vGrad;
    for(size_t k = 0; k < nNodeCount; k++)
    {
      vGrad = pTess->get_grad(k, vPhi);
      m_vField[k] = Vector3F(-(float)vGrad.x, -(float)vGrad.y, -(float)vGrad.z);
      m_vPhi[k] = vPhi.at(k);

// An attempt to get analytic field in the flatapole. Alpha version.
      if(m_bAnalytField && (m_nType == typeFieldRF))
        apply_analytic_field(vNodes.at(k)->pos, m_vField[k], m_vPhi[k]);
    }

    output_convergence_history(info.vMaxErrHist);
  }
  catch(const std::exception& ex)
  {
    AfxMessageBox(ex.what());
  }

  m_bNeedRecalc = m_nType == typeMirror;
  notify_scene();

  return get_result();
}

bool CElectricFieldData::calc_eigen_lib_lap()
{
	try
	{
		terminate(false);
		EvaporatingParticle::CTracker * pObj = CParticleTrackingApp::Get()->GetTracker();
		CMeshAdapter mesh(pObj->get_elems(), pObj->get_nodes());
		mesh.progressBar()->set_handlers(m_hJobNameHandle, m_hProgressBarHandle, m_hDlgWndHandle);

		check_regions();
		if (!set_default_boundary_conditions(mesh) || !set_boundary_conditions(mesh))
			return false;

		size_t nNodesCount = pObj->get_nodes().size();
		std::vector<double> field(nNodesCount, 0.0);

		// Laplacian Solver:
		mesh.boundaryMesh()->applyBoundaryVals(field);
		DCMeshAdapter::PScalFieldOp pOp = mesh.createOperator(CMeshAdapter::EigenLibLaplacian, NULL);
		if (pOp == NULL)
			return false;

		set_job_name(job_name(jobCalcField));
		double fErr;
		set_progress(50);
		pOp->applyToField(field, &fErr);
		set_progress(100);

		// Gradient calculation:
		CMeshAdapter::PScalFieldOp pGrad = mesh.createOperator(CMeshAdapter::GradX, NULL);
		std::vector<double> dPhiDx(field);
		pGrad->applyToField(dPhiDx);

		pGrad = mesh.createOperator(CMeshAdapter::GradY, NULL);
		std::vector<double> dPhiDy(field);
		pGrad->applyToField(dPhiDy);

		pGrad = mesh.createOperator(CMeshAdapter::GradZ, NULL);
		std::vector<double> dPhiDz(field);
		pGrad->applyToField(dPhiDz);

		CNode3D* pNode = NULL;
		m_vPhi.resize(nNodesCount, 0);
		m_vField.resize(nNodesCount);
		for (size_t i = 0; i < nNodesCount; i++)
		{
			m_vPhi[i] = (float)field[i];
			m_vField[i] = Vector3F(-(float)dPhiDx[i], -(float)dPhiDy[i], -(float)dPhiDz[i]);

		}

		m_bNeedRecalc = m_nType == typeMirror;

	}
	catch (const std::exception& ex)
	{
		AfxMessageBox(ex.what());
	}

	notify_scene();
	return get_result();
}

void CElectricFieldData::output_convergence_history(const std::vector<float>& vMaxRelErr)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  COutputEngine& out_engine = pObj->get_output_engine();
  std::string cPath = out_engine.get_full_path(pObj->get_filename());
  std::string cName("Convergence_Hist");
  std::string cExt(".csv");
  std::string cFileName = cPath + cName + cExt;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if((nErr == 0) && (pStream != NULL))
  {
    UINT nIterCount = vMaxRelErr.size();
    for(UINT i = 0; i < nIterCount; i++)
      fprintf(pStream, "%d,   %f\n", i, vMaxRelErr.at(i));

    fclose(pStream);
  }
}

bool CElectricFieldData::set_boundary_conditions(CFiniteVolumesSolver& solver)
{
  set_default_boundary_conditions(solver);

  set_job_name(job_name(jobUserCond));
  set_progress(0);

  float fVal;
  CRegion* pReg = NULL;
  CPotentialBoundCond* pBC = NULL;
  size_t nCount = m_vBoundCond.size(), nRegCount;
  for(size_t i = 0; i < nCount; i++)
  {
    pBC = m_vBoundCond.at(i);
    nRegCount = pBC->get_region_names().size();
    int nPart = 100 * i / nCount;
    for(size_t j = 0; j < nRegCount; j++)
    {
      set_progress(nPart + int(0.5 + 100. * (j + 1) / (nRegCount * nCount)));
      if(get_terminate_flag())
        return false;

      pReg = CAnsysMesh::get_region(pBC->get_region_names().at(j));
      if(pReg == NULL)  // if the geometry file has been changed in the project, some region names can be not relevant.
        continue;

      CIndexVector vRegNodeInd = get_reg_nodes(pReg);   // global indices of nodes, belonging to region pReg.

      if(pBC->get_bc_type() == BoundaryMesh::ZERO_GRAD)
      {
        solver.set_boundary_conditions(vRegNodeInd);    // boundary conditions of the 2-nd type.
      }
      else if(pBC->get_bc_type() == BoundaryMesh::FIXED_VAL)
      {
        switch(pBC->get_fixed_val_type())
        {
          case CPotentialBoundCond::fvPlusUnity:
          {
            solver.set_boundary_conditions(vRegNodeInd, 1);
            break;
          }
          case CPotentialBoundCond::fvMinusUnity:
          {
            solver.set_boundary_conditions(vRegNodeInd, -1);
            break;
          }
          case CPotentialBoundCond::fvLinearStepsX:
          case CPotentialBoundCond::fvLinearStepsY:
          {
            fVal = (float)(linear_step_potential(pBC, pReg->vFaces.at(0)->p0->pos));
            solver.set_boundary_conditions(vRegNodeInd, fVal);
            break;
          }
          case CPotentialBoundCond::fvQuadricStepsX:
          case CPotentialBoundCond::fvQuadricStepsY:
          {
            fVal = (float)(quadric_step_potential(pBC, pReg->vFaces.at(0)->p0->pos));
            solver.set_boundary_conditions(vRegNodeInd, fVal);
            break;
          }
          case CPotentialBoundCond::fvCoulomb:
          {
            CFloatArray vClmbPhi;
            if(!coulomb_potential(vRegNodeInd, vClmbPhi))
              break;

            solver.set_boundary_conditions(vRegNodeInd, vClmbPhi);
            break;
          }
        }
      }
    }
  }

  return true;
}

bool CElectricFieldData::set_default_boundary_conditions(CFiniteVolumesSolver& solver)
{
  set_job_name(job_name(jobDfltCond));
  set_progress(0);

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CRegionsCollection& vRegions = pObj->get_regions(false);
  size_t nRegCount = vRegions.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    if(get_terminate_flag())
      return false;

    CRegion* pReg = vRegions.at(i);
    if(is_selected(pReg) || pReg->bCrossSection)
      continue;

    CIndexVector vRegNodeInd = get_reg_nodes(pReg);
    solver.set_boundary_conditions(vRegNodeInd, 0);

    set_progress(100 * i / nRegCount);
  }

  return true;
}

bool CElectricFieldData::get_result() const
{
  if(m_nType != typeFieldDC)
    return true;

  CNodesCollection& vNodes = CParticleTrackingApp::Get()->GetTracker()->get_nodes();
  size_t nNodeCount = vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
    vNodes.at(i)->field += Vector3D(m_vField[i].x, m_vField[i].y, m_vField[i].z) * m_fScale;

  return true;
}

void CElectricFieldData::apply_analytic_field(const Vector3D& vPos, Vector3F& vField, float& fPhi)
{
  if(vPos.x < m_fLowLimX || vPos.x > m_fHighLimX)
    return;

  double fR = sqrt(vPos.y * vPos.y + vPos.z * vPos.z);
  if(fR > m_fRadius)
    return;

  double fCoeff = 2. / (m_fRadius * m_fRadius);
  Vector3F vAnalytRF(0, -fCoeff * vPos.y, fCoeff * vPos.z);
  double fAnalytePhi = 0.5 * (vPos.y * vPos.y - vPos.z * vPos.z) * fCoeff;

  const double scfTransWidth = 0.2;
  double fLowLimMax = m_fLowLimX + scfTransWidth; // the field changes gradually from numeric to analytic; transition interval is 2 mm long.
  if((vPos.x >= m_fLowLimX) && (vPos.x <= fLowLimMax))
  {
    float fKsi = float((vPos.x - m_fLowLimX) / scfTransWidth);
    vField = vAnalytRF * fKsi + vField * (1 - fKsi);
    fPhi = float(fAnalytePhi) * fKsi + fPhi * (1 - fKsi);
    return;
  }

  vField = vAnalytRF;
  fPhi = fAnalytePhi;
}

bool CElectricFieldData::set_boundary_conditions(CMeshAdapter& mesh)
{
  set_job_name(job_name(jobUserCond));
  set_progress(0);

  CRegion* pReg = NULL;
  size_t i, j, nBoundCondCount = m_vBoundCond.size(), nRegCount, nSumRegCount = 0;
  for(i = 0; i < nBoundCondCount; i++)
    nSumRegCount += m_vBoundCond.at(i)->get_region_names().size();

  size_t k = 0;
  for(i = 0; i < nBoundCondCount; i++)
  {
    CPotentialBoundCond* pBC = m_vBoundCond.at(i);
    nRegCount = pBC->get_region_names().size();
    for(j = 0; j < nRegCount; j++)
    {
      if(get_terminate_flag())
        return false;

      k++;
      const std::string& sName = pBC->get_region_names().at(j);
      pReg = CAnsysMesh::get_region(sName);
      if(pReg == NULL)
        continue;

      set_boundary_values(mesh, pReg, pBC);

      set_progress(100 * k / nSumRegCount);
    }
  }

  return true;
}

bool CElectricFieldData::set_default_boundary_conditions(CMeshAdapter& mesh)
{
  set_job_name(job_name(jobDfltCond));
  set_progress(0);

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CRegionsCollection& vRegions = pObj->get_regions(false);
  size_t nRegCount = vRegions.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    if(get_terminate_flag())
      return false;

    CRegion* pReg = vRegions.at(i);
    if(is_selected(pReg) || pReg->bCrossSection)
      continue;

    set_boundary_values(mesh, pReg);

    set_progress(100 * i / nRegCount);
  }

  return true;
}

void CElectricFieldData::set_boundary_values(CMeshAdapter& mesh, CRegion* pReg, CPotentialBoundCond* pBC)
{
  CFace* pFace = NULL;
  std::vector<UINT> vNodeIds;
  std::vector<Vector3D> vNorms;
  size_t nBndFacesCount = pReg->vFaces.size();
  if(nBndFacesCount == 0)
    return;

  for(size_t k = 0; k < nBndFacesCount; k++)
  {
    pFace = pReg->vFaces.at(k);
    vNodeIds.push_back(pFace->p0->nInd);
    vNorms.push_back(calc_norm(pFace->p0));

    vNodeIds.push_back(pFace->p1->nInd);
    vNorms.push_back(calc_norm(pFace->p1));

    vNodeIds.push_back(pFace->p2->nInd);
    vNorms.push_back(calc_norm(pFace->p2));
  }

  const std::string& sName = pReg->sName;
  mesh.boundaryMesh()->addBoundary(sName, vNodeIds, vNorms);

  BoundaryMesh::BoundaryType nType = pBC != NULL ? (BoundaryMesh::BoundaryType)(pBC->get_bc_type()) : BoundaryMesh::FIXED_VAL;
  mesh.boundaryMesh()->boundaryType(sName, nType);

  if(nType != BoundaryMesh::FIXED_VAL)
    return;

// There are two cases: 
// 1) Ordinary field specified by a single constant potential value at all nodes of a region (metallic electrode);
// 2) Mirror Coulomb field, boundary potential values for which can be different at different nodes of a region.
  double fVal = 0;
  if((pBC != NULL) && (nType == BoundaryMesh::FIXED_VAL))
  {
    if(pBC->get_fixed_val_type() == CPotentialBoundCond::fvCoulomb)  // mirror Coulomb field.
    {
      CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
      CBarnesHut* pBHObj = pObj->get_BH_object();

      const BoundaryMesh::SetLabels& vNodeLabels = mesh.boundaryMesh()->boundaryLabels(sName);
      size_t nRegNodeCount = vNodeLabels.size();
      std::vector<double> vBoundPhi;
      BoundaryMesh::SetLabels::const_iterator citer;
      UINT nNodeId = 0;
      Vector3D vPos;
      for(citer = vNodeLabels.cbegin(); citer != vNodeLabels.cend(); citer++)
      {
        nNodeId = *citer;
        vPos = pObj->get_nodes().at(nNodeId)->pos;
        fVal = -pBHObj->coulomb_phi(vPos);
        vBoundPhi.push_back(fVal);
      }

      mesh.boundaryMesh()->boundaryVals(sName, vBoundPhi);
      return;
    }

    switch(pBC->get_fixed_val_type())
    {
      case CPotentialBoundCond::fvPlusUnity: fVal = 1.0; break;
      case CPotentialBoundCond::fvMinusUnity: fVal = -1.0; break;
      case CPotentialBoundCond::fvLinearStepsX:
      case CPotentialBoundCond::fvLinearStepsY: fVal = linear_step_potential(pBC, pReg->vFaces.at(0)->p0->pos); break;
      case CPotentialBoundCond::fvQuadricStepsX:
      case CPotentialBoundCond::fvQuadricStepsY: fVal = quadric_step_potential(pBC, pReg->vFaces.at(0)->p0->pos); break;
    }
  }
    
  mesh.boundaryMesh()->boundaryVals(
	  sName, 
	  std::vector<double>(mesh.boundaryMesh()->patchSize(sName), fVal));
}

double CElectricFieldData::linear_step_potential(CPotentialBoundCond* pBC, const Vector3D& vPos) const
{
  double x = pBC->get_fixed_val_type() == CPotentialBoundCond::fvLinearStepsX ? vPos.x : vPos.y;

  bool bCorrect = pBC->get_start_coord() < pBC->get_end_coord();

  double x0 = bCorrect ? pBC->get_start_coord() : pBC->get_end_coord();
  double x1 = bCorrect ? pBC->get_end_coord() : pBC->get_start_coord();
  if(fabs(x1 - x0) < Const_Almost_Zero)
    return 0;

  double phi0 = bCorrect ? 1.0 : pBC->get_end_phi();
  double phi1 = bCorrect ? pBC->get_end_phi() : 1.0;
  double grad = (phi1 - phi0) / (x1 - x0);
  double dx = fabs(pBC->get_step_coord());
  if(dx < Const_Almost_Zero)
    return 0;

  if(pBC->get_steps_count() == 0)
    return 0;

  double xs = pBC->get_center_first_electr();
  double xe = xs + (pBC->get_steps_count() - 1) * dx;

  if(x < xs)
    return pBC->linear_potential(xs);
  else if(x < xe)
    return pBC->linear_potential(xs) + grad * floor((x - xs) / dx + 0.5)* dx;

  return pBC->linear_potential(xe);
}

double CElectricFieldData::quadric_step_potential(CPotentialBoundCond* pBC, const Vector3D& vPos) const
{
  double x = pBC->get_fixed_val_type() == CPotentialBoundCond::fvQuadricStepsX ? vPos.x : vPos.y;

  bool bCorrect = pBC->get_start_coord() < pBC->get_end_coord();
  if(!bCorrect)
    return 0;

  double x0 = pBC->get_start_coord();
  double x1 = pBC->get_end_coord();
  if(fabs(x1 - x0) < Const_Almost_Zero)
    return 0;

  double dx = fabs(pBC->get_step_coord());
  if(dx < Const_Almost_Zero)
    return 0;

  double hdx = 0.5 * dx;
  if(x < x0 - hdx || x > x1 + hdx)
    return 0;

  double dPhi0 = pBC->get_first_dphi() / dx; // V/cm.
  double dPhi1 = pBC->get_last_dphi() / dx;

  double A = 0.5 * (dPhi1 - dPhi0) / (x1 - x0);
  double B = (x1 * dPhi0 - x0 * dPhi1) / (x1 - x0);
  double C = pBC->get_start_phi() - A * x0 * x0 - B * x0;

  double xs = x0 + floor((x - x0) / dx + 0.5) * dx;
  return A * xs * xs + B * xs + C;
}

bool CElectricFieldData::coulomb_potential(const CIndexVector& vNodeIds, std::vector<float>& vPhi) const
{
  size_t nLocCount = vNodeIds.size();
  if(nLocCount == 0)
    return false;

  Vector3D vPos;
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CNodesCollection& vNodes = pObj->get_nodes();
  CBarnesHut* pBHObj = pObj->get_BH_object();
  vPhi.assign(nLocCount, 0);

  for(size_t i = 0; i < nLocCount; i++)
  {
    vPos = vNodes.at(vNodeIds.at(i))->pos;
    vPhi[i] = -(float)(pBHObj->coulomb_phi(vPos));
  }

  return true;
}

Vector3D CElectricFieldData::calc_norm(CNode3D* pNode) const
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  const CRegionsCollection& vRegs = pObj->get_regions(false);
  size_t nRegCount = vRegs.size();

  UINT nReg, nFace;
  CFace* pFace = NULL;
  Vector3D vNorm(0, 0, 0);
  size_t nNbrCount = pNode->vNbrFaces.size();
  for(size_t i = 0; i < nNbrCount; i++)
  {
    nReg = pNode->vNbrFaces.at(i).nReg;
    if(nReg >= nRegCount)
      continue;

    nFace = pNode->vNbrFaces.at(i).nFace;
    if(nFace >= vRegs.at(nReg)->vFaces.size())
      continue;

    pFace = vRegs.at(nReg)->vFaces.at(nFace);
    vNorm += pFace->norm;
  }

  vNorm.normalize();
  return -vNorm;
}

bool CElectricFieldData::is_selected(CRegion* pReg) const
{
  size_t nBoundCondCount = m_vBoundCond.size();
  for(size_t i = 0; i < nBoundCondCount; i++)
  {
    CPotentialBoundCond* pBC = m_vBoundCond.at(i);
    size_t nRegCount = pBC->get_region_names().size();
    for(size_t j = 0; j < nRegCount; j++)
    {
      const std::string& sName = pBC->get_region_names().at(j);
      if(pReg->sName == sName)
        return true;
    }
  }

  return false;
}

CIndexVector CElectricFieldData::get_reg_nodes(CRegion* pReg) const
{
  CIndexVector vInd;
  CFace* pFace = NULL;
  size_t nFaceCount = pReg->vFaces.size();
  for(size_t i = 0; i < nFaceCount; i++)
  {
    pFace = pReg->vFaces.at(i);
    if(std::find(vInd.begin(), vInd.end(), pFace->p0->nInd) == vInd.end())
      vInd.push_back(pFace->p0->nInd);
    if(std::find(vInd.begin(), vInd.end(), pFace->p1->nInd) == vInd.end())
      vInd.push_back(pFace->p1->nInd);
    if(std::find(vInd.begin(), vInd.end(), pFace->p2->nInd) == vInd.end())
      vInd.push_back(pFace->p2->nInd);
  }

  return vInd;
}

const char* CElectricFieldData::get_field_type_name(int nType)
{
  switch(nType)
  {
    case typeFieldDC: return _T("DC Electric Field");
    case typeFieldRF: return _T("Radio-Frequency Field");
    case typeMirror:  return _T("Mirror Coulomb Field");
  }

  return _T("None");
}

const char* CElectricFieldData::get_calc_method_name(int nCalcMethod)
{
  switch(nCalcMethod)
  {
    case cmLaplacian3: return _T("Laplacian Solver #3");
    case cmDirTessLap3: return _T("Advanced Laplacian Solver");
    case cmFinVolJacobi:  return _T("Finite Volume (Jacobi)");
	//Eigen lib solver
	case cmEigenLibSolver: return _T("Eigen lib laplacian");
  }

  return _T("None");
}

void CElectricFieldData::clear_bc()
{
  size_t nCount = m_vBoundCond.size();
  for(size_t i = 0; i < nCount; i++)
  {
    CPotentialBoundCond* pObj = m_vBoundCond.at(i);
    delete pObj;
  }

  m_vBoundCond.clear();
}

void CElectricFieldData::add_bc()
{
  char buff[4];
  size_t nId = m_vBoundCond.size();
  CString cName = CString("Boundary Conditions # ") + CString(itoa(nId + 1, buff, 10));

  CPotentialBoundCond* pBC = new CPotentialBoundCond();
  pBC->set_name(cName);
  m_vBoundCond.push_back(pBC);
}

void CElectricFieldData::remove_bc(size_t nId)
{
  if(nId < m_vBoundCond.size())
  {
    CPotentialBoundCond* pBC = m_vBoundCond.at(nId);
    m_vBoundCond.erase(m_vBoundCond.begin() + nId);
    delete pBC;
  }
}

void CElectricFieldData::notify_scene()
{
  CCrossSectColl* pColl = CParticleTrackingApp::Get()->GetPlanes();
  for(size_t i = 0; i < pColl->size(); i++)
    pColl->at(i)->invalidate();
}

CString CElectricFieldData::default_name()
{
  CString sName;
  char buff[4];
  size_t nFieldsCount = CParticleTrackingApp::Get()->GetFields()->size();
  bool bNameExists = true;
  int nNmbr = 0;
  while(bNameExists)
  {
    nNmbr++;
    sName = CString("Electric Field # ") + CString(itoa(nNmbr, buff, 10));
    bool bFound = false;
    for(size_t i = 0; i < nFieldsCount; i++)
    {
      bFound = CParticleTrackingApp::Get()->GetFields()->at(i)->get_field_name() == sName;
      if(bFound)
        break;
    }

    bNameExists = bFound;
  }

  return sName;
}

const char* CElectricFieldData::job_name(int nJobType) const
{
  CString sJobName = get_field_name();
  CString sSuffix = m_nType == typeMirror ? CString(_T(" (Mirror): ")) : CString(_T(": "));
  sJobName += sSuffix;
  switch(nJobType)
  {
    case jobDfltCond:  sJobName += CString(_T("Default boundary conditions...")); break;
    case jobUserCond:  sJobName += CString(_T("User-defined boundary conditions...")); break;
    case jobCalcField: sJobName += CString(_T("Solving...")); break;
  }

  return (const char*)sJobName;
}

void CElectricFieldData::check_regions()
{
  CRegion* pReg = NULL;
  CPotentialBoundCond* pBC = NULL;
  size_t nBndCondCount = m_vBoundCond.size(), nRegCount;
  for(size_t i = 0; i < nBndCondCount; i++)
  {
    pBC = m_vBoundCond.at(i);
    nRegCount = pBC->get_region_names().size();
    for(int j = nRegCount - 1; j >= 0; j--)
    {
      pReg = CAnsysMesh::get_region(pBC->get_region_names().at(j));
      if(pReg == NULL)
      {
        pBC->get_region_names().erase(pBC->get_region_names().begin() + j);
      }
    }
  }
}

double CElectricFieldData::get_stab_param(double fAmpl, double fFreq, double fInscR)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  double fIonMass = pObj->get_ion_mass();
  double fIonCharge = pObj->get_particle_charge();
  double fOmega = Const_2PI * fFreq;

  return 4 * fIonCharge * fAmpl / (fIonMass * fOmega * fOmega * fInscR * fInscR);
}

void CElectricFieldData::save(CArchive& ar)
{
  UINT nVersion = 7;  // 7 - m_vClmbPhi; 6 -  m_bEnableVis; 5 - m_vPhi; 4 - Tolerance; 3 - Calculation method; 2 - Analytic RF field kitchen, alpha version; 1 - m_bEnable and calculated
  ar << nVersion;

  ar << m_bEnable;
  ar << m_nType;
  ar << m_fScale;
  ar << m_fOmega;
  ar << m_nIterCount;
  ar << m_nCalcMethod;

  size_t nCount = m_vBoundCond.size();
  ar << nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    CPotentialBoundCond* pObj = m_vBoundCond.at(i);
    pObj->save(ar);
  }

  size_t nNodesCount = m_bNeedRecalc ? 0 : m_vField.size();
  ar << nNodesCount;
  for(size_t j = 0; j < nNodesCount; j++)
  {
    ar << m_vField[j].x;
    ar << m_vField[j].y;
    ar << m_vField[j].z;
    ar << m_vPhi[j];    // since version 5.
  }

  ar << m_bAnalytField;
  ar << m_fRadius;      // inscribed radius of the flatapole electrodes.
  ar << m_fLowLimX;     // an analytic formula will be used if m_fLowLimX < x < m_fHighLimX;
  ar << m_fHighLimX;
  ar << m_sName;

  ar << m_fTol;
  ar << m_bEnableVis;

  size_t nClmbPhiSize = m_vClmbPhi.size();  // it must be either nNodesCount or 0.
  ar << nClmbPhiSize;
  for(size_t k = 0; k < nClmbPhiSize; k++)
    ar << m_vClmbPhi[k];
}

void CElectricFieldData::load(CArchive& ar)
{
  clear_bc();

  UINT nVersion;
  ar >> nVersion;

  if(nVersion >= 1)
    ar >> m_bEnable;

  ar >> m_nType;
  ar >> m_fScale;
  ar >> m_fOmega;
  ar >> m_nIterCount;

  if(nVersion >= 3)
    ar >> m_nCalcMethod;

  size_t nCount;
  ar >> nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    add_bc();
    CPotentialBoundCond* pObj = m_vBoundCond.at(m_vBoundCond.size() - 1);
    pObj->load(ar);
  }

  if(nVersion >= 1)
  {
    size_t nNodesCount;
    ar >> nNodesCount;
    if(nNodesCount > 0)
    {
      m_vField.resize(nNodesCount, Vector3F(0, 0, 0));
      m_vPhi.resize(nNodesCount, 0.0f);
    }

    for(size_t j = 0; j < nNodesCount; j++)
    {
      ar >> m_vField[j].x;
      ar >> m_vField[j].y;
      ar >> m_vField[j].z;
      if(nVersion >= 5)
        ar >> m_vPhi[j];
    }

    if(nNodesCount > 0)
      m_bNeedRecalc = false;
  }

  if(nVersion >= 2)
  {
    ar >> m_bAnalytField;
    ar >> m_fRadius;      // inscribed radius of the flatapole electrodes.
    ar >> m_fLowLimX;     // an analytic formula will be used if m_fLowLimX < x < m_fHighLimX;
    ar >> m_fHighLimX;
    ar >> m_sName;
  }

  if(nVersion >= 4)
    ar >> m_fTol;

  if(nVersion >= 6)
    ar >> m_bEnableVis;

  if(nVersion >= 7)
  {
    size_t nClmbPhiSize; // it must be either nNodesCount or 0.
    ar >> nClmbPhiSize;
    if(nClmbPhiSize > 0)
    {
      m_vClmbPhi.resize(nClmbPhiSize, 0.0f);
      for(size_t k = 0; k < nClmbPhiSize; k++)
        ar >> m_vClmbPhi[k];
    }
  }
}

};  // namespace EvaporatingParticle.