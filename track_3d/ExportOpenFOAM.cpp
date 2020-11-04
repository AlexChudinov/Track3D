
#include "stdafx.h"


#include <stdio.h>
#include "float.h"

#include "ExportOpenFOAM.h"
#include "ParticleTracking.h"
#include "constant.hpp"

#include <algorithm>


namespace EvaporatingParticle
{

//---------------------------------------------------------------------------------------
// CBoundaryConditions
//---------------------------------------------------------------------------------------
void CBoundaryConditions::save(CArchive& ar)
{
  const UINT nVersion = 0;
  ar << nVersion;

  ar << fPress;
  ar << fTemp;
  ar << nType;

  CString cStr;
  size_t nCount = vRegNames.size();
  ar << nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    cStr = CString(vRegNames.at(i).c_str());
    ar << cStr;
  }
}

void CBoundaryConditions::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> fPress;
  ar >> fTemp;
  ar >> nType;

  CString cStr;
  size_t nCount;
  ar >> nCount;
  vRegNames.reserve(nCount);

  for(size_t i = 0; i < nCount; i++)
  {
    ar >> cStr;
    std::string sRegName((const char*)cStr);
    vRegNames.push_back(sRegName);
  }
}

//---------------------------------------------------------------------------------------
// COpenFoamFace
//---------------------------------------------------------------------------------------
bool COpenFoamFace::operator == (const COpenFoamFace& face) const
{
  if(nNodesCount != face.nNodesCount)
    return false;

  bool bFound = false;
// If any 3 nodes of the two faces coincide then the faces are equal.
  for(size_t j = 0; j < 3; j++)
  {
    bFound = false;
    size_t n0 = nodes[j];
    for(size_t i = 0; i < nNodesCount; i++)
    {
      if(n0 == face.nodes[i])
      {
        bFound = true;
        break;
      }
    }

// If at least one node of this face does not have its equivalent at the probe face, the faces are different:
    if(!bFound)
      return false;
  }

  return true;
}

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
CExportOpenFOAM::CExportOpenFOAM()
  : m_pAuxFaces(NULL)
{
  set_default();
}

CExportOpenFOAM::~CExportOpenFOAM()
{
  clear();
  clear_bound_cond();
}

void CExportOpenFOAM::set_default()
{
  m_bEnableBoundCond = false;
  m_bExportInternal = false;

  m_fDefPress = 0.4;    // Pa, this approximately corresponds to 4e-6 atm, a character press in Q00.
  set_def_temp(300.);   // K, this function sets also m_fDefDens.
  m_fDefVx = 100.;      // m/s
  m_fShiftX = 0;

  m_nSymPlanes = CAnsysMesh::spXY;
}

void CExportOpenFOAM::prepare()
{
  clear();
  terminate(false);

  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  if(!pObj->is_ready())
  {
    pObj->set_handlers(m_hJobNameHandle, m_hProgressBarHandle);
    pObj->terminate(false);
    bool bOK = pObj->read_data();
    pObj->set_handlers(NULL, NULL);
    if(!bOK)
      return;
  }

  m_pNodes = &(pObj->get_nodes());
  m_pElems = &(pObj->get_elems());

// The element centers will be used in export_internal_data() and export_boundary_data(),
// when m_pElems and m_pNodes point to elements and nodes of a different CAnsysMesh object.
  backup_element_centers();

  size_t nAuxSize = 6 * m_pElems->size();
  m_pAuxFaces = new COpenFoamFace*[nAuxSize];
  for(size_t i = 0; i < nAuxSize; i++)
    m_pAuxFaces[i] = NULL;

  m_nBoundFaceCount = 0;
  m_nUnits = euMeters;
}

void CExportOpenFOAM::clear()
{
  if(m_pAuxFaces != NULL)
  {
    size_t nAuxSize = 6 * m_pElems->size();
    for(size_t i = 0; i < nAuxSize; i++)
      delete m_pAuxFaces[i];

    delete[] m_pAuxFaces;
    m_pAuxFaces = NULL;
  }

  size_t nPatchCount = m_Patches.size();
  for(size_t j = 0; j < nPatchCount; j++)
  {
    CExportRegion* pReg = m_Patches.at(j);
    size_t nBndFaceCount = pReg->vFaces.size();
    for(size_t k = 0; k < nBndFaceCount; k++)
      delete pReg->vFaces.at(k);

    delete pReg;
  }

  m_Patches.clear();
  m_vElemCenters.clear();

  m_pNodes = NULL;
  m_pElems = NULL;
}

void CExportOpenFOAM::clear_bound_cond()
{
  size_t nBoundCondCount = m_vBoundConditions.size();
  for(size_t l = 0; l < nBoundCondCount; l++)
    delete m_vBoundConditions.at(l); 

  m_vBoundConditions.clear();
}

bool CExportOpenFOAM::do_export()
{
  prepare();

  if(!export_vertices())
    return false;

// The regions (patches) must be created before all faces are collected. Note that the owners
// of the boundary faces are still unknown.
  if(!read_2D_regions())
    return false;

// The internal faces will be placed into m_Faces, the boundary faces will be identified with those
// collected earlier in read_2D_regions() and the owner will be assigned to them.
  collect_internal_faces();

  add_boundary_faces();

// First the internal faces must be exported, then the boundary faces patch by patch.
  if(!export_faces())
    return false;

// The list of the owners must be in the order of the faces.
  if(!export_owners())
    return false;

// The list of the neighbors must be in the order of the faces.
  if(!export_neighbors())
    return false;

  merge_2D_regions();   // maybe, merging should be done only optionally?
  if(!export_boundary())
    return false;

  if(m_bEnableBoundCond)
    export_internal_and_boundary_data();

  print_bound_names();
  return true;
}

static const char sLastStr[] = { "//*******************************************************************//\n" };

bool CExportOpenFOAM::export_vertices()
{
  if(m_pNodes == NULL || m_pNodes->size() == 0)
    return false;

  std::string sFileName = m_sPath + "\\points.txt";
  
  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, sFileName.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return false;

  set_job_name("Exporting vertices...");
  set_progress(0);

// Common header.
  print_header(pStream);

// Header.
  fputs("FoamFile\n", pStream);
  fputs("{\n", pStream);
  fputs("    version     2.0;\n", pStream);
  fputs("    format      ascii;\n", pStream);
  fputs("    class       vectorField;\n", pStream);
  fputs("    location    ", pStream);
  fputc('"', pStream); fputs("constant/polyMesh", pStream); fputc('"', pStream); fputs(";\n", pStream);
  fputs("    object      points;\n", pStream);
  fputs("}\n\n", pStream);
  fputs("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n", pStream);

// Units.
  double fCoeff = m_nUnits == euMeters ? 0.01 : (m_nUnits == euCentimeters ? 1.0 : 10.0);

  size_t nNodeCount = m_pNodes->size();
  fprintf(pStream, "%zd\n(\n", nNodeCount);
  for(size_t i = 0; i < nNodeCount; i++)
  {
    const CNode3D& node = m_pNodes->at(i);
// Important: This must be the only place where m_fShiftX is used.
    fprintf(pStream, "(%f %f %f)\n", fCoeff * (node.pos.x - m_fShiftX), fCoeff * node.pos.y, fCoeff * node.pos.z);

    set_progress(int(0.5 + 100. * (i + 1) / nNodeCount));
    if(get_terminate_flag())
    {
      fclose(pStream);
      return false;
    }
  }

  fputs(")\n\n\n", pStream);
  fputs(sLastStr, pStream);

  fclose(pStream);
  return true;
}

bool CExportOpenFOAM::export_faces()
{
  std::string sFileName = m_sPath + "\\faces.txt";
  
  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, sFileName.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return false;

  set_job_name("Exporting faces...");
  set_progress(0);

// Common header.
  print_header(pStream);

// Header.
  fputs("FoamFile\n", pStream);
  fputs("{\n", pStream);
  fputs("    version     2.0;\n", pStream);
  fputs("    format      ascii;\n", pStream);
  fputs("    class       faceList;\n", pStream);
  fputs("    location    ", pStream);
  fputc('"', pStream); fputs("constant/polyMesh", pStream); fputc('"', pStream); fputs(";\n", pStream);
  fputs("    object      faces;\n", pStream);
  fputs("}\n\n", pStream);
  fputs("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n", pStream);

  size_t nFaceCount = m_Faces.size();
  fprintf(pStream, "%zd\n(\n", nFaceCount);
  for(size_t i = 0; i < nFaceCount; i++)
  {
    COpenFoamFace* pFace = m_Faces.at(i);
    size_t nNodeCount = pFace->nNodesCount;
    fprintf(pStream, "%zd", nNodeCount);
    switch(nNodeCount)
    {
      case 3: fprintf(pStream, "(%zd %zd %zd)\n", pFace->nodes[0], pFace->nodes[1], pFace->nodes[2]); break;
      case 4: fprintf(pStream, "(%zd %zd %zd %zd)\n", pFace->nodes[0], pFace->nodes[1], pFace->nodes[2], pFace->nodes[3]); break;
    }

    set_progress(int(0.5 + 100. * (i + 1) / nFaceCount));
    if(get_terminate_flag())
    {
      fclose(pStream);
      return false;
    }
  }

  fputs(")\n\n\n", pStream);
  fputs(sLastStr, pStream);

  fclose(pStream);
  return true;
}

bool CExportOpenFOAM::export_owners()
{
  std::string sFileName = m_sPath + "\\owner.txt";
  
  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, sFileName.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return false;

  set_job_name("Exporting owners...");
  set_progress(0);

// Common header.
  print_header(pStream);

// Header.
  fputs("FoamFile\n", pStream);
  fputs("{\n", pStream);
  fputs("    version     2.0;\n", pStream);
  fputs("    format      ascii;\n", pStream);
  fputs("    class       labelList;\n", pStream);

  int nNodeCount = m_pNodes->size();
  int nCellCount = m_pElems->size();
  int nFaceCount = m_Faces.size();
  int nIntCount = nFaceCount - m_nBoundFaceCount;
  fputs("    note        ", pStream); 
  fputc('"', pStream); 
  fprintf(pStream, "nPoints: %d nCells: %d nFaces: %d nInternalFaces: %d", nNodeCount, nCellCount, nFaceCount, nIntCount);
  fputc('"', pStream); fputs(";\n", pStream);

  fputs("    location    ", pStream); fputc('"', pStream); fputs("constant/polyMesh", pStream); fputc('"', pStream); fputs(";\n", pStream);
  fputs("    object      owner;\n", pStream);
  fputs("}\n\n", pStream);
  fputs("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n", pStream);

  fprintf(pStream, "%d\n(\n", nFaceCount);
  for(size_t i = 0; i < nFaceCount; i++)
  {
    COpenFoamFace* pFace = m_Faces.at(i);
    fprintf(pStream, "%d\n", pFace->nOwnerIndex);

    set_progress(int(0.5 + 100. * (i + 1) / nFaceCount));
    if(get_terminate_flag())
    {
      fclose(pStream);
      return false;
    }
  }

  fputs(")\n\n\n", pStream);
  fputs(sLastStr, pStream);

  fclose(pStream);
  return true;
}

bool CExportOpenFOAM::export_neighbors()
{
  std::string sFileName = m_sPath + "\\neighbour.txt";
  
  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, sFileName.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return false;

  set_job_name("Exporting neighbours...");
  set_progress(0);

// Common header.
  print_header(pStream);

// Header.
  fputs("FoamFile\n", pStream);
  fputs("{\n", pStream);
  fputs("    version     2.0;\n", pStream);
  fputs("    format      ascii;\n", pStream);
  fputs("    class       labelList;\n", pStream);

  int nNodeCount = m_pNodes->size();
  int nCellCount = m_pElems->size();
  int nFaceCount = m_Faces.size();
  int nIntCount = nFaceCount - m_nBoundFaceCount;
  fputs("    note        ", pStream); 
  fputc('"', pStream); 
  fprintf(pStream, "nPoints: %d nCells: %d nFaces: %d nInternalFaces: %d", nNodeCount, nCellCount, nFaceCount, nIntCount);
  fputc('"', pStream); fputs(";\n", pStream);

  fputs("    location    ", pStream); fputc('"', pStream); fputs("constant/polyMesh", pStream); fputc('"', pStream); fputs(";\n", pStream);
  fputs("    object      neighbour;\n", pStream);
  fputs("}\n\n", pStream);
  fputs("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n", pStream);

  fprintf(pStream, "%d\n(\n", nIntCount);
  for(size_t i = 0; i < nFaceCount; i++)
  {
    COpenFoamFace* pFace = m_Faces.at(i);
    if(pFace->nNbrIndex == -1)
      break;    // the internal faces are over and there will not be any in the rest of the collection.

    fprintf(pStream, "%d\n", pFace->nNbrIndex);

    set_progress(int(0.5 + 100. * (i + 1) / nFaceCount));
    if(get_terminate_flag())
    {
      fclose(pStream);
      return false;
    }
  }

  fputs(")\n\n\n", pStream);
  fputs(sLastStr, pStream);

  fclose(pStream);
  return true;
}

bool CExportOpenFOAM::export_boundary()
{
  std::string sFileName = m_sPath + "\\boundary.txt";
  
  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, sFileName.c_str(), (const char*)("w"));
  if(nErr != 0 || pStream == 0)
    return false;

  set_job_name("Exporting boundaries...");
  set_progress(0);

// Common header.
  print_header(pStream);

// Header.
  fputs("FoamFile\n", pStream);
  fputs("{\n", pStream);
  fputs("    version     2.0;\n", pStream);
  fputs("    format      ascii;\n", pStream);
  fputs("    class       polyBoundaryMesh;\n", pStream);
  fputs("    location    ", pStream);
  fputc('"', pStream); fputs("constant/polyMesh", pStream); fputc('"', pStream); fputs(";\n", pStream);
  fputs("    object      boundary;\n", pStream);
  fputs("}\n\n", pStream);
  fputs("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n", pStream);

  size_t nPatchCount = m_Patches.size();
  fprintf(pStream, "%zd\n(\n", nPatchCount);

  size_t nFirstFaceInd = m_Faces.size() - m_nBoundFaceCount;
  for(size_t i = 0; i < nPatchCount; i++)
  {
    CExportRegion* pReg = m_Patches.at(i);
    fputs(pReg->sTitle.c_str(), pStream); fputc('\n', pStream);
    fputs("{\n", pStream);

    int nType = get_bc_type(pReg);
    switch(nType)
    {
      case bcSymm:
      {
        fputs("    type        symmetryPlane;\n", pStream);
        break;
      }
      case bcPatch:
      case bcInlet:
      {
        fputs("    type        patch;\n", pStream);
        break;
      }
      case bcWall:
      {
        fputs("    type        wall;\n", pStream);
        break;
      }
    }

    size_t nFaceCount = pReg->vFaces.size();
    fprintf(pStream, "    nFaces      %zd;\n", nFaceCount);
    fprintf(pStream, "    startFace   %zd;\n", nFirstFaceInd);

    fputs("}\n", pStream);

    nFirstFaceInd += nFaceCount;

    set_progress(int(0.5 + 100. * (i + 1) / nPatchCount));
    if(get_terminate_flag())
    {
      fclose(pStream);
      return false;
    }
  }

  fputs(")\n\n\n", pStream);
  fputs(sLastStr, pStream);

  fclose(pStream);
  return true;
}

bool CExportOpenFOAM::read_2D_regions()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  std::string cBase = COutputEngine::get_base_name(pObj->get_filename());

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (cBase + "rgn").c_str(), (const char*)("r"));
  if(nErr != 0 || pStream == 0)
    return false;

  set_job_name("Reading 2D regions...");
  set_progress(0);

  char cName[32];
  int nRes, n0, n1, n2, n3;
  UINT nRegCount, nTriCount, nQuadCount, i, j;
  COpenFoamFace* pFace = NULL;

  nRes = fscanf_s(pStream, "%d", &nRegCount);
  for(i = 0; i < nRegCount; i++)
  {
    nRes = fscanf(pStream, "%s", cName);
    CExportRegion* pReg = new CExportRegion(cName);

    nRes = fscanf_s(pStream, "%d", &nTriCount);
    for(j = 0; j < nTriCount; j++)
    {
// Important! Nodes in ANSYS output files are enumerated from 1 to nNodeCount. In connectivity data, too!
// In our arrays we enumerate nodes from 0 to nNodeCount - 1, i.e. subtract unity from the original node index.
      nRes = fscanf_s(pStream, "%d %d %d", &n0, &n1, &n2);
      pFace = new COpenFoamFace(m_pNodes->at(n0 - 1).nInd, m_pNodes->at(n1 - 1).nInd, m_pNodes->at(n2 - 1).nInd);
      pReg->vFaces.push_back(pFace);
    }

    nRes = fscanf_s(pStream, "%d", &nQuadCount);
    for(j = 0; j < nQuadCount; j++)
    {
      nRes = fscanf_s(pStream, "%d %d %d %d", &n0, &n1, &n2, &n3);
      pFace = new COpenFoamFace(m_pNodes->at(n0 - 1).nInd, m_pNodes->at(n1 - 1).nInd, m_pNodes->at(n2 - 1).nInd, m_pNodes->at(n3 - 1).nInd);
      pReg->vFaces.push_back(pFace);
    }

    m_Patches.push_back(pReg);

    set_progress(int(0.5 + 100. * (i + 1) / nRegCount));
    if(get_terminate_flag())
    {
      fclose(pStream);
      return false;
    }
  }

  fclose(pStream);
  return true;
}

void CExportOpenFOAM::merge_2D_regions()
{
  size_t nRegCount = m_Patches.size();
  if(nRegCount == 0)
    return;

  set_job_name("Merging 2D regions...");
  set_progress(0);

  CPatchCollection vMergedRegions;
  CExportRegion* pPrevReg = m_Patches.at(0);
  const CBoundaryConditions* pPrevBoundCond = get_bc(pPrevReg);
  vMergedRegions.push_back(pPrevReg);   // the first region must be added in any way.
  
  UINT nWallsCount = 0;
  char sWallsCount[8];
  bool bIsPrevWall = get_bc_type(pPrevReg) == bcWall;
  if(bIsPrevWall)
  {
    nWallsCount++;
    pPrevReg->sTitle = "Wall" + std::string(_itoa(nWallsCount, sWallsCount, 10));  // rename only walls.
  }

  for(size_t i = 1; i < nRegCount; i++)
  {
    CExportRegion* pCurrReg = m_Patches.at(i);
    const CBoundaryConditions* pCurrBoundCond = get_bc(pCurrReg);
    int nBcType = pCurrBoundCond != NULL ? pCurrBoundCond->nType : bcWall;

    bool bIsCurrWall = nBcType == bcWall;
    bool bSimilarWalls = bIsCurrWall && (pCurrBoundCond == pPrevBoundCond);

    if(bSimilarWalls)
    {
      merge(pPrevReg, pCurrReg);
    }
    else if(bIsPrevWall && !bIsCurrWall)
    {
      vMergedRegions.push_back(pCurrReg);
      bIsPrevWall = bIsCurrWall;
      pPrevBoundCond = pCurrBoundCond;
      pPrevReg = pCurrReg;
    }
    else if(!bIsPrevWall && bIsCurrWall || bIsPrevWall && bIsCurrWall && !bSimilarWalls)
    {
      nWallsCount++;
      pCurrReg->sTitle = "Wall" + std::string(_itoa(nWallsCount, sWallsCount, 10));
      vMergedRegions.push_back(pCurrReg);
      bIsPrevWall = bIsCurrWall;
      pPrevBoundCond = pCurrBoundCond;
      pPrevReg = pCurrReg;
    }
    else // !bIsPrevWall && !bIsCurrWall
    {
      vMergedRegions.push_back(pCurrReg);
      bIsPrevWall = bIsCurrWall;
      pPrevBoundCond = pCurrBoundCond;
      pPrevReg = pCurrReg;
    }

    set_progress(int(0.5 + 100. * (i + 1) / nRegCount));
    if(get_terminate_flag())
      return;
  }

  m_Patches.clear();
  nRegCount = vMergedRegions.size();
  for(size_t j = 0; j < nRegCount; j++)
    m_Patches.push_back(vMergedRegions.at(j));
}

void CExportOpenFOAM::merge(CExportRegion* pDest, CExportRegion* pSrc)
{
  size_t nFacesCount = pSrc->vFaces.size();
  for(size_t i = 0; i < nFacesCount; i++)
    pDest->vFaces.push_back(pSrc->vFaces.at(i));

  pSrc->vFaces.clear();
}

void CExportOpenFOAM::collect_internal_faces()
{
  set_job_name("Collecting internal faces...");
  set_progress(0);

  size_t nElemCount = m_pElems->size();
  for(size_t i = 0; i < nElemCount; i++)
  {
    CElem3D* pElem = m_pElems->at(i);
    size_t nNodeCount = pElem->get_node_count();
    switch(nNodeCount)
    {
      case 4: // tetrahedron.
      {
        CNode3D* p0 = pElem->get_node(0);
        CNode3D* p1 = pElem->get_node(1);
        CNode3D* p2 = pElem->get_node(2);
        CNode3D* p3 = pElem->get_node(3);

        COpenFoamFace* pf0 = new COpenFoamFace(p0->nInd, p2->nInd, p1->nInd);
        process_face(pElem, pf0);

        COpenFoamFace* pf1 = new COpenFoamFace(p1->nInd, p2->nInd, p3->nInd);
        process_face(pElem, pf1);

        COpenFoamFace* pf2 = new COpenFoamFace(p0->nInd, p3->nInd, p2->nInd);
        process_face(pElem, pf2);
       
        COpenFoamFace* pf3 = new COpenFoamFace(p0->nInd, p1->nInd, p3->nInd);
        process_face(pElem, pf3);

        break;
      }
      case 5: // pyramid.
      {
        CNode3D* p0 = pElem->get_node(0);
        CNode3D* p1 = pElem->get_node(1);
        CNode3D* p2 = pElem->get_node(2);
        CNode3D* p3 = pElem->get_node(3);
        CNode3D* p4 = pElem->get_node(4);

        COpenFoamFace* pf0 = new COpenFoamFace(p0->nInd, p1->nInd, p4->nInd);
        process_face(pElem, pf0);

        COpenFoamFace* pf1 = new COpenFoamFace(p1->nInd, p2->nInd, p4->nInd);
        process_face(pElem, pf1);

        COpenFoamFace* pf2 = new COpenFoamFace(p2->nInd, p3->nInd, p4->nInd);
        process_face(pElem, pf2);
       
        COpenFoamFace* pf3 = new COpenFoamFace(p3->nInd, p0->nInd, p4->nInd);
        process_face(pElem, pf3);

        COpenFoamFace* pf4 = new COpenFoamFace(p0->nInd, p3->nInd, p2->nInd, p1->nInd);
        process_face(pElem, pf4);

        break;
      }
      case 6: // wedge or prism.
      {
        CNode3D* p0 = pElem->get_node(0);
        CNode3D* p1 = pElem->get_node(1);
        CNode3D* p2 = pElem->get_node(2);
        CNode3D* p3 = pElem->get_node(3);
        CNode3D* p4 = pElem->get_node(4);
        CNode3D* p5 = pElem->get_node(5);

        COpenFoamFace* pf0 = new COpenFoamFace(p0->nInd, p2->nInd, p1->nInd);
        process_face(pElem, pf0);

        COpenFoamFace* pf1 = new COpenFoamFace(p3->nInd, p4->nInd, p5->nInd);
        process_face(pElem, pf1);

        COpenFoamFace* pf2 = new COpenFoamFace(p0->nInd, p1->nInd, p4->nInd, p3->nInd);
        process_face(pElem, pf2);
       
        COpenFoamFace* pf3 = new COpenFoamFace(p1->nInd, p2->nInd, p5->nInd, p4->nInd);
        process_face(pElem, pf3);

        COpenFoamFace* pf4 = new COpenFoamFace(p0->nInd, p3->nInd, p5->nInd, p2->nInd);
        process_face(pElem, pf4);

        break;
      }
      case 8: // hexahedron.
      {
        CNode3D* p0 = pElem->get_node(0);
        CNode3D* p1 = pElem->get_node(1);
        CNode3D* p2 = pElem->get_node(2);
        CNode3D* p3 = pElem->get_node(3);
        CNode3D* p4 = pElem->get_node(4);
        CNode3D* p5 = pElem->get_node(5);
        CNode3D* p6 = pElem->get_node(6);
        CNode3D* p7 = pElem->get_node(7);

        COpenFoamFace* pf0 = new COpenFoamFace(p0->nInd, p1->nInd, p5->nInd, p4->nInd);
        process_face(pElem, pf0);

        COpenFoamFace* pf1 = new COpenFoamFace(p1->nInd, p3->nInd, p7->nInd, p5->nInd);
        process_face(pElem, pf1);

        COpenFoamFace* pf2 = new COpenFoamFace(p3->nInd, p2->nInd, p6->nInd, p7->nInd);
        process_face(pElem, pf2);
       
        COpenFoamFace* pf3 = new COpenFoamFace(p2->nInd, p0->nInd, p4->nInd, p6->nInd);
        process_face(pElem, pf3);

        COpenFoamFace* pf4 = new COpenFoamFace(p0->nInd, p2->nInd, p3->nInd, p1->nInd);
        process_face(pElem, pf4);

        COpenFoamFace* pf5 = new COpenFoamFace(p4->nInd, p5->nInd, p7->nInd, p6->nInd);
        process_face(pElem, pf5);

        break;
      }
    }

    set_progress(int(0.5 + 100. * (i + 1) / nElemCount));
    if(get_terminate_flag())
      return;
  }
}

bool CExportOpenFOAM::process_face(CElem3D* pElem, COpenFoamFace* pFace)
{
  int nNbrIndex = find_neighbour(pElem, pFace);

  if(nNbrIndex == -1) // there is no neighbour => pFace is a boundary face.
  {
    COpenFoamFace* pBndFace = identify_boundary_face(pFace);
    if(pBndFace == NULL)
      return false;

    correct_normal(pElem, pBndFace);  // correct the sense of rotation in the face if necessary.

    pBndFace->nOwnerIndex = pElem->nInd;
    pBndFace->vFaceCenter = get_face_center(pBndFace);
    m_nBoundFaceCount++;
    delete pFace;
    return true;
  }

  size_t nAuxInd = 6 * nNbrIndex;  // approximate location of pFace in the auxiliary array m_pAuxFaces.

  COpenFoamFace* pExistingFace = m_pAuxFaces[nAuxInd];
  while(pExistingFace != NULL)
  {
    if(*pExistingFace == *pFace)
      break;
    else
    {
      nAuxInd++;
      pExistingFace = m_pAuxFaces[nAuxInd];
    }
  }

  if(pExistingFace != NULL) // this face already exists in the collection, and pElem is its neighbour.
  {
    pExistingFace->nNbrIndex = pElem->nInd;
    delete pFace;
  }
  else  // this is a new face, insert it into the collection and assign pElem->nInd to it as its owner index.
  {
    correct_normal(pElem, pFace); // correct the sense of rotation in the face if necessary.

    pFace->nOwnerIndex = pElem->nInd;
    pFace->vFaceCenter = get_face_center(pFace);  // this vector will be needed during exporting boundary values.
    m_Faces.push_back(pFace);

// Find a proper place for the new face in the auxiliary array:
    nAuxInd = 6 * pElem->nInd;
    pExistingFace = m_pAuxFaces[nAuxInd];
    while(pExistingFace != NULL)
    {
      nAuxInd++;
      pExistingFace = m_pAuxFaces[nAuxInd];
    }

    m_pAuxFaces[nAuxInd] = pFace;
  }

  return true;
}

int CExportOpenFOAM::find_neighbour(CElem3D* pElem, COpenFoamFace* pFace)
{
  const CElementsCollection& vElems = *m_pElems;
// It is enough to try any three vertices of the face. If all three of them belong to an element, 
// this element is what we are looking for.
  CNode3D* p0 = &(m_pNodes->at(pFace->nodes[0]));
  CNode3D* p1 = &(m_pNodes->at(pFace->nodes[1]));
  CNode3D* p2 = &(m_pNodes->at(pFace->nodes[2]));

// The element to be found must be among vNbrElems of any vertex of pFace. Try the neighbours of p0.
  size_t nNbrCount = p0->vNbrElems.size();
  for(size_t i = 0; i < nNbrCount; i++)
  {
    CElem3D* pNbrElem = vElems[p0->vNbrElems.at(i)];
    if(pNbrElem == pElem) // vNbrElems contains also pElem, skip it. 
      continue;

    if(node_in_elem(p0, pNbrElem) && node_in_elem(p1, pNbrElem) && node_in_elem(p2, pNbrElem))
      return pNbrElem->nInd;
  }

  return -1;  // the face hasn't got a neighbour => it is a boundary face.
}

bool CExportOpenFOAM::node_in_elem(CNode3D* pNode, CElem3D* pElem)
{
  size_t nCount = pElem->get_node_count();
  for(size_t i = 0; i < nCount; i++)
    if(pElem->get_node(i) == pNode)
      return true;

  return false;
}

COpenFoamFace* CExportOpenFOAM::identify_boundary_face(COpenFoamFace* pFace)
{
  size_t nRegCount = m_Patches.size();
  for(size_t i = 0; i < nRegCount; i++)
  {
    CExportRegion* pReg = m_Patches.at(i);
    size_t nFaceCount = pReg->vFaces.size();
    for(size_t j = 0; j < nFaceCount; j++)
    {
      if(*pFace == *(pReg->vFaces.at(j)))
        return pReg->vFaces.at(j);
    }
  }

  return NULL;
}

void CExportOpenFOAM::add_boundary_faces()
{
// Add the boundary faces in the order as they are stored in m_Patches:
  size_t nRegCount = m_Patches.size();
  for(size_t j = 0; j < nRegCount; j++)
  {
    CExportRegion* pReg = m_Patches.at(j);
    size_t nBndFaceCount = pReg->vFaces.size();
    for(size_t k = 0; k < nBndFaceCount; k++)
    {
      COpenFoamFace* pFace = pReg->vFaces.at(k);
      m_Faces.push_back(pFace);
    }
  }
}

void CExportOpenFOAM::correct_normal(CElem3D* pOwnerElem, COpenFoamFace* pFace)
{
// Normal must point always to the cell with the larger index, i.e. from the owner to the neighbour.
// In the case of the boundary face the normal must point from the domain, i.e. again from the owner.
  Vector3D vCellCenter(0, 0, 0);
  size_t nNodeElemCount = pOwnerElem->get_node_count();
  for(size_t i = 0; i < nNodeElemCount; i++)
    vCellCenter += pOwnerElem->get_node(i)->pos;

  vCellCenter /= (double)nNodeElemCount;

  Vector3D vFaceCenter(0, 0, 0);
  size_t j, nFaceNodeCount = pFace->nNodesCount;
  for(j = 0; j < nFaceNodeCount; j++)
    vFaceCenter += m_pNodes->at(pFace->nodes[j]).pos;

  vFaceCenter /= (double)nFaceNodeCount;

  Vector3D vDir = (vFaceCenter - vCellCenter).normalized(); // direction out of the cell.
  Vector3D vFaceNorm = (m_pNodes->at(pFace->nodes[1]).pos - m_pNodes->at(pFace->nodes[0]).pos) * (m_pNodes->at(pFace->nodes[2]).pos - m_pNodes->at(pFace->nodes[0]).pos);
  double fDotProd = vDir & vFaceNorm;
  if(fDotProd > 0)
    return; // everything is OK, the normal points from the owner cell.

// Change the order of the nodes in the face:
  size_t vTmp[4];
  for(j = 0; j < nFaceNodeCount; j++)
    vTmp[j] = pFace->nodes[j];

  for(j = 0; j < nFaceNodeCount; j++)
    pFace->nodes[j] = vTmp[nFaceNodeCount - j - 1];
}

void CExportOpenFOAM::print_header(FILE* pStream)
{
  fputs("/*--------------------------------*- C++ -*----------------------------------*\\", pStream); fputc('\n', pStream);
  fputs("| =========                 |                                                 |", pStream); fputc('\n', pStream);
  fputs("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |", pStream); fputc('\n', pStream);
  fputs("|  \\\\    /   O peration     | Version:  2.3.0                                 |", pStream); fputc('\n', pStream);
  fputs("|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |", pStream); fputc('\n', pStream);
  fputs("|    \\\\/     M anipulation  |                                                 |", pStream); fputc('\n', pStream);
  fputs("\\*---------------------------------------------------------------------------*/", pStream); fputc('\n', pStream);
}

// Boundary conditions.
void CExportOpenFOAM::export_internal_and_boundary_data()
{
  if(!PathFileExists(m_sBoundCondPath.c_str()))
    return;

// Instead of using a temporary CTracker object, we re-initialize the global CTracker.
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();

// Backup the data filename of the current CTracker.
  CString sFileName(pObj->get_filename());
  bool bConv2CGS = pObj->get_convert_to_cgs();
  int nSymPlanes = pObj->get_sym_plane();

  pObj->set_convert_to_cgs(false); // make sure that all the variables are in SI.

  pObj->clear();
  pObj->set_sym_plane(m_nSymPlanes);
  pObj->set_filename(m_sBoundCondPath.c_str());
  pObj->set_handlers(m_hJobNameHandle, m_hProgressBarHandle);
  if(!pObj->read_geometry() || !pObj->read_gasdyn_data())
    return;

  std::string sOutName = m_sPath + "\\boundaryT.dat";
  FILE* pFileT;
  errno_t nErr = fopen_s(&pFileT, sOutName.c_str(), (const char*)("w"));
  if(nErr != 0 || pFileT == 0)
    return;

  sOutName = m_sPath + "\\boundaryU.dat";
  FILE* pFileU;
  nErr = fopen_s(&pFileU, sOutName.c_str(), (const char*)("w"));
  if(nErr != 0 || pFileU == 0)
  {
    fclose(pFileT);
    return;
  }

  sOutName = m_sPath + "\\rhoN.dat";
  FILE* pFileN;
  nErr = fopen_s(&pFileN, sOutName.c_str(), (const char*)("w"));
  if(nErr != 0 || pFileN == 0)
  {
    fclose(pFileT);
    fclose(pFileU);
    return;
  }

  print_header_T(pFileT);
  print_header_U(pFileU);
  print_header_N(pFileN);

  if(m_bExportInternal)
    export_internal_data(pObj, pFileT, pFileU, pFileN);
  else
    write_uniform_internal_fields(pFileT, pFileU, pFileN);

  export_boundary_data(pObj, pFileT, pFileU, pFileN);

  fclose(pFileT);
  fclose(pFileU);
  fclose(pFileN);

// Restore the global CTracker to its original condition:
  pObj->clear();
  pObj->set_sym_plane(nSymPlanes);
  pObj->set_filename((const char*)sFileName);
  pObj->set_convert_to_cgs(bConv2CGS);
  pObj->read_data();
  pObj->set_handlers(NULL, NULL);
}

void CExportOpenFOAM::write_uniform_internal_fields(FILE* pFileT, FILE* pFileU, FILE* pFileN)
{
// Temperature
  fprintf(pFileT, "internalField uniform %f;\n\n", m_fDefTemp);

// Velocity
  Vector3D vVel(m_fDefVx, 0, 0);
  fprintf(pFileU, "internalField uniform (%f %f %f);\n\n", vVel.x, vVel.y, vVel.z);

// Number density
  fprintf(pFileN, "internalField uniform %e;\n\n", m_fDefDen);
}

void CExportOpenFOAM::backup_element_centers()
{
  Vector3D vNull(0, 0, 0);
  size_t nElemCount = m_pElems->size();
  m_vElemCenters.resize(nElemCount, vNull);
  for(size_t i = 0; i < nElemCount; i++)
    m_vElemCenters[i] = get_elem_center(m_pElems->at(i));
}

void CExportOpenFOAM::export_internal_data(CTracker* pObj, FILE* pFileT, FILE* pFileU, FILE* pFileN)
{
  set_job_name("Exporting internal data...");
  set_progress(0);

// DEBUG
  bool bAuxOK = true;
  std::string sAuxName = m_sPath + "\\__internal_report.txt";
  FILE* pAuxFile;
  errno_t nErr = fopen_s(&pAuxFile, sAuxName.c_str(), (const char*)("w"));
  if(nErr != 0 || pAuxFile == 0)
    bAuxOK = false;

  if(bAuxOK)
    fputs("Element number  Element center (mm)\n\n", pAuxFile);
// END DEBUG

  size_t nElemCount = m_vElemCenters.size();  // elements count in the DSMC mesh.

// Temperature
  fputs("internalField nonuniform List<scalar>\n", pFileT);
  fprintf(pFileT, "%zd\n", nElemCount);
  fputs("(\n", pFileT);
// Velocity
  fputs("internalField nonuniform List<vector>\n", pFileU);
  fprintf(pFileU, "%zd\n", nElemCount);
  fputs("(\n", pFileU);
// Number density
  fputs("internalField nonuniform List<scalar>\n", pFileN);
  fprintf(pFileN, "%zd\n", nElemCount);
  fputs("(\n", pFileN);

  Vector3D vPos, vVel, vAccel;
  const CElem3D* pSearchElem = NULL;
  CNode3D node;   // this is just a container for interpolated data.

  CBox box = pObj->get_box();
  double fNumDens;
  for(size_t i = 0; i < nElemCount; i++)
  {
    vPos = m_vElemCenters.at(i);
    CSymCorrData data = pObj->sym_corr_forward(vPos, vVel);

    if(!box.inside(vPos) || !pObj->interpolate(vPos, 0, 0, node, pSearchElem) && bAuxOK)
    {
// DEBUG
      fprintf(pAuxFile, "%zd      (%f %f %f)\n", i, 10*vPos.x, 10*vPos.y, 10*vPos.z);
// END DEBUG
      node.temp = m_fDefTemp;
      node.vel = Vector3D(m_fDefVx, 0, 0);
      fNumDens = m_fDefDen;
    }
    else
    {
      fNumDens = node.press / (node.temp * Const_Boltzmann_SI);
    }

    vVel = node.vel;
    pObj->sym_corr_back(vPos, vVel, vAccel, data);
    node.vel = vVel;

    fprintf(pFileT, "%f\n", node.temp);
    fprintf(pFileU, "(%f %f %f)\n", node.vel.x, node.vel.y, node.vel.z);
    fprintf(pFileN, "%e\n", fNumDens);

    set_progress(int(0.5 + 100. * (i + 1) / nElemCount));
    if(get_terminate_flag())
      return;
  }

  fputs(");\n\n", pFileT);
  fputs(");\n\n", pFileU);
  fputs(");\n\n", pFileN);

// DEBUG
  if(bAuxOK)
    fclose(pAuxFile);
// END DEBUG
}

void CExportOpenFOAM::export_boundary_data(CTracker* pObj, FILE* pFileT, FILE* pFileU, FILE* pFileN)
{
  set_job_name("Exporting boundary data...");
  set_progress(0);

  fputs("boundaryField\n", pFileT); fputs("{\n", pFileT);
  fputs("boundaryField\n", pFileU); fputs("{\n", pFileU);
  fputs("boundaryField\n", pFileN); fputs("{\n", pFileN);

  CBox box = pObj->get_box();

  Vector3D vFaceCenter, vElemCenter, vPos, vVel, vAccel, vDefVel(m_fDefVx, 0, 0);

  size_t nElemCount = m_vElemCenters.size();
  size_t nPatchCount = m_Patches.size();

  for(size_t i = 0; i < nPatchCount; i++)
  {
    CExportRegion* pReg = m_Patches.at(i);
    size_t nFaceCount = pReg->vFaces.size();

    set_progress(int(0.5 + 100. * (i + 1) / nPatchCount));
    if(get_terminate_flag())
      return;

    CBoundaryConditions* pBC = get_bc(pReg);
    int nType = (pBC != NULL) ? pBC->nType : bcWall;

    double fNumDens, fTemp;
    switch(nType)
    {
      case bcSymm:
      {
        print_symm_name(pReg->sName.c_str(), pFileT);
        print_symm_name(pReg->sName.c_str(), pFileU);
        print_symm_name(pReg->sName.c_str(), pFileN);
        continue;
      }
      case bcPatch:
      case bcInlet:
      {
        print_region_name(pReg->sName.c_str(), nFaceCount, true, pFileT);
        print_region_name(pReg->sName.c_str(), nFaceCount, false, pFileU);
        print_region_name(pReg->sName.c_str(), nFaceCount, true, pFileN);

        CNode3D node;   // this is just a container for interpolated data.
        const CElem3D* pSearchElem = NULL;
        for(size_t j = 0; j < nFaceCount; j++)
        {
          COpenFoamFace* pFace = pReg->vFaces.at(j);
          vFaceCenter = pFace->vFaceCenter;

          vElemCenter = m_vElemCenters.at(pFace->nOwnerIndex);
          vPos = vFaceCenter + 0.1 * (vElemCenter - vFaceCenter);  // shift a bit from the center of the boundary face inside the volume.
          CSymCorrData data = pObj->sym_corr_forward(vPos, vVel);

          if(!box.inside(vPos) || !pObj->interpolate(vPos, 0, 0, node, pSearchElem))
          {
            vPos = vElemCenter;
            data = pObj->sym_corr_forward(vPos, vVel);
            if(!box.inside(vPos) || !pObj->interpolate(vPos, 0, 0, node, pSearchElem))
            {
              node.temp = m_fDefTemp;
              node.vel = Vector3D(m_fDefVx, 0, 0);
              fNumDens = m_fDefDen;
            }
          }
          else
          {
            fNumDens = node.press / (node.temp * Const_Boltzmann_SI);
          }

          vVel = node.vel;
          pObj->sym_corr_back(vPos, vVel, vAccel, data);
          node.vel = vVel;

          if(nType == bcInlet)
          {
            fprintf(pFileT, "    %f\n", node.temp);
            fprintf(pFileU, "    (%f %f %f)\n", node.vel.x, node.vel.y, node.vel.z);
            fprintf(pFileN, "    %e\n", fNumDens);
          }
          else  // bcPatch.
          {
            fprintf(pFileT, "    %f\n", pBC->fTemp);
            fprintf(pFileU, "    (%f %f %f)\n", node.vel.x, node.vel.y, node.vel.z);
            fNumDens = (pBC->fTemp > 1.0) ? pBC->fPress / (pBC->fTemp * Const_Boltzmann_SI) : m_fDefDen;
            fprintf(pFileN, "    %e\n", fNumDens);
          }
        }
      
        fputs("    );\n", pFileT);
        fputs("    );\n", pFileU);
        fputs("    );\n", pFileN);
        break;
      }
      case bcWall:
      {
        fTemp = (pBC != NULL) ? pBC->fTemp : m_fDefTemp;
        fprintf(pFileT, "  %s\n  {\n    type fixedValue;\n    value uniform %f;\n", pReg->sTitle.c_str(), fTemp);
        fprintf(pFileU, "  %s\n  {\n    type fixedValue;\n    value uniform (%f %f %f);\n", pReg->sTitle.c_str(), 0., 0., 0.);
        fprintf(pFileN, "  %s\n  {\n    type calculated;\n    value uniform %d;\n", pReg->sTitle.c_str(), 0);
        break;
      }
    }

    fputs("  ", pFileT); fputs("}\n", pFileT);
    fputs("  ", pFileU); fputs("}\n", pFileU);
    fputs("  ", pFileN); fputs("}\n", pFileN);
  }

  fputs("}\n\n", pFileT); fputs(sLastStr, pFileT);
  fputs("}\n\n", pFileU); fputs(sLastStr, pFileU);
  fputs("}\n\n", pFileN); fputs(sLastStr, pFileN);
}

CBoundaryConditions* CExportOpenFOAM::get_bc(CExportRegion* pReg) const
{
  const std::string& sRegName = pReg->sName;
  size_t nBoundCondCount = m_vBoundConditions.size();
  for(size_t i = 0; i < nBoundCondCount; i++)
  {
    CBoundaryConditions* pBoundCond = m_vBoundConditions.at(i);
    size_t nRegCount = pBoundCond->vRegNames.size();
    for(size_t j = 0; j < nRegCount; j++)
    {
      const std::string& sName = pBoundCond->vRegNames.at(j);
      if(sRegName == sName)
        return pBoundCond;
    }
  }

  return NULL;
}

int CExportOpenFOAM::get_bc_type(CExportRegion* pReg) const
{
  const CBoundaryConditions* pBoundCond = get_bc(pReg);
  if(pBoundCond != NULL)
    return pBoundCond->nType;

  return bcWall;
}

void CExportOpenFOAM::add_bound_cond(CBoundaryConditions* pBC)
{
  m_vBoundConditions.push_back(pBC);
}

void CExportOpenFOAM::remove_bound_cond(CBoundaryConditions* pBoundCond)
{
  int nPos = -1;
  size_t nBoundCondCount = m_vBoundConditions.size();
  for(size_t i = 0; i < nBoundCondCount; i++)
  {
    CBoundaryConditions* pBC = m_vBoundConditions.at(i);
    if(pBC == pBoundCond)
    {
      nPos = i;
      break;
    }
  }

  if(nPos >= 0)
    m_vBoundConditions.erase(m_vBoundConditions.begin() + nPos);
}

Vector3D CExportOpenFOAM::get_face_center(COpenFoamFace* pFace) const
{
  Vector3D vC(0, 0, 0);
  size_t nNodeCount = pFace->nNodesCount;  // can be either 3 or 4.
  for(size_t k = 0; k < nNodeCount; k++)
    vC += m_pNodes->at(pFace->nodes[k]).pos;

  return (vC / (double)nNodeCount);
}

Vector3D CExportOpenFOAM::get_elem_center(const CElem3D* pElem) const
{
  Vector3D vC(0, 0, 0);
  size_t nNodeCount = pElem->get_node_count();
  for(size_t k = 0; k < nNodeCount; k++)
    vC += pElem->get_node(k)->pos;

  return (vC / (double)nNodeCount);
}

void CExportOpenFOAM::print_header_T(FILE* pOutFile)
{
// Common header.
  print_header(pOutFile);

// Header.
  fputs("FoamFile\n", pOutFile);
  fputs("{\n", pOutFile);
  fputs("    version     2.0;\n", pOutFile);
  fputs("    format      ascii;\n", pOutFile);
  fputs("    class       volScalarField;\n", pOutFile);
  fputs("    object      T;\n", pOutFile);
  fputs("}\n\n", pOutFile);
  fputs("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n", pOutFile);
  fputs("dimensions [0 0 0 1 0 0 0];\n\n", pOutFile);
}

void CExportOpenFOAM::print_header_U(FILE* pOutFile)
{
// Common header.
  print_header(pOutFile);

// Header.
  fputs("FoamFile\n", pOutFile);
  fputs("{\n", pOutFile);
  fputs("    version     2.0;\n", pOutFile);
  fputs("    format      ascii;\n", pOutFile);
  fputs("    class       volVectorField;\n", pOutFile);
  fputs("    object      U;\n", pOutFile);
  fputs("}\n\n", pOutFile);
  fputs("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n", pOutFile);
  fputs("dimensions [0 1 -1 0 0 0 0];\n\n", pOutFile);
}

void CExportOpenFOAM::print_header_N(FILE* pOutFile)
{
// Common header.
  print_header(pOutFile);

// Header.
  fputs("FoamFile\n", pOutFile);
  fputs("{\n", pOutFile);
  fputs("    version     2.0;\n", pOutFile);
  fputs("    format      ascii;\n", pOutFile);
  fputs("    class       volScalarField;\n", pOutFile);
  fputs("    object      rhoN;\n", pOutFile);
  fputs("}\n\n", pOutFile);
  fputs("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n", pOutFile);
  fputs("dimensions [0 -3 0 0 0 0 0];\n\n", pOutFile);
}

void CExportOpenFOAM::print_region_name(const char* pName, size_t nCount, bool bScalar, FILE* pOutFile)
{
  fputs("  ", pOutFile); fputs(pName, pOutFile); fputc('\n', pOutFile);
  fputs("  ", pOutFile); fputs("{\n", pOutFile);
  fputs("    ", pOutFile); fputs("type fixedValue;", pOutFile); fputc('\n', pOutFile);

  fputs("    ", pOutFile); fputs("value nonuniform List", pOutFile);
  if(bScalar)
    fputs("<scalar> ", pOutFile);
  else
    fputs("<vector> ", pOutFile);

  fprintf(pOutFile, "%zd\n", nCount);
  fputs("    (\n", pOutFile);
}

void CExportOpenFOAM::print_symm_name(const char* pName, FILE* pOutFile)
{
  fputs("  ", pOutFile); fputs(pName, pOutFile); fputc('\n', pOutFile);
  fputs("  ", pOutFile); fputs("{\n", pOutFile);
  fputs("  ", pOutFile); fputs("type symmetryPlane;", pOutFile); fputc('\n', pOutFile);
  fputs("  ", pOutFile); fputs("}\n", pOutFile);
}

void CExportOpenFOAM::print_bound_names()
{
  std::string sOutName = m_sPath + "\\boundary_names.txt";
  FILE* pFile;
  errno_t nErr = fopen_s(&pFile, sOutName.c_str(), (const char*)("w"));
  if(nErr != 0 || pFile == 0)
    return;

  size_t nBcCount = m_vBoundConditions.size(), nRegCount;
  for(size_t i = 0; i < nBcCount; i++)
  {
    const CBoundaryConditions* pBC = m_vBoundConditions.at(i);
    switch(pBC->nType)
    {
      case bcInlet:
      {
        fputs("Inlet:\n", pFile);
        break;
      }
      case bcPatch:
      {
        fprintf(pFile, "Patch:   P = %f,  T = %f\n", pBC->fPress, pBC->fTemp);
        break;
      }
      case bcSymm:
      {
        fputs("Symmetry:\n", pFile);
        break;
      }
      case bcWall:
      {
        fprintf(pFile, "Wall:   T = %f\n", pBC->fTemp);
        break;
      }
    }

    nRegCount = pBC->vRegNames.size();
    for(size_t j = 0; j < nRegCount; j++)
    {
      std::string sName = pBC->vRegNames.at(j);
      if(j == nRegCount - 1)
        fprintf(pFile, "%s\n", sName.c_str());
      else
        fprintf(pFile, "%s, ", sName.c_str());
    }

    fputc('\n', pFile);
  }

  fclose(pFile);
}

//-------------------------------------------------------------------------------------------------
// Streaming:
//-------------------------------------------------------------------------------------------------
void CExportOpenFOAM::save(CArchive& ar)
{
  const UINT nVersion = 1;  // CBoundaryConditions.
  ar << nVersion;

// Output path:
  CString cFileName(m_sPath.c_str());
  ar << cFileName;

// Boundary conditions path (if any)
  CString cBoundFileName(m_sBoundCondPath.c_str());
  ar << cBoundFileName;

  ar << m_bEnableBoundCond;
  ar << m_fShiftX;
  ar << m_fDefPress;
  ar << m_fDefTemp;
  ar << m_fDefVx;

  size_t nBoundCondCount = m_vBoundConditions.size();
  ar << nBoundCondCount;
  for(size_t i = 0; i < nBoundCondCount; i++)
  {
    CBoundaryConditions* pBC = m_vBoundConditions.at(i);
    pBC->save(ar);
  }
}

void CExportOpenFOAM::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CString cFileName;
  ar >> cFileName;
  set_path((const char*)cFileName);

  CString cBoundFileName;
  ar >> cBoundFileName;
  set_boundary_cond_file((const char*)cBoundFileName);

  ar >> m_bEnableBoundCond;
  ar >> m_fShiftX;
  ar >> m_fDefPress;
  ar >> m_fDefTemp;
  ar >> m_fDefVx;

  clear_bound_cond();

  if(nVersion < 1)
  {
// Boundary conditions of symmetry type:
    CBoundaryConditions* pSymmCond = new CBoundaryConditions(0., 0., bcSymm);

    size_t nSymmCount;
    ar >> nSymmCount;
    pSymmCond->vRegNames.reserve(nSymmCount);
    for(size_t i = 0; i < nSymmCount; i++)
    {
      CString cSymmName;
      ar >> cSymmName;
      std::string sName((const char*)cSymmName);
      pSymmCond->vRegNames.push_back(sName);
    }

    m_vBoundConditions.push_back(pSymmCond);

// Boundary conditions of patch type:
    CBoundaryConditions* pPatchCond = new CBoundaryConditions(m_fDefPress, m_fDefTemp, bcPatch);

    size_t nPatchCount;
    ar >> nPatchCount;
    pPatchCond->vRegNames.reserve(nPatchCount);
    for(size_t j = 0; j < nPatchCount; j++)
    {
      CString cPatchName;
      ar >> cPatchName;
      std::string sName((const char*)cPatchName);
      pPatchCond->vRegNames.push_back(sName);
    }

    m_vBoundConditions.push_back(pPatchCond);
  }
  else
  {
    size_t nBoundCondCount;
    ar >> nBoundCondCount;
    m_vBoundConditions.reserve(nBoundCondCount);
    for(size_t i = 0; i < nBoundCondCount; i++)
    {
      CBoundaryConditions* pBC = new CBoundaryConditions();
      pBC->load(ar);
      m_vBoundConditions.push_back(pBC);
    }
  }
}


//-------------------------------------------------------------------------------------------------
// CExternalGridExport - export to grid nodes
//-------------------------------------------------------------------------------------------------
bool CExternalGridExport::do_export()
{
  FILE* pInFile = open_input_file();
  if(pInFile == NULL)
    return false;

  set_job_name("Exporting ANSYS Data...");
  set_progress(0);

  CLocVector vPoints;
  read_locations(pInFile, vPoints);
  if(vPoints.size() == 0)
    return false;

  fclose(pInFile);

  FILE* pOutFile = open_output_file();
  if(pOutFile == NULL)
    return false;

  interpolate_and_write_data(pOutFile, vPoints);

  fclose(pOutFile);
  set_progress(100);
  return true;
}

void CExternalGridExport::read_locations(FILE* pInFile, CLocVector& vLoc)
{
  int nRes = 0;
  Vector3D vPos;
  double x, y, z;

  switch(m_nFmtType)
  {
    case fmtCSV:
    {
      int nLocCount = 0;
      int nRes = fscanf_s(pInFile, "%d", &nLocCount);
      if(nRes != 1)
        return;

      for(int i = 0; i < nLocCount; i++)
      {
        nRes = fscanf_s(pInFile, "%lf, %lf, %lf", &vPos.x, &vPos.y, &vPos.z);
        if(nRes != 3)
          continue;

// I expect that if user set m_bUseSI to "true" the coordinates are in meters. 
// Transfer them into CGS for future use in the interpolation procedure.
        if(m_bUseSI)
          vPos *= SI_to_CGS_Len;

        vLoc.push_back(vPos);
      }
      break;
    }
    case fmtDAT4:
    {
      int n, nCount = 0;
      while(nRes != EOF)
      {
        nRes = fscanf_s(pInFile, "%d %lf %lf %lf", &n, &vPos.x, &vPos.y, &vPos.z);
        if(nRes != 4)
          continue;

// I expect that if user set m_bUseSI to "true" the coordinates are in meters. 
// Transfer them into CGS for future use in the interpolation procedure.
        if(m_bUseSI)
          vPos *= SI_to_CGS_Len;

        vLoc.push_back(vPos);
      }
      break;
    }
    case fmtDAT3:
    {
      int n, nCount = 0;
      while(nRes != EOF)
      {
        nRes = fscanf_s(pInFile, "%lf %lf %lf", &vPos.x, &vPos.y, &vPos.z);
        if(nRes != 3)
          continue;

// I expect that if user set m_bUseSI to "true" the coordinates are in meters. 
// Transfer them into CGS for future use in the interpolation procedure.
        if(m_bUseSI)
          vPos *= SI_to_CGS_Len;

        vLoc.push_back(vPos);
      }
      break;
    }
  }
}

static const Vector3F vNull(0, 0, 0);

static const char sccHeader[] = 
"    Node,        X,            Y,            Z,           Vx,           Vy,           Vz,        Press,       Density,        Temp,           Cp,       Therm Cond,     Dyn Visc\n";

void CExternalGridExport::interpolate_and_write_data(FILE* pOutFile, const CLocVector& vLoc)
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  fputs(sccHeader, pOutFile);

// Interpolation errors handling:
  std::string cPath = COutputEngine::get_full_path(pObj->get_filename());
  CString sFailedLocFile = CString(cPath.c_str()) + CString("failed_points.csv");

  FILE* pErrFile;
  errno_t nErr = fopen_s(&pErrFile, (const char*)sFailedLocFile, (const char*)("w"));
  bool bErrFileOk = (nErr == 0) && (pErrFile != NULL);

  CNode3D node;
  const CElem3D* pElem = NULL;
  Vector3D vPos, vVel, vAccel;
  size_t nLocCount = vLoc.size();
  for(size_t i = 0; i < nLocCount; i++)
  {
    set_progress(int(100 * double(i) / nLocCount));
    if(get_terminate_flag())
      return;

    vPos = vLoc.at(i);
    CSymCorrData data = pObj->sym_corr_forward(vPos, vVel);
    if(pObj->interpolate(vPos, 0, 0, node, pElem))
    {
      vVel = node.vel;
      pObj->sym_corr_back(vPos, vVel, vAccel, data);
      node.pos = vPos;
      node.vel = vVel;
    }
    else
    {
      node.pos = vPos;
      node.set_data(0, 0, 0, 0, 0, 0, vNull, vNull, vNull);
      if(bErrFileOk)
        fprintf(pErrFile, "%8d\n", (int)(i+1));   // write down the node numbers in which interpolation failed.
    }

    if(m_bUseSI)
      convert_to_SI(node);

    switch(m_nFmtType)
    {
      case fmtCSV:
      {
        fprintf(pOutFile, "%8d, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e, %12.5e\n", 
          (int)(i+1), node.pos.x, node.pos.y, node.pos.z, node.vel.x, node.vel.y, node.vel.z, node.press, node.dens, node.temp, node.cp, node.cond, node.visc);
        break;
      }
      case fmtDAT4:
      case fmtDAT3:
      {
        fprintf(pOutFile, "%8d  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n", 
          (int)(i+1), node.pos.x, node.pos.y, node.pos.z, node.vel.x, node.vel.y, node.vel.z, node.press, node.dens, node.temp, node.cp, node.cond, node.visc);
        break;
      }
    }
  }

  if(bErrFileOk)
    fclose(pErrFile);
}

void CExternalGridExport::convert_to_SI(CNode3D& node) const
{
  node.pos *= CGS_to_SI_Len;
  node.vel *= CGS_to_SI_Len;

  node.press *= CGS_to_SI_Press;
  node.dens  *= CGS_to_SI_Dens;
  node.cp    /= SI_to_CGS_Cp;
  node.cond  *= CGS_to_SI_ThermCond;
  node.visc  *= CGS_to_SI_DynVisc;
}

FILE* CExternalGridExport::open_input_file() const
{
  FILE* pFile;
  errno_t nErr = fopen_s(&pFile, (const char*)m_sInputFile, (const char*)("r"));
  if(nErr != 0 || pFile == 0)
  {
    AfxMessageBox("Can not open input file for reading data.");
    return NULL;
  }

  return pFile;
}

FILE* CExternalGridExport::open_output_file() const
{
  FILE* pFile;
  errno_t nErr = fopen_s(&pFile, (const char*)m_sOutputFile, (const char*)("w"));
  if(nErr != 0 || pFile == 0)
  {
    AfxMessageBox("Can not open output file for writing.");
    return NULL;
  }

  return pFile;
}

const char* CExternalGridExport::get_format_name(int nFmtType)
{
  switch(nFmtType)
  {
    case fmtCSV: return _T(" CSV (x, y, z)");
    case fmtDAT4: return _T(" DAT-4 (n x y z)");
    case fmtDAT3: return _T(" DAT-3 (x y z)");
  }
  return _T("");
}

};