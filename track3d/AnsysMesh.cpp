
#include "stdafx.h"

#include "AnsysMesh.h"
#include "ParticleTracking.h"

namespace EvaporatingParticle
{
//-------------------------------------------------------------------------------------------------
// CAnsysMesh - the base class for mesh data reading from the ANSYS data file
//-------------------------------------------------------------------------------------------------
CAnsysMesh::CAnsysMesh()
{
  m_bReady = false;
  m_bConv2CGS = true;
  m_nSymPlanes = spXY | spXZ;
}

CAnsysMesh::~CAnsysMesh()
{
  clear();
}

void CAnsysMesh::clear()
{
  size_t nElemCount = m_vElems.size();
  for(size_t j = 0; j < nElemCount; j++)
  {
    set_status("Deleting elements", 100 * j / nElemCount);
    delete m_vElems.at(j);
  }

  size_t nRegCount = m_vRegions.size();
  for(size_t k = 0; k < nRegCount; k++)
  {
    CRegion* pReg = m_vRegions.at(k);
    size_t nFaceCount = pReg->vFaces.size();
    for(size_t l = 0; l < nFaceCount; l++)
      delete pReg->vFaces.at(l);

    set_status("Deleting regions", 100 * k / nRegCount);
    delete pReg;
  }

  size_t nNodeCount = m_vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    set_status("Deleting nodes", 100 * i / nNodeCount);
    delete m_vNodes.at(i);
  }

  m_vNodes.clear();
  m_vElems.clear();
  m_vRegions.clear();
  m_vExtRegions.clear();

  m_Transform.set_default();

  CTrackDraw* pDrawObj = CParticleTrackingApp::Get()->GetDrawObj();
  pDrawObj->clear();
  pDrawObj->invalidate_all();

  set_status("Ready", -1);
}

//-------------------------------------------------------------------------------------------------
//  Data reading methods
//-------------------------------------------------------------------------------------------------
bool CAnsysMesh::read_data()
{
  clear();

  if(!read_geometry())
    return false;

  if(!read_gasdyn_data()) // read the gas-dynamic variables in the nodes of the mesh.
    return false; // wrong gas-dynamic data format.

  if(!read_2D_regions())
    return false;

  invalidate_calculators();
  m_bReady = true;
  return true;
}

bool CAnsysMesh::read_geometry()
{
  set_job_name("Reading points...");
  set_progress(0);

  std::string cBase = COutputEngine::get_base_name(get_filename());

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (cBase + "geom").c_str(), (const char*)("r"));
  if(nErr != 0 || pStream == 0)
    return false;

  char cStr[32];
  int nNode, nNodeCount, nRes = 0;
  float fX, fY, fZ;

  nRes = fscanf_s(pStream, "%s", cStr, 32);
  nRes = fscanf_s(pStream, "%d", &nNodeCount);
// Nodes:
  m_vNodes.clear();
  m_vNodes.reserve(nNodeCount);

  CNode3D* pNode = NULL;
// Important! Nodes in ANSYS output are enumerated from 1 to nNodeCount. In connectivity data, too!
// In our arrays we enumerate nodes from 0 to nNodeCount - 1, i.e. subtract unity from the original node index.
  for(UINT i = 0; i < nNodeCount; i++)
  {
    nRes = fscanf_s(pStream, "%d %e %e %e", &nNode, &fX, &fY, &fZ);
    if(nRes == EOF)
      return false;

    fX *= SI_to_CGS_Len;
    fY *= SI_to_CGS_Len;
    fZ *= SI_to_CGS_Len;

    Vector3D vPos(fX, fY, fZ);
// Mesh transformation:
    if(m_Transform.get_enable())
      m_Transform.transform(vPos);

    pNode = new CNode3D(vPos);
    pNode->nInd = i;
    m_vNodes.push_back(pNode);

// Progress bar and termination support:
    if(i % 100 == 0)
      set_progress(int(0.5 + 100. * i / nNodeCount));
    if(m_bTerminate)
      return abort(pStream);
  }

// Connectivity information.
  m_vElems.clear();

// Tetrahedra (4 nodes, 4 planes):
  int nTetraCount, n0, n1, n2, n3;
  nRes = fscanf_s(pStream, "%s", cStr, 32);
  nRes = fscanf_s(pStream, "%d", &nTetraCount);
  m_vElems.reserve(nTetraCount);

  if(nTetraCount > 0)
  {
    set_job_name("Adding tetra...");

    for(UINT i = 0; i < nTetraCount; i++)
    {
      nRes = fscanf_s(pStream, "%d %d %d %d", &n0, &n1, &n2, &n3);
      if(nRes == EOF)
        return false;

      add_tetra(m_vNodes.at(n0 - 1), m_vNodes.at(n1 - 1), m_vNodes.at(n2 - 1), m_vNodes.at(n3 - 1));

// Progress bar and termination support:
      if(i % 100 == 0)
        set_progress(int(0.5 + 100. * i / nTetraCount));
      if(m_bTerminate)
        return abort(pStream);
    }
  }

// Pyramids (5 nodes, 5 planes):
  int nPyrCount, n4;
  nRes = fscanf_s(pStream, "%s", cStr, 32);
  nRes = fscanf_s(pStream, "%d", &nPyrCount);
  m_vElems.reserve(nTetraCount + nPyrCount);

  if(nPyrCount > 0)
  {
    set_job_name("Adding pyramids...");

    for(UINT i = 0; i < nPyrCount; i++)
    {
      nRes = fscanf_s(pStream, "%d %d %d %d %d", &n0, &n1, &n2, &n3, &n4);
      if(nRes == EOF)
        return false;

      add_pyramid(m_vNodes.at(n0 - 1), m_vNodes.at(n1 - 1), m_vNodes.at(n2 - 1), m_vNodes.at(n3 - 1), m_vNodes.at(n4 - 1));

// Progress bar and termination support:
      if(i % 100 == 0)
        set_progress(int(0.5 + 100. * i / nPyrCount));
      if(m_bTerminate)
        return abort(pStream);
    }
  }

// Wedges (or prisms, 6 nodes, 5 planes):
  int nWedgeCount, n5;
  nRes = fscanf_s(pStream, "%s", cStr, 32);
  nRes = fscanf_s(pStream, "%d", &nWedgeCount);
  m_vElems.reserve(nTetraCount + nPyrCount + nWedgeCount);

  if(nWedgeCount > 0)
  {
    set_job_name("Adding wedges...");

    for(UINT i = 0; i < nWedgeCount; i++)
    {
      nRes = fscanf_s(pStream, "%d %d %d %d %d %d", &n0, &n1, &n2, &n3, &n4, &n5);
      if(nRes == EOF)
        return false;

      add_wedge(m_vNodes.at(n0 - 1), m_vNodes.at(n1 - 1), m_vNodes.at(n2 - 1), m_vNodes.at(n3 - 1), m_vNodes.at(n4 - 1), m_vNodes.at(n5 - 1));

// Progress bar and termination support:
      if(i % 100 == 0)
        set_progress(int(0.5 + 100. * i / nWedgeCount));
      if(m_bTerminate)
        return abort(pStream);
    }
  }

// Hexahedra (8 nodes, 6 planes):
  int nHexaCount, n6, n7;
  nRes = fscanf_s(pStream, "%s", cStr, 32);
  nRes = fscanf_s(pStream, "%d", &nHexaCount);
  m_vElems.reserve(nTetraCount + nPyrCount + nWedgeCount + nHexaCount);

  if(nHexaCount > 0)
  {
    set_job_name("Adding hexa...");

    for(UINT i = 0; i < nHexaCount; i++)
    {
      nRes = fscanf_s(pStream, "%d %d %d %d %d %d %d %d", &n0, &n1, &n2, &n3, &n4, &n5, &n6, &n7);
      if(nRes == EOF)
        return false;

      add_hexa(m_vNodes.at(n0 - 1), m_vNodes.at(n1 - 1), m_vNodes.at(n3 - 1), m_vNodes.at(n2 - 1), 
               m_vNodes.at(n4 - 1), m_vNodes.at(n5 - 1), m_vNodes.at(n7 - 1), m_vNodes.at(n6 - 1));

// Progress bar and termination support:
      if(i % 100 == 0)
        set_progress(int(0.5 + 100. * i / nHexaCount));
      if(m_bTerminate)
        return abort(pStream);
    }
  }

  for(UINT i = 0; i < nNodeCount; i++)
    m_vNodes.at(i)->shrink_to_fit();

  fclose(pStream);
  bounding_box();
  return true;
}

bool CAnsysMesh::read_2D_regions()
{
  set_job_name("Reading 2D regions...");
  set_progress(0);

  std::string cBase = COutputEngine::get_base_name(get_filename());

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (cBase + "rgn").c_str(), (const char*)("r"));
  if(nErr != 0 || pStream == 0)
    return false;

  char cName[32];
  int nRes, n0, n1, n2, n3;
  UINT nRegCount, nTriCount, nQuadCount, i, j, nFace;
  CFace* pFace = NULL;

  nRes = fscanf_s(pStream, "%d", &nRegCount);
  for(i = 0; i < nRegCount; i++)
  {
    nRes = fscanf(pStream, "%s", cName);
    CRegion* pReg = new CRegion(cName);

    nRes = fscanf_s(pStream, "%d", &nTriCount);
    for(j = 0; j < nTriCount; j++)
    {
      nRes = fscanf_s(pStream, "%d %d %d", &n0, &n1, &n2);
      n0--;
      n1--;
      n2--;
      pFace = new CFace(m_vNodes.at(n0), m_vNodes.at(n1), m_vNodes.at(n2));
      pReg->vFaces.push_back(pFace);
      nFace = pReg->vFaces.size() - 1;
      CRegFacePair pair(i, nFace);
      m_vNodes.at(n0)->vNbrFaces.push_back(pair);
      m_vNodes.at(n1)->vNbrFaces.push_back(pair);
      m_vNodes.at(n2)->vNbrFaces.push_back(pair);
    }

    nRes = fscanf_s(pStream, "%d", &nQuadCount);
    for(j = 0; j < nQuadCount; j++)
    {
      nRes = fscanf_s(pStream, "%d %d %d %d", &n0, &n1, &n2, &n3);
      n0--;
      n1--;
      n2--;
      n3--;
      pFace = new CFace(m_vNodes.at(n0), m_vNodes.at(n1), m_vNodes.at(n2));
      pReg->vFaces.push_back(pFace);
      nFace = pReg->vFaces.size() - 1;
      CRegFacePair pair1(i, nFace);
      m_vNodes.at(n0)->vNbrFaces.push_back(pair1);
      m_vNodes.at(n1)->vNbrFaces.push_back(pair1);
      m_vNodes.at(n2)->vNbrFaces.push_back(pair1);

      pFace = new CFace(m_vNodes.at(n0), m_vNodes.at(n2), m_vNodes.at(n3));
      pReg->vFaces.push_back(pFace);
      nFace = pReg->vFaces.size() - 1;
      CRegFacePair pair2(i, nFace);
      m_vNodes.at(n0)->vNbrFaces.push_back(pair2);
      m_vNodes.at(n2)->vNbrFaces.push_back(pair2);
      m_vNodes.at(n3)->vNbrFaces.push_back(pair2);
    }

    pReg->bounding_box();
    m_vRegions.push_back(pReg);

// Progress bar and termination support:
    if(i % 100 == 0)
      set_progress(int(0.5 + 100. * i / nRegCount));
    if(m_bTerminate)
      return abort(pStream);
  }

  fclose(pStream);
  return true;
}

bool CAnsysMesh::read_gasdyn_data(bool bFieldsOnly)
{
  set_job_name("Reading gas-dynamic data...");
  set_progress(0);

  std::string cBase = COutputEngine::get_base_name(get_filename());;

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, (cBase + "var").c_str(), (const char*)("r"));
  if(nErr != 0 || pStream == 0)
    return false;

  int nRes = 0;
  char sHeader[16];
  for(UINT j = 0; j < 24; j++)
    nRes = fscanf(pStream, "%s", sHeader);

  float fPress, fDens, fTemp, fDynVisc, fThermCond, fCp, fVx, fVy, fVz, fEx, fEy, fEz, fRFEx, fRFEy, fRFEz;

  CNode3D* pNode = NULL;
  Vector3D vVel, vFieldDC, vFieldRF;
  size_t nNodeCount = m_vNodes.size();
  for(size_t i = 0; i < nNodeCount; i++)
  {
    nRes = fscanf_s(pStream, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",
      &fPress, &fDens, &fTemp, &fDynVisc, &fThermCond, &fCp, &fVx, &fVy, &fVz, &fEx, &fEy, &fEz, &fRFEx, &fRFEy, &fRFEz);

    if(nRes == EOF)
      return false;

    if(m_bConv2CGS)
      conv_to_cgs(fPress, fDens, fDynVisc, fThermCond, fCp, fVx, fVy, fVz, fEx, fEy, fEz, fRFEx, fRFEy, fRFEz);

    vVel = Vector3D(fVx, fVy, fVz);
    vFieldDC = Vector3D(fEx, fEy, fEz);
    vFieldRF = Vector3D(fRFEx, fRFEy, fRFEz);
// Mesh transformation:
    if(m_Transform.get_enable())
    {
      m_Transform.transform(vVel, false);
      m_Transform.transform(vFieldDC, false);
      m_Transform.transform(vFieldRF, false);
    }

    pNode = m_vNodes.at(i);
    if(bFieldsOnly) // [MS] 29-06-2017 workaround to allow ANSYS fields together with DSMC gas dynamics.
    {
      pNode->field = vFieldDC;
      pNode->rf = vFieldRF;
    }
    else
    {
      pNode->set_data(fPress, fDens, fTemp, fDynVisc, fThermCond, fCp, vVel, vFieldDC, vFieldRF);
    }

// Progress bar and termination support:
    if(i % 100 == 0)
      set_progress(int(0.5 + 100. * i / nNodeCount));
    if(m_bTerminate)
      return abort(pStream);
  }

  fclose(pStream);
  return true;
}

void CAnsysMesh::add_tetra(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3)
{
  CTetra* pTetra = new CTetra(m_vNodes, p0->nInd, p1->nInd, p2->nInd, p3->nInd);
  pTetra->nInd = m_vElems.size(); // this index is used in integrators as well as in export of the mesh to the OpenFOAM format.
  m_vElems.push_back((CElem3D*)pTetra);

  p0->vNbrElems.push_back(pTetra->nInd);
  p1->vNbrElems.push_back(pTetra->nInd);
  p2->vNbrElems.push_back(pTetra->nInd);
  p3->vNbrElems.push_back(pTetra->nInd);
}

void CAnsysMesh::add_pyramid(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4)
{
  CPyramid* pPyr = new CPyramid(m_vNodes, p0->nInd, p1->nInd, p2->nInd, p3->nInd, p4->nInd);
  pPyr->nInd = m_vElems.size();
  m_vElems.push_back((CElem3D*)pPyr);

  p0->vNbrElems.push_back(pPyr->nInd);
  p1->vNbrElems.push_back(pPyr->nInd);
  p2->vNbrElems.push_back(pPyr->nInd);
  p3->vNbrElems.push_back(pPyr->nInd);
  p4->vNbrElems.push_back(pPyr->nInd);
}

void CAnsysMesh::add_wedge(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4, CNode3D* p5)
{
  CWedge* pWedge = new CWedge(m_vNodes, p0->nInd, p1->nInd, p2->nInd, p3->nInd, p4->nInd, p5->nInd);
  pWedge->nInd = m_vElems.size();
  m_vElems.push_back((CElem3D*)pWedge);

  p0->vNbrElems.push_back(pWedge->nInd);
  p1->vNbrElems.push_back(pWedge->nInd);
  p2->vNbrElems.push_back(pWedge->nInd);
  p3->vNbrElems.push_back(pWedge->nInd);
  p4->vNbrElems.push_back(pWedge->nInd);
  p5->vNbrElems.push_back(pWedge->nInd);
}

void CAnsysMesh::add_hexa(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4, CNode3D* p5, CNode3D* p6, CNode3D* p7)
{
  CHexa* pHexa = new CHexa(m_vNodes, p0->nInd, p1->nInd, p2->nInd, p3->nInd, p4->nInd, p5->nInd, p6->nInd, p7->nInd);
  pHexa->nInd = m_vElems.size();
  m_vElems.push_back((CElem3D*)pHexa);

  p0->vNbrElems.push_back(pHexa->nInd);
  p1->vNbrElems.push_back(pHexa->nInd);
  p2->vNbrElems.push_back(pHexa->nInd);
  p3->vNbrElems.push_back(pHexa->nInd);
  p4->vNbrElems.push_back(pHexa->nInd);
  p5->vNbrElems.push_back(pHexa->nInd);
  p6->vNbrElems.push_back(pHexa->nInd);
  p7->vNbrElems.push_back(pHexa->nInd);
}

void CAnsysMesh::bounding_box()
{
  size_t nNodesCount = m_vNodes.size();
  if(nNodesCount == 0)
    return;

  m_Box.vMin = m_vNodes.at(0)->pos;
  m_Box.vMax = m_vNodes.at(0)->pos;
  
  for(size_t i = 1; i < nNodesCount; i++)
  {
    CNode3D* pNode = m_vNodes.at(i);

    if(pNode->pos.x < m_Box.vMin.x)
      m_Box.vMin.x = pNode->pos.x;
    if(pNode->pos.x > m_Box.vMax.x)
      m_Box.vMax.x = pNode->pos.x;

    if(pNode->pos.y < m_Box.vMin.y)
      m_Box.vMin.y = pNode->pos.y;
    if(pNode->pos.y > m_Box.vMax.y)
      m_Box.vMax.y = pNode->pos.y;

    if(pNode->pos.z < m_Box.vMin.z)
      m_Box.vMin.z = pNode->pos.z;
    if(pNode->pos.z > m_Box.vMax.z)
      m_Box.vMax.z = pNode->pos.z;
  }
}

void CAnsysMesh::conv_to_cgs(float& fPress, float& fDens, float& fDynVisc, float& fThermCond, float& fCp,
  float& fVx, float& fVy, float& fVz, float& fEx, float& fEy, float& fEz, float& fRFEx, float& fRFEy, float& fRFEz)
{
  fPress *= (float)SI_to_CGS_Press;
  fDens *= (float)SI_to_CGS_Dens;
  fDynVisc *= (float)SI_to_CGS_DynVisc;
  fThermCond *= (float)SI_to_CGS_ThermCond;

  fEx *= (float)SI_to_CGS_ElecField;
  fEy *= (float)SI_to_CGS_ElecField;
  fEz *= (float)SI_to_CGS_ElecField;

  fRFEx *= (float)SI_to_CGS_ElecField;
  fRFEy *= (float)SI_to_CGS_ElecField;
  fRFEz *= (float)SI_to_CGS_ElecField;

  fVx *= (float)SI_to_CGS_Vel;
  fVy *= (float)SI_to_CGS_Vel;
  fVz *= (float)SI_to_CGS_Vel;

  fCp *= (float)SI_to_CGS_Cp;
}

//-------------------------------------------------------------------------------------------------
//  Searching elements through the mesh methods
//-------------------------------------------------------------------------------------------------
CElem3D* CAnsysMesh::find_elem(CElem3D* pPrevElem, const Vector3D& vPos) const
{
  CElem3D* pElem = NULL;
// First, try the nearest neighbors of the previous face (if it is not zero), including the previous face itself:
  if(pPrevElem != NULL)
  {
    pElem = try_neighbors(pPrevElem, vPos);
    if(pElem != NULL)
      return pElem;
  }

  size_t nCount = m_vElems.size();  // run over all the elements collection:
  for(size_t i = 0; i < nCount; i++)
  {
    pElem = m_vElems.at(i);
    if(pElem->inside(vPos))
      return pElem;
  }

  return NULL;
}

CElem3D* CAnsysMesh::try_neighbors(CElem3D* pElem, const Vector3D& vPos) const
{
// First, try the input element itself.
  if(pElem->inside(vPos))
    return pElem;

  UINT nNbrInd = 0;
  CElem3D* pNbrElem = NULL;
// Try the neighbors of pElem.
  CIndexVector vUsedElems;
  size_t nNodeCount = pElem->get_node_count();
  for(size_t j = 0; j < nNodeCount; j++)
  {
    CNode3D* pNode = pElem->get_node(j);
    size_t nNbrCount = pNode->vNbrElems.size();
    for(size_t i = 0; i < nNbrCount; i++)
    {
      nNbrInd = pNode->vNbrElems.at(i);
      if(std::find(vUsedElems.begin(), vUsedElems.end(), nNbrInd) != vUsedElems.end())
        continue;

      pNbrElem = m_vElems.at(nNbrInd);
      if(pNbrElem->inside(vPos))
        return pNbrElem;

      vUsedElems.push_back(nNbrInd);
    }
  }

  return NULL;
}

//-------------------------------------------------------------------------------------------------
// Reflection of positions and velocities due to symmetry. If any of XY or XZ symmetry planes exist
// POSITIVE z and y are supposed to be inside the domain.
//-------------------------------------------------------------------------------------------------
bool CAnsysMesh::sym_corr(Vector3D& vPos, Vector3D& vVel, Vector3D& vAccel, bool bForceReflect) const
{
  bool bReflDone = false;
  if((m_nSymPlanes & spXY) && (vPos.z < 0 || bForceReflect))
  {
    vPos.z = -vPos.z;
    vVel.z = -vVel.z;
    vAccel.z = -vAccel.z;
    bReflDone = true;
  }

  if((m_nSymPlanes & spXZ) && (vPos.y < 0 || bForceReflect))
  {
    vPos.y = -vPos.y;
    vVel.y = -vVel.y;
    vAccel.y = -vAccel.y;
    bReflDone = true;
  }

  if((m_nSymPlanes & spYZ) && (vPos.x < 0 || bForceReflect))
  {
    vPos.x = -vPos.x;
    vVel.x = -vVel.x;
    vAccel.x = -vAccel.x;
    bReflDone = true;
  }

  return bReflDone;
}

//-------------------------------------------------------------------------------------------------
// Streaming:
//-------------------------------------------------------------------------------------------------
void CAnsysMesh::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  CString cFileName(m_sDataFile.c_str());
  ar << cFileName;

  ar << m_nSymPlanes;

  m_Transform.save(ar);
}

void CAnsysMesh::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  CString cFileName;
  ar >> cFileName;
  set_filename((const char*)cFileName);

  ar >> m_nSymPlanes;

  m_Transform.load(ar);
}

//-------------------------------------------------------------------------------------------------
// Auxiliary functions:
//-------------------------------------------------------------------------------------------------
CRegionsCollection& CAnsysMesh::get_regions()
{
  CCrossSectColl* pColl = CParticleTrackingApp::Get()->GetPlanes();
  size_t nPlanesCount = pColl->size();
  if(nPlanesCount == 0)
    return m_vRegions;

  m_vExtRegions.clear();
  size_t nRegCount = m_vRegions.size();
  size_t nExtRegCount = nRegCount + nPlanesCount;
  m_vExtRegions.reserve(nExtRegCount);

  for(size_t i = 0; i < nRegCount; i++)
    m_vExtRegions.push_back(m_vRegions.at(i));

  for(size_t j = 0; j < nPlanesCount; j++)
  {
    CDomainCrossSection* pPlane = pColl->at(j);
    CRegion* pReg = pPlane->get_region();
    if(pReg != NULL)
      m_vExtRegions.push_back(pReg);
  }

  return m_vExtRegions;
}

void CAnsysMesh::invalidate_calculators()
{
  CParticleTrackingApp::Get()->GetCalcs()->invalidate_calcs();
}

bool CAnsysMesh::abort(FILE* pStream)
{
  if(pStream != NULL)
    fclose(pStream);

  m_bReady = false;

  return false;
}

}; // namespace EvaporatingParticle
