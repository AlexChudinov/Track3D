
#include "stdafx.h"

#include "float.h"
#include "DomainCrossSection.h"
#include "ParticleTracking.h"


namespace EvaporatingParticle
{

//---------------------------------------------------------------------------------------
//  CCrossSectColl
//---------------------------------------------------------------------------------------
void CCrossSectColl::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  size_t nCount = size();
  ar << nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    CDomainCrossSection* pObj = at(i);
    pObj->save(ar);
  }
}

void CCrossSectColl::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  size_t nCount;
  ar >> nCount;

  for(size_t i = 0; i < nCount; i++)
  {
    CDomainCrossSection* pObj = new CDomainCrossSection();
    pObj->load(ar);
    push_back(pObj);
  }
}

//---------------------------------------------------------------------------------------
//  CDomainCrossSection
//---------------------------------------------------------------------------------------
CDomainCrossSection::CDomainCrossSection(const Vector3D& vOrigin, const Vector3D& vNorm)
: m_Region(_T("Cross-Section # ")), m_Plane(vOrigin, vNorm), m_bReady(false)
{
  set_default();
}

CDomainCrossSection::~CDomainCrossSection()
{
  clear();
}

void CDomainCrossSection::set_default()
{
  set_plane_type(ptPlaneYZ);
  m_Plane.pos.x = 0.1;

  m_Region.bEnabled = true;
  m_Region.bCrossSection = true;

  char buff[4];
  CCrossSectColl* pColl = CParticleTrackingApp::Get()->GetPlanes();
  size_t nId = pColl->size() + 1;
  m_Region.sName += std::string((const char*)itoa(nId, buff, 10));
}

void CDomainCrossSection::clear()
{
  size_t nNodeCount = m_vNodes.size();
  for(size_t j = 0; j < nNodeCount; j++)
    delete m_vNodes.at(j);

  m_vNodes.clear();

  size_t nFaceCount = m_Region.vFaces.size();
  for(size_t i = 0; i < nFaceCount; i++)
    delete m_Region.vFaces.at(i);

  m_Region.vFaces.clear();
}

void CDomainCrossSection::build_mesh()
{
  clear();
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  CElementsCollection& elems = pObj->get_elems();

  CElem3D* pElem = NULL;
  CFacesCollection vFaces;
  size_t nElemCount = elems.size();
  for(size_t i = 0; i < nElemCount; i++)
  {
    pElem = elems.at(i);
    if(!intersect_element(vFaces, pElem))
      continue;

    size_t nFaceCount = vFaces.size();
    for(size_t j = 0; j < nFaceCount; j++)
      m_Region.vFaces.push_back(vFaces.at(j));

    vFaces.clear();
  }

  if(m_Region.vFaces.size() > 0)
  {
    m_Region.bounding_box();
    m_bReady = true;
  }
}

bool CDomainCrossSection::intersect_element(CFacesCollection& vFaces, CElem3D* pElem)
{
  const size_t nNodeCount = pElem->get_node_count();
  if(nNodeCount == 0)
    return false;

// First, inquire whether or not the plane intersects the element.
  bool bIntersect = false;
  int n1, n0 = m_Plane.inside(pElem->get_node(0)->pos) ? 1 : -1;
  for(size_t i = 1; i < nNodeCount; i++)
  {
    n1 = m_Plane.inside(pElem->get_node(i)->pos) ? 1 : -1;
    if(n0 * n1 < 0)
    {
      bIntersect = true;
      break;
    }
  }

  if(!bIntersect)
    return false;

// Second, process elements depending on their types:
  switch(nNodeCount)
  {
    case 4: return intersect_tetra(vFaces, (CTetra*)pElem);
    case 5: return intersect_piramid(vFaces, (CPyramid*)pElem);
    case 6: return intersect_wedge(vFaces, (CWedge*)pElem);
    case 8: return intersect_hexa(vFaces, (CHexa*)pElem);
  }

  return false;
}

static const UINT scnTetEdgeCount = 6;
static const CTetraEdge scTetEdges[scnTetEdgeCount] =
  { CTetraEdge(0, 1, 3), CTetraEdge(1, 2, 4), CTetraEdge(2, 0, 5), CTetraEdge(2, 3, 0), CTetraEdge(0, 3, 1), CTetraEdge(1, 3, 2) };

static CNode3D* svTetraSection[4];

bool CDomainCrossSection::intersect_tetra(CFacesCollection& vFaces, CTetra* pElem)
{
  Vector3D vRes;
  CNode3D* pNode = NULL;
  double fEdgeLen, fDist, fKsi;
  CRay ray(cvDefOrigin, cvDefNorm);

  std::vector<UINT> vEdgeIds;
  vEdgeIds.reserve(scnTetEdgeCount);
  for(UINT k = 0; k < scnTetEdgeCount; k++)
    vEdgeIds.push_back(k);

  int nLastIntersect = -1;
  size_t nEdgeCount = vEdgeIds.size();  // variable count of the edges to be processed.
  UINT nId, nNodeCount = 0, nAttemptsCount = 0;
  while((nEdgeCount != 0) && (nAttemptsCount < scnTetEdgeCount))
  {
    for(UINT i = 0; i < nEdgeCount; i++)
    {
      nId = vEdgeIds.at(i);
      if((nLastIntersect >= 0) && (nLastIntersect == scTetEdges[nId].nExcl))
        continue;   // skip this edge, but do not erase its id from the vEdgeIds vector - it will be tried later.

      ray.orig = pElem->get_node(scTetEdges[nId].n0)->pos;
      ray.dir = pElem->get_node(scTetEdges[nId].n1)->pos - ray.orig;
      fEdgeLen = ray.dir.length();
      if(fEdgeLen < Const_Almost_Zero)
        return false;

      ray.dir /= fEdgeLen;
      if(!m_Plane.intersect(ray, vRes, fDist) || fDist > fEdgeLen)
      {
        vEdgeIds.erase(vEdgeIds.begin() + i); // erase the processed edge.
        break;
      }

      fKsi = fDist / fEdgeLen;  // I expect fKsi to be 0 <= fKsi <= 1.
      pNode = interpolate(vRes, pElem->get_node(scTetEdges[nId].n0), pElem->get_node(scTetEdges[nId].n1), fKsi);
  // Remember the last intersection edge index.
      nLastIntersect = nId;
      svTetraSection[nNodeCount] = pNode;
      nNodeCount++;

      vEdgeIds.erase(vEdgeIds.begin() + i); // erase the processed edge.
      break;
    }

    nEdgeCount = vEdgeIds.size();
    nAttemptsCount++;
  }

  if(nNodeCount < 3)
    return false;

  CFace* pFace = new CFace(svTetraSection[0], svTetraSection[1], svTetraSection[2]);
  vFaces.push_back(pFace);
  if(nNodeCount == 3)
    return true;

  pFace = new CFace(svTetraSection[0], svTetraSection[2], svTetraSection[3]);
  vFaces.push_back(pFace);
  return true;
}

static const UINT scnPyrEdgeCount = 8;
static const CPyramidEdge scPyrEdges[scnPyrEdgeCount] = 
  { CPyramidEdge(0, 1, 4, 5), CPyramidEdge(1, 2, 5, 6), CPyramidEdge(2, 3, 6, 7), CPyramidEdge(3, 0, 4, 7),
    CPyramidEdge(2, 4, 0, 3), CPyramidEdge(3, 4, 0, 1), CPyramidEdge(0, 4, 1, 2), CPyramidEdge(1, 4, 2, 3) };

static CNode3D* svPyrSection[5];

bool CDomainCrossSection::intersect_piramid(CFacesCollection& vFaces, CPyramid* pElem)
{
  Vector3D vRes;
  CNode3D* pNode = NULL;
  double fEdgeLen, fDist, fKsi;
  CRay ray(cvDefOrigin, cvDefNorm);

  std::vector<UINT> vEdgeIds;
  vEdgeIds.reserve(scnPyrEdgeCount);
  for(UINT k = 0; k < scnPyrEdgeCount; k++)
    vEdgeIds.push_back(k);

  int nLastIntersect = -1;
  size_t nEdgeCount = vEdgeIds.size();  // variable count of the edges to be processed.
  UINT nId, nNodeCount = 0, nAttemptsCount = 0;
  while((nEdgeCount != 0) && (nAttemptsCount < scnPyrEdgeCount))
  {
    for(UINT i = 0; i < nEdgeCount; i++)
    {
      nId = vEdgeIds.at(i);
      if((nLastIntersect >= 0) && (nLastIntersect == scPyrEdges[nId].nExcl1 || nLastIntersect == scPyrEdges[nId].nExcl2))
        continue;   // skip this edge, but do not erase its id from the vEdgeIds vector - it will be tried later.

      ray.orig = pElem->get_node(scPyrEdges[nId].n0)->pos;
      ray.dir = pElem->get_node(scPyrEdges[nId].n1)->pos - ray.orig;
      fEdgeLen = ray.dir.length();
      if(fEdgeLen < Const_Almost_Zero)
        return false;

      ray.dir /= fEdgeLen;
      if(!m_Plane.intersect(ray, vRes, fDist) || fDist > fEdgeLen)
      {
        vEdgeIds.erase(vEdgeIds.begin() + i); // erase the processed edge.
        break;
      }

      fKsi = fDist / fEdgeLen;  // I expect fKsi to be 0 <= fKsi <= 1.
      pNode = interpolate(vRes, pElem->get_node(scPyrEdges[nId].n0), pElem->get_node(scPyrEdges[nId].n1), fKsi);
  // Remember the last intersection edge index.
      nLastIntersect = nId;
      svPyrSection[nNodeCount] = pNode;
      nNodeCount++;

      vEdgeIds.erase(vEdgeIds.begin() + i); // erase the processed edge.
      break;
    }

    nEdgeCount = vEdgeIds.size();
    nAttemptsCount++;
  }

  if(nNodeCount < 3)
    return false;

  CFace* pFace = new CFace(svPyrSection[0], svPyrSection[1], svPyrSection[2]);
  vFaces.push_back(pFace);
  if(nNodeCount == 3)
    return true;

  pFace = new CFace(svPyrSection[0], svPyrSection[2], svPyrSection[3]);
  vFaces.push_back(pFace);
  if(nNodeCount == 4)
    return true;

  pFace = new CFace(svPyrSection[0], svPyrSection[3], svPyrSection[4]);
  vFaces.push_back(pFace);
  return true;
}

static const UINT scnWedgeEdgeCount = 9;
static const CWedgeEdge scWedgeEdges[scnWedgeEdgeCount] =
{ CWedgeEdge(0, 1, 3, 7, 8), CWedgeEdge(1, 2, 4, 6, 8), CWedgeEdge(2, 0, 5, 6, 7),
  CWedgeEdge(2, 5, 0, 6, 9), CWedgeEdge(0, 3, 1, 7, 9), CWedgeEdge(1, 4, 2, 8, 9),    // 9 is a dummy index, these edges have 2 exclusions only.
  CWedgeEdge(3, 4, 1, 2, 3), CWedgeEdge(4, 5, 0, 2, 4), CWedgeEdge(5, 3, 0, 1, 5) };

static CNode3D* svWedgeSection[5];

bool CDomainCrossSection::intersect_wedge(CFacesCollection& vFaces, CWedge* pElem)
{
  Vector3D vRes;
  CNode3D* pNode = NULL;
  double fEdgeLen, fDist, fKsi;
  CRay ray(cvDefOrigin, cvDefNorm);

  std::vector<UINT> vEdgeIds;
  vEdgeIds.reserve(scnWedgeEdgeCount);
  for(UINT k = 0; k < scnWedgeEdgeCount; k++)
    vEdgeIds.push_back(k);

  int nLastIntersect = -1;
  size_t nEdgeCount = vEdgeIds.size();  // variable count of the edges to be processed.
  UINT nId, nNodeCount = 0, nAttemptsCount = 0;
  while((nEdgeCount != 0) && (nAttemptsCount < scnWedgeEdgeCount))
  {
    for(UINT i = 0; i < nEdgeCount; i++)
    {
      nId = vEdgeIds.at(i);
      bool bExcl = nLastIntersect == scWedgeEdges[nId].nExcl1 || nLastIntersect == scWedgeEdges[nId].nExcl2 || nLastIntersect == scWedgeEdges[nId].nExcl3;
      if((nLastIntersect >= 0) && bExcl)
        continue;   // skip this edge, but do not erase its id from the vEdgeIds vector - it will be tried later.

      ray.orig = pElem->get_node(scWedgeEdges[nId].n0)->pos;
      ray.dir = pElem->get_node(scWedgeEdges[nId].n1)->pos - ray.orig;
      fEdgeLen = ray.dir.length();
      if(fEdgeLen < Const_Almost_Zero)
        return false;

      ray.dir /= fEdgeLen;
      if(!m_Plane.intersect(ray, vRes, fDist) || fDist > fEdgeLen)
      {
        vEdgeIds.erase(vEdgeIds.begin() + i); // erase the processed edge.
        break;
      }

      fKsi = fDist / fEdgeLen;  // I expect fKsi to be 0 <= fKsi <= 1.
      pNode = interpolate(vRes, pElem->get_node(scWedgeEdges[nId].n0), pElem->get_node(scWedgeEdges[nId].n1), fKsi);
  // Remember the last intersection edge index.
      nLastIntersect = nId;
      svWedgeSection[nNodeCount] = pNode;
      nNodeCount++;

      vEdgeIds.erase(vEdgeIds.begin() + i); // erase the processed edge.
      break;
    }

    nEdgeCount = vEdgeIds.size();
    nAttemptsCount++;
  }

  if(nNodeCount < 3)
    return false;

  CFace* pFace = new CFace(svWedgeSection[0], svWedgeSection[1], svWedgeSection[2]);
  vFaces.push_back(pFace);
  if(nNodeCount == 3)
    return true;

  pFace = new CFace(svWedgeSection[0], svWedgeSection[2], svWedgeSection[3]);
  vFaces.push_back(pFace);
  if(nNodeCount == 4)
    return true;

  pFace = new CFace(svWedgeSection[0], svWedgeSection[3], svWedgeSection[4]);
  vFaces.push_back(pFace);
  return true;
}

static const UINT scnHexaEdgeCount = 12;
static const CHexaEdge scHexaEdges[scnHexaEdgeCount] =
{ CHexaEdge(0, 1, 4, 5, 9, 10, 11), CHexaEdge(1, 2, 5, 6, 8, 10, 11), CHexaEdge(2, 3, 6, 7, 8,  9, 11), CHexaEdge(3, 0, 4, 7, 8, 9, 10),
  CHexaEdge(2, 6, 0, 3, 6,  8, 11), CHexaEdge(3, 7, 0, 1, 7,  8,  9), CHexaEdge(0, 4, 1, 2, 4,  9, 10), CHexaEdge(1, 5, 2, 3, 5, 10, 11),
  CHexaEdge(4, 5, 1, 2, 3,  4,  5), CHexaEdge(5, 6, 0, 2, 3,  5,  6), CHexaEdge(6, 7, 0, 1, 3,  6,  7), CHexaEdge(7, 4, 0, 1, 2,  4,  7) };

static CNode3D* svHexaSection[6];

bool CDomainCrossSection::intersect_hexa(CFacesCollection& vFaces, CHexa* pElem)
{
  Vector3D vRes;
  CNode3D* pNode = NULL;
  double fEdgeLen, fDist, fKsi;
  CRay ray(cvDefOrigin, cvDefNorm);

  std::vector<UINT> vEdgeIds;
  vEdgeIds.reserve(scnHexaEdgeCount);
  for(UINT k = 0; k < scnHexaEdgeCount; k++)
    vEdgeIds.push_back(k);

  int nLast = -1;
  size_t nEdgeCount = vEdgeIds.size();  // variable count of the edges to be processed.
  UINT nId, nNodeCount = 0, nAttemptsCount = 0;
  while((nEdgeCount != 0) && (nAttemptsCount < scnHexaEdgeCount))
  {
    for(UINT i = 0; i < nEdgeCount; i++)
    {
      nId = vEdgeIds.at(i);

      bool bExcl = nLast < 0 ? false :
                   nLast == scHexaEdges[nId].nExcl1 || nLast == scHexaEdges[nId].nExcl2 || nLast == scHexaEdges[nId].nExcl3 ||
                   nLast == scHexaEdges[nId].nExcl4 || nLast == scHexaEdges[nId].nExcl5;

      if(bExcl)
        continue;   // skip this edge, but do not erase its id from the vEdgeIds vector - it will be tried later.

      ray.orig = pElem->get_node(scHexaEdges[nId].n0)->pos;
      ray.dir = pElem->get_node(scHexaEdges[nId].n1)->pos - ray.orig;
      fEdgeLen = ray.dir.length();
      if(fEdgeLen < Const_Almost_Zero)
        return false;

      ray.dir /= fEdgeLen;
      if(!m_Plane.intersect(ray, vRes, fDist) || fDist > fEdgeLen)
      {
        vEdgeIds.erase(vEdgeIds.begin() + i); // erase the processed edge.
        break;
      }

      fKsi = fDist / fEdgeLen;  // I expect fKsi to be 0 <= fKsi <= 1.
      pNode = interpolate(vRes, pElem->get_node(scHexaEdges[nId].n0), pElem->get_node(scHexaEdges[nId].n1), fKsi);
  // Remember the last intersection edge index.
      nLast = nId;
      svHexaSection[nNodeCount] = pNode;
      nNodeCount++;

      vEdgeIds.erase(vEdgeIds.begin() + i); // erase the processed edge.
      break;
    }

    nEdgeCount = vEdgeIds.size();
    nAttemptsCount++;
  }

  if(nNodeCount < 3)
    return false;

  CFace* pFace = new CFace(svHexaSection[0], svHexaSection[1], svHexaSection[2]);
  vFaces.push_back(pFace);
  if(nNodeCount == 3)
    return true;

  pFace = new CFace(svHexaSection[0], svHexaSection[2], svHexaSection[3]);
  vFaces.push_back(pFace);
  if(nNodeCount == 4)
    return true;

  pFace = new CFace(svHexaSection[0], svHexaSection[3], svHexaSection[4]);
  vFaces.push_back(pFace);
  if(nNodeCount == 5)
    return true;

  pFace = new CFace(svHexaSection[0], svHexaSection[4], svHexaSection[5]);
  vFaces.push_back(pFace);
  return true;
}

CNode3D* CDomainCrossSection::interpolate(const Vector3F& vPos, CNode3D* p0, CNode3D* p1, float ksi)
{
  CNode3D* pNode = new CNode3D(vPos);
// Scalars:
  pNode->dens  = p0->dens + ksi * (p1->dens - p0->dens);
  pNode->press = p0->press + ksi * (p1->press - p0->press);
  pNode->temp  = p0->temp + ksi * (p1->temp - p0->temp);
  pNode->visc  = p0->visc + ksi * (p1->visc - p0->visc);
  pNode->cond  = p0->cond + ksi * (p1->cond - p0->cond);
  pNode->cp    = p0->cp + ksi * (p1->cp - p0->cp);
  pNode->phi   = p0->phi + ksi * (p1->phi - p0->phi);

// Vectors:
  pNode->vel   = p0->vel + ksi * (p1->vel - p0->vel);
  pNode->field = p0->field + ksi * (p1->field - p0->field);
  pNode->clmb = p0->clmb + ksi * (p1->clmb - p0->clmb);
  pNode->rf    = p0->rf + ksi * (p1->rf - p0->rf);

  m_vNodes.push_back(pNode);
  return pNode;
}

bool CDomainCrossSection::set_plane_type(int nType)
{
  if(m_nPlaneType != nType)
  {
    m_nPlaneType = nType;
    m_bReady = false;

    switch(m_nPlaneType)
    {
      case ptPlaneXY: m_Plane.norm = Vector3D(0, 0, 1); break;
      case ptPlaneXZ: m_Plane.norm = Vector3D(0, 1, 0); break;
      case ptPlaneYZ: m_Plane.norm = Vector3D(1, 0, 0); break;
    }

    return true;
  }

  return false;
}

const char* CDomainCrossSection::get_type_name(int nPlaneType) const
{
  switch(nPlaneType)
  {
    case ptPlaneXY: return _T("Plane XY");
    case ptPlaneXZ: return _T("Plane XZ");
    case ptPlaneYZ: return _T("Plane YZ");
  }

  return _T("None");
}

void CDomainCrossSection::save(CArchive& ar)
{
  UINT nVersion = 0;
  ar << nVersion;

  ar << m_Region.bEnabled;
  ar << m_nPlaneType;

  Vector3D vOrig = m_Plane.pos;
  ar << vOrig.x;
  ar << vOrig.y;
  ar << vOrig.z;

  Vector3D vNorm = m_Plane.norm;
  ar << vNorm.x;
  ar << vNorm.y;
  ar << vNorm.z;

  CString cName(m_Region.sName.c_str());
  ar << cName;
}

void CDomainCrossSection::load(CArchive& ar)
{
  UINT nVersion;
  ar >> nVersion;

  ar >> m_Region.bEnabled;
  ar >> m_nPlaneType;

  Vector3D vOrig;
  ar >> vOrig.x;
  ar >> vOrig.y;
  ar >> vOrig.z;
  m_Plane.pos = vOrig;

  Vector3D vNorm;
  ar >> vNorm.x;
  ar >> vNorm.y;
  ar >> vNorm.z;
  m_Plane.norm = vNorm;

  CString cName;
  ar >> cName;
  m_Region.sName = std::string(cName);

  m_bReady = false;
}

};  // namespace EvaporatingParticle.
