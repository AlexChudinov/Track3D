
#pragma once
#ifndef _ELEMENTS_
#define _ELEMENTS_

#include <string>
#include <vector>
#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "constant.hpp"

namespace EvaporatingParticle
{

//-------------------------------------------------------------------------------------------------
// Intersection of an arbitrary ray with the scene.
//-------------------------------------------------------------------------------------------------
struct CRay
{
  CRay(const Vector3D& p, const Vector3D& d);

  Vector3D    orig,
              dir;
};

struct CRegFacePair
{
  CRegFacePair(UINT nr, UINT nf)
    : nReg(nr), nFace(nf)
  {
  }

  UINT  nReg, nFace;
};

struct CFace;
struct CNode3D;
struct CElem3D;
typedef std::vector<CFace*> CFacesCollection;
typedef std::vector<CElem3D*> CElementsCollection;
typedef std::vector<CNode3D*> CNodesCollection;
typedef std::vector<CRegFacePair> CFaceIndices;
//-------------------------------------------------------------------------------------------------
// CNode3D - a structure containing all the gas-dynamic parameters at the node of the input mesh.
//           In addition, it contains information on neighboring elements.
//-------------------------------------------------------------------------------------------------
struct CNode3D
{
  CNode3D()
    : dens(0), press(0), temp(0), visc(0), cond(0), cp(0), nInd(0)
  {
  }

  CNode3D(const Vector3D& p)
    : pos(p), dens(0), press(0), temp(0), visc(0), cond(0), cp(0), nInd(0)
  {
  }

  CNode3D(const Vector3D& p, const Vector3D& v, const Vector3D& dcE, const Vector3D& rfE, float den, float pres, float tmp, float vis, float con, float hcp)
    : pos(p), vel(v), field(dcE), rf(rfE), dens(den), press(pres), temp(tmp), visc(vis), cond(con), cp(hcp), nInd(0)
  {
  }

  size_t    nInd;   // index of the node in the global nodes collection.

  Vector3D  pos,
            vel,
            field,  // DC field.
            rf;     // RF field.

  float     dens,   // density.
            press,  // absolute pressure.
            temp,   // temperature, K.
            visc,   // dynamic viscosity.
            cond,   // thermal conductivity.
            cp;     // specific heat capacity at constant pressure.
// DEBUG
  float     phi;    // electric potential.
// END DEBUG

  void      set_data(float fPress, float fDens, float fTemp, float fVisc, float fCond, float fCp,
              const Vector3D& vVel, const Vector3D& vDCField, const Vector3D& vRFField);

// Neighbors:
  CElementsCollection vNbrElems;
  CFaceIndices        vNbrFaces;    // contains indices of boundary faces for boundary nodes; empty for inner nodes.
};

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
struct CBox
{
  CBox()
  {
  }

  CBox(const Vector3D& v0, const Vector3D& v1)
    : vMin(v0), vMax(v1)
  {
  }

  Vector3D  vMin,
            vMax;

  bool      inside(const Vector3D& p) const;
  Vector3D  get_center() const;

  bool      intersect(const CRay& ray) const;
  bool      intersect_plane(const CRay& ray, int nface) const;
};

//-------------------------------------------------------------------------------------------------
// CFace - a structure containing normalized vector "norm" in addition to pointers to
//         the three nodes of the triangle.
//-------------------------------------------------------------------------------------------------
struct CFace
{
  CFace(CNode3D* pn0, CNode3D* pn1, CNode3D* pn2);

  CNode3D*  p0;
  CNode3D*  p1;
  CNode3D*  p2;

  Vector3D  norm;
  CBox      box;    // bounding box is initialized in init().

  void      init();
  bool      intersect(const CRay& ray, double& dist) const;
};

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
struct CIntersectPoint
{
  CIntersectPoint(const Vector3D& v, double d) { pos = v; dist = d; } 

  Vector3D    pos;
  double      dist;

  bool operator < (const CIntersectPoint& v) const { return dist < v.dist; }
};

typedef std::vector<CIntersectPoint> CIntersectColl;

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
struct CRegion
{
  CRegion(const char* pName)
    : sName(pName), bEnabled(true), bSelected(false), bCrossSection(false)
  {
  }

  CBox              box;      // bounding box.
  CFacesCollection  vFaces;
  std::string       sName;

  bool              bEnabled;
  bool              bSelected;
  bool              bCrossSection;

  void              bounding_box();

// Use this function to find a single intersection with the region:
  bool              intersect(const CRay& ray, double& dist) const;

// Use this function to find all intersections with the region along a given ray:
  bool              intersect(const CRay& ray, CIntersectColl& results) const;
};

typedef std::vector<CRegion*> CRegionsCollection;
//-------------------------------------------------------------------------------------------------
// CTreeBox - is an element of the octo-tree of the CTracker.
//-------------------------------------------------------------------------------------------------
struct CTreeBox : public CBox
{
  CTreeBox()
    : CBox(), pChild(NULL), pParent(NULL), nLevel(0)
  {
  }

  CTreeBox(const Vector3D& v0, const Vector3D& v1)
    : CBox(v0, v1), pChild(NULL), pParent(NULL), nLevel(0)
  {
  }

  virtual ~CTreeBox();

  void                create_child_boxes();
  void                delete_child_boxes();

  void                collect_elems_inside();

  CTreeBox*           pParent;

// Either zero or array of 8 child boxes.
  CTreeBox**          pChild;

// If at least one node of an element is inside the box, this element must be collected here.
  CElementsCollection vElems;

  UINT                nLevel;
};

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
struct CPlane
{
  CPlane(const Vector3D& v, const Vector3D& n)
    : pos(v), norm(n)
  {
  }

  Vector3D  pos;
  Vector3D  norm;

  bool      inside(const Vector3D& p) const;
  bool      intersect(const CRay& ray, Vector3D& pnt, double& ksi) const;
};

typedef std::vector<CPlane> CPlanesSet;

//-------------------------------------------------------------------------------------------------
// CElem3D - a base abstract class for all spatial elements (tetrahedra, pyramides, wedges and hexes).
//-------------------------------------------------------------------------------------------------
struct CElem3D
{
  CElem3D() {};
  virtual ~CElem3D() {};

  size_t              nInd;   // index of this element in the global elements collection.

  CNodesCollection    vNodes;     // nodes.

  CPlanesSet          vFaces;     // faces of the element.

  CBox                box;        // bounding box:

  bool                valid;  // true if initialization is OK, false otherwise.

  void                add_plane(const Vector3D& p0, const Vector3D& p1, const Vector3D& p2);

  void                interpolate(const Vector3D& p, CNode3D& node) const;

  void                bounding_box();

  virtual bool        inside(const Vector3D& p) const;

// Purely virtual functions which must be overridden in the descendants:
  virtual bool        init() = 0;

// Input: parametric coordinates (s, t, u) and index of the shape function nfunc.
  virtual double      shape_func(double s, double t, double u, size_t nfunc) const = 0;

// Input: parametric coordinates (s, t, u), index of the shape function nfunc and index of the variable.
  virtual double      shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const = 0;

// Newton method support.
// Discrepancy of positions: computed for a k-th iteration and given point.
  virtual Vector3D    func(double s, double t, double u, const Vector3D& pos) const;

//                            dFx/ds, dFx/dt, dFx/du
// Returned matrix contains:  dFy/ds, dFy/dt, dFy/du
//                            dFz/ds, dFz/dt, dFz/du
//
  virtual Matrix3D    deriv(double s, double t, double u) const;

// Parametric coordinates for a point inside the element. If not overridden this function uses
// the Newton iteration method to obtain parametric coordinates by the position inside the element.
  virtual bool        param(const Vector3D& p, double& s, double& t, double& u) const;
};

//-------------------------------------------------------------------------------------------------
//                    CTetra - a tetrahedron object of 4 nodes.
//-------------------------------------------------------------------------------------------------
struct CTetra : public CElem3D
{
  CTetra(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3);

  virtual bool        init();
  virtual double      shape_func(double s, double t, double u, size_t node_id) const;
  virtual double      shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const;
};

//-------------------------------------------------------------------------------------------------
// CPyramid - a pyramid object of 5 nodes. To determine whether or not a point is inside the element,
//            the pyramid is subdivided in two tetras.
//-------------------------------------------------------------------------------------------------
struct CPyramid : public CElem3D
{
  CPyramid(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4);
  virtual ~CPyramid();

  virtual bool      init();
  virtual bool      inside(const Vector3D& p) const;
  virtual double    shape_func(double s, double t, double u, size_t node_id) const;
  virtual double    shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const;

  CTetra*           pTet1;
  CTetra*           pTet2;
};

typedef std::vector<CTetra*> CTetraColl;
//-------------------------------------------------------------------------------------------------
// CWedge - has 6 nodes. To determine whether or not a point is inside the element, the wedge is
//          subdivided in eight tetras with a common node at the center of the wedge.
//-------------------------------------------------------------------------------------------------
struct CWedge : public CElem3D
{
  CWedge(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4, CNode3D* p5);
  virtual ~CWedge();

  virtual bool      init();
  virtual bool      inside(const Vector3D& p) const;
  virtual double    shape_func(double s, double t, double u, size_t node_id) const;
  virtual double    shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const;

  CTetraColl        vTetras;
  CNode3D*          pCenter;
};

//-------------------------------------------------------------------------------------------------
// CHexa - has 8 nodes.
//-------------------------------------------------------------------------------------------------
struct CHexa : public CElem3D
{
  CHexa(CNode3D* p0, CNode3D* p1, CNode3D* p2, CNode3D* p3, CNode3D* p4, CNode3D* p5, CNode3D* p6, CNode3D* p7);

  virtual bool      init();
  virtual bool      inside(const Vector3D& p) const;
  virtual double    shape_func(double s, double t, double u, size_t node_id) const;
  virtual double    shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const;

  CTetraColl        vTetras;
  CNode3D*          pCenter;
};

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
struct CEdge
{
  CEdge(CNode3D* p0, CNode3D* p1, CFace* pf0)
    : pNode0(p0), pNode1(p1), pFace0(pf0), pFace1(NULL)
  {
  }

  CNode3D*  pNode0;
  CNode3D*  pNode1;
// Neighboring faces. If one of them is NULL the edge is a boundary one.
  CFace*    pFace0;
  CFace*    pFace1;

  bool operator == (const CEdge& edge)
  { return (pNode0 == edge.pNode0) && (pNode1 == edge.pNode1) || (pNode0 == edge.pNode1) && (pNode1 == edge.pNode0); }
};

typedef std::vector<CEdge> CEdgesVector;

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
class CSimpleStack
{
public:
  CSimpleStack(UINT nSize);
  ~CSimpleStack();

  void    add(double fVal);
  double  get(UINT nOffset = 0) const;
  double  average(UINT nLen) const;

protected:
  void    init();

private:
  double* m_pArr;
  UINT    m_nSize,
          m_nPos;
};

//-------------------------------------------------------------------------------------------------
// CAverValue
//-------------------------------------------------------------------------------------------------
struct CAverValue
{
  CAverValue(double x = 0, double y1 = 0, double y2 = 0)
    : fX(x), fY1(y1), fY2(y2)
  {
  }

  double    fX,
            fY1,
            fY2;
};

enum CAverType
{
  atAverage = 0,
  atMax     = 1
};

//-------------------------------------------------------------------------------------------------
// CAverBin
//-------------------------------------------------------------------------------------------------
class CAverBin
{
public:
  CAverBin()
    : m_fArgMin(0), m_fArgMax(1), m_nCount(0)
  {
  }

  bool        add_val(const CAverValue& vVal);
  void        set_bounds(double fMin, double fMax);
  void        set_type(CAverType nType);

  CAverValue  get() const;
  UINT        get_count() const;

private:
  double      m_fArgMin,
              m_fArgMax;

  UINT        m_nCount;   // count of values that are currently in the bin.

// Type of processing CAverValue::fY1 and CAverValue::fY2. CAverValue::fX is always averaged.
  CAverType   m_nAverType;

  CAverValue  m_vAver;    // averaged values stored in the bin.
};

typedef std::vector<CAverBin> CBinVector;
//-------------------------------------------------------------------------------------------------
// CAveragingEngine
//-------------------------------------------------------------------------------------------------
class CAveragingEngine
{
public:
  CAveragingEngine();
  ~CAveragingEngine();

  void        clear();

  void        set_type(CAverType nType);

  void        get_range(double& fArgMin, double& fArgMax);
  void        set_range(double fArgMin, double fArgMax, UINT nBinsCount);

  void        add_value(const CAverValue& vVal);

  UINT        get_count() const;

  CAverBin    get_aver_bin(UINT nPos) const;

private:
  void        init(UINT nBinsCount);

  double      m_fArgMin,
              m_fArgMax,
              m_fStep;

  CAverType   m_nAverType;
  CBinVector  m_vAverBins;
};


//-------------------------------------------------------------------------------------------------
// CBox: inline implementation.
//-------------------------------------------------------------------------------------------------
inline bool CBox::inside(const Vector3D& p) const
{
  return !(p.x < vMin.x || p.x > vMax.x || p.y < vMin.y || p.y > vMax.y || p.z < vMin.z || p.z > vMax.z);
}

inline Vector3D CBox::get_center() const
{
  return 0.5 * (vMin + vMax);
}

//-------------------------------------------------------------------------------------------------
// CPlane: inline implementation.
//-------------------------------------------------------------------------------------------------
inline bool CPlane::inside(const Vector3D& p) const
{
  return ((p - pos) & norm) <= 0;
}

//-------------------------------------------------------------------------------------------------
// CNode3D: inline implementation.
//-------------------------------------------------------------------------------------------------
inline void CNode3D::set_data(float fPress, float fDens, float fTemp, float fVisc, float fCond, float fCp,
              const Vector3D& vVel, const Vector3D& vDCField, const Vector3D& vRFField)
{
  press = fPress;
  dens = fDens;
  temp = fTemp;
  visc = fVisc;
  cond = fCond;
  cp = fCp;

  vel = vVel;
  field = vDCField;
  rf = vRFField;
}

//-------------------------------------------------------------------------------------------------
// CAverBin: inline implementation.
//-------------------------------------------------------------------------------------------------
inline void CAverBin::set_bounds(double fMin, double fMax)
{
  m_fArgMin = fMin;
  m_fArgMax = fMax;
}

inline void CAverBin::set_type(CAverType nType)
{
  m_nAverType = nType;
}

inline CAverValue CAverBin::get() const
{
  return m_vAver;
}

inline UINT CAverBin::get_count() const
{
  return m_nCount;
}

//-------------------------------------------------------------------------------------------------
// CAveragingEngine: inline implementation.
//-------------------------------------------------------------------------------------------------
inline void CAveragingEngine::set_type(CAverType nType)
{
  m_nAverType = nType;
}

inline UINT CAveragingEngine::get_count() const
{
  return m_vAverBins.size();
}

inline CAverBin CAveragingEngine::get_aver_bin(UINT nPos) const
{
  return m_vAverBins.at(nPos);
}

inline void CAveragingEngine::get_range(double& fArgMin, double& fArgMax)
{
  fArgMin = m_fArgMin;
  fArgMax = m_fArgMax;
}

typedef std::vector<CElem3D*> Elems;
typedef std::vector<CNode3D*> Nodes;

void sortElems(Elems& elems);

};  // namespace EvaporatingParticle

#endif