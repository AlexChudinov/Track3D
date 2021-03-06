
#pragma once
#ifndef _ELEMENTS_
#define _ELEMENTS_

#include <string>
#include <vector>
#include "vector3d.hpp"
#include "matrix3d.hpp"
#include "constant.hpp"

//[AC 24/03/2017]
//Memory manager
#include <../utilities/MemoryPool.h>
//[/AC]

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
  CRegFacePair(UINT nr = UINT_MAX, UINT nf = UINT_MAX)
    : nReg(nr), nFace(nf)
  {
  }

  UINT  nReg, nFace;

  bool operator == (const CRegFacePair& face) const
  { 
	  return (nReg == face.nReg) && (nFace == face.nFace); 
  }

  bool operator < (const CRegFacePair& face) const
  {
	  return nReg < face.nReg || (nReg == face.nReg && nFace < face.nFace);
  }
};

struct CFace;
struct CNode3D;
struct CElem3D;
typedef std::vector<CFace*> CFacesCollection;
typedef std::vector<CElem3D*> CElementsCollection;

typedef std::vector<CNode3D*> CNodesCollection;
typedef std::vector<CNode3D> CNodesVector;

typedef std::vector<CRegFacePair> CFaceIndices;
typedef std::vector<UINT> CIndexVector;
//-------------------------------------------------------------------------------------------------
// CNode3D - a structure containing all the gas-dynamic parameters at the node of the input mesh.
//           In addition, it contains information on neighboring elements.
//-------------------------------------------------------------------------------------------------
struct CNode3D : public BlockAllocator<CNode3D>
{
  CNode3D()
    : nInd(0), dens(0), press(0), temp(0), visc(0), cond(0), cp(0), phi(0)
  {
  }

  CNode3D(const Vector3F& p)
    : nInd(0), pos(p), dens(0), press(0), temp(0), visc(0), cond(0), cp(0), phi(0)
  {
  }

  CNode3D(const Vector3F& p, const Vector3F& v, const Vector3F& dcE, const Vector3F& rfE, float den, float pres, float tmp, float vis, float con, float hcp)
    : nInd(0), pos(p), vel(v), field(dcE), rf(rfE), dens(den), press(pres), temp(tmp), visc(vis), cond(con), cp(hcp), phi(0)
  {
  }

  size_t    nInd;   // index of the node in the global nodes collection.

  Vector3F  pos,
            vel,
            field,  // DC field.
            clmb,   // Coulomb field.
            rf;     // RF field.

  float     dens,   // density.
            press,  // absolute pressure.
            temp,   // temperature, K.
            visc,   // dynamic viscosity.
            cond,   // thermal conductivity.
            cp;     // specific heat capacity at constant pressure.

  float     phi;    // electric potential.

  void      set_data(float fPress, float fDens, float fTemp, float fVisc, float fCond, float fCp,
              const Vector3F& vVel, const Vector3F& vDCField, const Vector3F& vRFField, const Vector3F& vClmb = Vector3F(0, 0, 0));

// Neighbors:
  CIndexVector        vNbrElems;
  CIndexVector        vNbrNodes;    // for Dirichlet cell creation around this node.
  CFaceIndices        vNbrFaces;    // contains indices of boundary faces for boundary nodes; empty for inner nodes.

  void                shrink_to_fit();

  CElementsCollection get_nbr_elems() const;

  //[AC 03/03/2017]
  //Returns references to a neighbor element indices array
  const CIndexVector& nbr_elems() const;
  //[/AC 03/03/2017]

// Returns true if: 1). The node is a boundary node; 2). The node belongs to at least two regions.
  bool                is_wireframe_node() const;
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
  bool      intersect_plane(const CRay& ray, int nface, Vector3D& vRes) const;
};

//-------------------------------------------------------------------------------------------------
// CFace - a structure containing normalized vector "norm" in addition to pointers to
//         the three nodes of the triangle.
//-------------------------------------------------------------------------------------------------
struct CFace : public BlockAllocator<CFace>
{
  CFace(CNode3D* pn0, CNode3D* pn1, CNode3D* pn2);

  CNode3D*  p0;
  CNode3D*  p1;
  CNode3D*  p2;

  Vector3D  norm;
  CBox      box;    // bounding box is initialized in init().

  void      init();
  bool      intersect(const CRay& ray, double& dist) const;
  double    square() const;

  bool      is_wireframe_face() const;  // returns true if at least one node of this face belongs to wireframe.
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
struct CRegion : public BlockAllocator<CRegion>
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
  bool              intersect(const CRay& ray, double& fDist, UINT& nFaceID) const;

// Use this function to find all intersections with the region along a given ray:
  bool              intersect(const CRay& ray, CIntersectColl& results) const;
};

typedef std::vector<CRegion*> CRegionsCollection;
//-------------------------------------------------------------------------------------------------
// CTreeBox - is an element of the octo-tree of the CTracker.
//-------------------------------------------------------------------------------------------------
struct CTreeBox : public CBox, public BlockAllocator<CTreeBox>
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
struct CPlane : public BlockAllocator<CPlane>
{
  CPlane(const Vector3D& v = Vector3D(0, 0, 0), const Vector3D& n = Vector3D(0, 0, 1))
    : pos(v), norm(n)
  {
  }

  Vector3D  pos;
  Vector3D  norm;

  bool      inside(const Vector3D& p) const;
  bool      intersect(const CRay& ray, Vector3D& pnt, double& ksi) const;

  bool      operator == (const CPlane& plane) const;
};

typedef std::vector<CPlane> CPlanesSet;

//-------------------------------------------------------------------------------------------------
// CElem3D - a base abstract class for all spatial elements (tetrahedra, pyramides, wedges and hexes).
//-------------------------------------------------------------------------------------------------
struct CElem3D
{
  size_t                  nInd;       // index of this element in the global elements collection.

  CBox                    box;        // bounding box:

  bool                    valid;  // true if initialization is OK, false otherwise.

  bool                    coeff(const Vector3D& p, double* pWeight) const;  // pWeight must be allocated in the calling function.
  void                    interpolate(const Vector3D& p, CNode3D& node) const;  // obsolete, call coeff(...) instead.

  void                    bounding_box();

// Purely virtual functions which must be overridden in the descendants:
  virtual bool            inside(const Vector3D& p) const = 0;

  virtual bool            init() = 0;

// Input: parametric coordinates (s, t, u) and index of the shape function nfunc.
  virtual double          shape_func(double s, double t, double u, size_t nfunc) const = 0;

// Input: parametric coordinates (s, t, u), index of the shape function nfunc and index of the variable.
  virtual double          shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const = 0;

// Newton method support.
// Initial guess for parametric coordinates (the middle of the element). Note that at this point, which is always inside the element,
// the shape functions of each node return 1/N, where N is the nodes count of the element.
  virtual void            elem_middle(double& s, double& t, double& u) const = 0;

// Discrepancy of positions: computed for a k-th iteration and given point.
  virtual Vector3D        func(double s, double t, double u, const Vector3D& pos) const;

//                            dFx/ds, dFx/dt, dFx/du
// Returned matrix contains:  dFy/ds, dFy/dt, dFy/du
//                            dFz/ds, dFz/dt, dFz/du
//
  virtual Matrix3D        deriv(double s, double t, double u) const;

// Parametric coordinates for a point inside the element. If not overridden this function uses
// the Newton iteration method to obtain parametric coordinates by the position inside the element.
  virtual bool            param(const Vector3D& p, double& s, double& t, double& u) const;

  virtual UINT            get_node_count() const = 0;
  virtual UINT            get_node_index(UINT nInd) const = 0;  // index of this node in the global nodes collection.
  virtual CNode3D*        get_node(UINT nInd) const = 0;

  CNodesCollection        get_nodes() const;

  //[AC 03/03/2017]
  //Returns pointer to element node c-like array
  virtual const UINT*     nodes() const = 0;
  //[/AC]

  //[AC 27/03/2017] memory manager
  static void operator delete(void* ptr, size_t n);
  virtual void            deleteObj() = 0;
  //[/AC]
};

//-------------------------------------------------------------------------------------------------
//    CElemTetra - a simple object that consists of 4 tetrahedron planes (only planes, no nodes).
//-------------------------------------------------------------------------------------------------
struct CElemTetra
{
  CPlane            vFaces[4];

  bool              contains(const Vector3D& p) const;

  void              elem_tetra(const Vector3D& p0, const Vector3D& p1, const Vector3D& p2, const Vector3D& p3);
  void              add_plane(UINT nInd, const Vector3D& p0, const Vector3D& p1, const Vector3D& p2);
};

//-------------------------------------------------------------------------------------------------
//                    CTetra - a tetrahedron object of 4 nodes.
//-------------------------------------------------------------------------------------------------
struct CTetra : public CElem3D, public BlockAllocator<CTetra>
{
	//[AC] memory manager
	using CElem3D::operator delete;
	//[/AC]
  CTetra(UINT i0, UINT i1, UINT i2, UINT i3);
  CTetra(const Vector3D& v0, const Vector3D& v1, const Vector3D& v2, const Vector3D& v3);

  UINT              vTetNodeIds[4]; // indices of the tetra nodes in the global collection.
  CElemTetra        vFaces;         // faces of the element.

  virtual bool      init();
  virtual bool      inside(const Vector3D& p) const;
  virtual double    shape_func(double s, double t, double u, size_t node_id) const;
  virtual double    shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const;

  virtual void      elem_middle(double& s, double& t, double& u) const { s = 0.25; t = 0.25; u = 0.25; }

  virtual UINT      get_node_count() const { return 4; }
  virtual UINT      get_node_index(UINT nInd) const { return vTetNodeIds[nInd]; }  // index of this node in the global nodes collection.
  virtual CNode3D*  get_node(UINT nInd) const;

  void              add_plane(UINT nInd, const Vector3D& p0, const Vector3D& p1, const Vector3D& p2);

  //[AC 03/03/2017]
  //Returns pointer to element node c-like array
  const UINT* nodes() const;
  //[/AC]

  //[AC 27/03/2017] memory manager
  virtual void deleteObj();
  //[/AC]
};

//-------------------------------------------------------------------------------------------------
// CPyramid - a pyramid object of 5 nodes. To determine whether or not a point is inside the element,
//            the pyramid is subdivided in two tetras.
//-------------------------------------------------------------------------------------------------
struct CPyramid : public CElem3D, public BlockAllocator<CPyramid>
{
	//[AC] memory manager
	using CElem3D::operator delete;
	//[/AC]
  CPyramid(UINT i0, UINT i1, UINT i2, UINT i3, UINT i4);
  virtual ~CPyramid();

  UINT              vPyrNodeIds[5]; // indices of the pyramid nodes in the global collection.

  virtual bool      init();
  virtual bool      inside(const Vector3D& p) const;
  virtual double    shape_func(double s, double t, double u, size_t node_id) const;
  virtual double    shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const;

  virtual void      elem_middle(double& s, double& t, double& u) const { s = 0.5; t = 0.5; u = 0.2; }

  virtual UINT      get_node_count() const { return 5; }
  virtual UINT      get_node_index(UINT nInd) const { return vPyrNodeIds[nInd]; }  // index of this node in the global nodes collection.
  virtual CNode3D*  get_node(UINT nInd) const;

  CElemTetra        vSubTetra[2];

  //[AC 03/03/2017]
  //Returns pointer to element node c-like array
  const UINT* nodes() const;
  //[/AC]

  //[AC 27/03/2017] memory manager
  virtual void deleteObj();
  //[/AC]
};

//-------------------------------------------------------------------------------------------------
// CWedge - has 6 nodes. To determine whether or not a point is inside the element, the wedge is
//          subdivided in eight tetras with a common node at the center of the wedge.
//-------------------------------------------------------------------------------------------------
struct CWedge : public CElem3D, public BlockAllocator<CWedge>
{
	//[AC] memory manager
	using CElem3D::operator delete;
	//[/AC]
  CWedge(UINT i0, UINT i1, UINT i2, UINT i3, UINT i4, UINT i5);
  virtual ~CWedge();

  UINT              vWdgNodeIds[6]; // indices of the wedge nodes in the global collection.

  virtual bool      init();
  virtual bool      inside(const Vector3D& p) const;
  virtual double    shape_func(double s, double t, double u, size_t node_id) const;
  virtual double    shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const;

  virtual void      elem_middle(double& s, double& t, double& u) const { s = Const_One_Third; t = Const_One_Third; u = 0.5; }

  virtual UINT      get_node_count() const { return 6; }
  virtual UINT      get_node_index(UINT nInd) const { return vWdgNodeIds[nInd]; }  // index of this node in the global nodes collection.
  virtual CNode3D*  get_node(UINT nInd) const;

  CElemTetra        vSubTetra[8];

  //[AC 03/03/2017]
  //Returns pointer to element node c-like array
  const UINT* nodes() const;
  //[/AC]

  //[AC 27/03/2017] memory manager
  virtual void deleteObj();
  //[/AC]
};

//-------------------------------------------------------------------------------------------------
// CHexa - has 8 nodes.
//-------------------------------------------------------------------------------------------------
struct CHexa : public CElem3D, public BlockAllocator<CHexa>
{
	//[AC] memory manager
	using CElem3D::operator delete;
	//[/AC]
  CHexa(UINT i0, UINT i1, UINT i2, UINT i3, UINT i4, UINT i5, UINT i6, UINT i7);

  UINT              vHexNodeIds[8]; // indices of the wedge nodes in the global collection.

  virtual bool      init();
  virtual bool      inside(const Vector3D& p) const;
  virtual double    shape_func(double s, double t, double u, size_t node_id) const;
  virtual double    shape_func_deriv(double s, double t, double u, size_t nfunc, size_t nvar) const;

  virtual void      elem_middle(double& s, double& t, double& u) const { s = 0.5; t = 0.5; u = 0.5; }

  virtual UINT      get_node_count() const { return 8; }
  virtual UINT      get_node_index(UINT nInd) const { return vHexNodeIds[nInd]; }  // index of this node in the global nodes collection.
  virtual CNode3D*  get_node(UINT nInd) const;

  CElemTetra        vSubTetra[12];

  //[AC 03/03/2017]
  //Returns pointer to element node c-like array
  const UINT* nodes() const;
  //[/AC]

  //[AC 27/03/2017] memory manager
  virtual void deleteObj();
  //[/AC]
};

//-------------------------------------------------------------------------------------------------
//
//-------------------------------------------------------------------------------------------------
struct CEdge : public BlockAllocator<CEdge>
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

inline double CFace::square() const
{
  return 0.5 * ((p1->pos - p0->pos) * (p2->pos - p0->pos)).length();
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
              const Vector3F& vVel, const Vector3F& vDCField, const Vector3F& vRFField, const Vector3F& vClmbField)
{
  press = fPress;
  dens = fDens;
  temp = fTemp;
  visc = fVisc;
  cond = fCond;
  cp = fCp;

  vel = vVel;
  field = vDCField;
  clmb = vClmbField;
  rf = vRFField;
}

inline void CNode3D::shrink_to_fit()
{
  vNbrElems.shrink_to_fit();
  vNbrFaces.shrink_to_fit();
}

//-------------------------------------------------------------------------------------------------
// CElem3D: inline implementation.
//-------------------------------------------------------------------------------------------------
/*
inline CNodesCollection CElem3D::get_nodes() const
{
  CNodesCollection vNodes(get_node_count());
  for(UINT i = 0; i < get_node_count(); i++)
    vNodes[i] = vGlobNodes.at(get_node_index(i));

  return vNodes;
}
*/
//[AC 27/03/2017] memory manager
inline void CElem3D::operator delete(void * ptr, size_t n)
{
	((CElem3D*)ptr)->deleteObj();
}
//[/AC]

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

};  // namespace EvaporatingParticle

#endif