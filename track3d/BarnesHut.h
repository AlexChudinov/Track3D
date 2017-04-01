
#include "Vector3D.hpp"
#include "CalcThread.h"
#include "Elements.h"
#include <vector>

namespace EvaporatingParticle
{

struct ChargeData
{
  ChargeData(const Vector3D& p, double q)
    : pos(p), charge(q)
  {
  }

  Vector3D      pos;
  double        charge;
};

typedef std::vector<ChargeData*> CChargesColl;
//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
struct OctoTreeCell
{
  OctoTreeCell(const Vector3D& c, double s);
  ~OctoTreeCell();

  Vector3D      box_center;     // geometrical center of the cell.
  Vector3D      charge_center;  // center of charges of this cell.

  double        size;           // length of the cell's edge, cm.

  CChargesColl  charges;        // collection of ions populating this cell, (cm, cm, cm, CGSE).

  double        qSum;           // total charge (CGSE) contained in this cell.

  double        Qxx,
                Qyy,
                Qzz,
                Qxy,
                Qyz,
                Qxz;

  OctoTreeCell* pChild[8];

  UINT          ndep;           // recursion depth, the depth of the main tree is zero.
};

typedef std::vector<OctoTreeCell*> CellsColl;

const Vector3D vNull = Vector3D(0, 0, 0);
//---------------------------------------------------------------------------
// CBarnesHut - this class approximately computes Coulomb interaction 
//              between charged particles. 
//---------------------------------------------------------------------------
class CBarnesHut
{
public:
  CBarnesHut();
  ~CBarnesHut();

  enum
  {
    symNone   = 0,
    symXYonly = 1,
    symXZonly = 2,
    symBoth   = 3
  };

// User interface:
  double          get_dist_coeff() const;
  void            set_dist_coeff(double fCoeff);

  double          get_crit_radius() const;
  void            set_crit_radius(double fR);

  UINT            get_max_rec_depth() const;
  void            set_max_rec_depth(UINT nDepth);

  bool            get_enable_quad_terms() const;
  void            set_enable_quad_terms(bool bEnable);

// Public interface:
  double          coulomb_phi(const Vector3D& pos);     // Coulomb potential at a given point, CGS.
  Vector3D        coulomb_force(const Vector3D& pos);   // Coulomb electric field at a given point, CGS.

  void            create_main_cell(const Vector3D& cube_center, double cube_edge);

  void            add_particle(const Vector3D& pos, double charge); // both in CGS.

// In multi-threading applications this function must be called before start of the threads.
  void            prepare(CalcThreadVector& vThreads, CNodesCollection& vNodes, UINT nIter);

  Vector3D        get_center() const;
  OctoTreeCell*   get_main_cell() const;
  const CellsColl& get_cells() const;

  int             get_sym_type() const;
  void            set_sym_type(int nSymType);

  void            scale_all_charges(double fCoeff);

protected:
  double          coulomb_phi_cell(const Vector3D& vPos, OctoTreeCell* pCell);
// Correction for symmetry:
  double          correct_phi_for_symm(const Vector3D& vPos, OctoTreeCell* pCell);
  double          coulomb_phi_cell_symm(const Vector3D& vPos, OctoTreeCell* pCell, int nSymType);
// vRelPos = vPos - pCell->charge_center; fR5 = |vRelPos|^5 is used for quadrupole terms computation.
  double          quad_phi_cell(const Vector3D vRelPos, double fR5, OctoTreeCell* pCell);

  Vector3D        coulomb_force_cell(const Vector3D& vPos, OctoTreeCell* pCell);  // this is not force but E - electric field strength.
// Correction for symmetry:
  Vector3D        correct_field_for_symm(const Vector3D& vPos, OctoTreeCell* pCell);
  Vector3D        coulomb_field_cell_symm(const Vector3D& vPos, OctoTreeCell* pCell, int nSymType);
// vRelPos = vPos - pCell->charge_center; fRR2 = 1 / |vRelPos|^2, fRR3 = 1 / |vRelPos|^3; these two are used for quadrupole terms computation.
  Vector3D        quad_field_cell(const Vector3D vRelPos, double fRR2, double fRR3, OctoTreeCell* pCell);

// When the octo-tree is populated, compute moments with respect to the centers of charge of every cell.
  void            compute_moments(CalcThreadVector& vThreads);

  static UINT __stdcall build_moments(LPVOID pCalcThread);

// Creation of new cells.
// Get the child index by relative position of a new particle with respect to the parent box center.
  int             child_cell_index(const Vector3D& d);  // d = pos - parent_center.

  OctoTreeCell*   create_child_cell(int index, OctoTreeCell* pParentCell);
  
private:
// Particles distribution in cells, the upmost level.
  OctoTreeCell*   m_pTree;

  UINT            m_nMaxRecDepth;

  bool            m_bEnableQuadTerms;

  int             m_nSymmType;

// Let D be a size of a cell. If distance R between a point and the center of charge of the cell
// is greater than m_fDistCoeff * D, use approximate expression for the potential. Otherwise use q/R.
  double          m_fDistCoeff;   // dimensionless.

// If r = |vPos - pCell->charge_center| < m_fCritRadius the field grows lineary with r; this is the conception of a distributed charge. 
  double          m_fCritRadius;

// Run-time:
  bool            m_bReady;
  double          m_fOne_ovr_R3;  // 1./ m_fCritRadius^3.
  
// Linear (without tree-structure) collection of all cells. For fast computation of moments.
  CellsColl       m_vCells;
};

//---------------------------------------------------------------------------
// CBarnesHut inlines:
//---------------------------------------------------------------------------
inline double CBarnesHut::get_dist_coeff() const
{
  return m_fDistCoeff;
}

inline void CBarnesHut::set_dist_coeff(double fCoeff)
{
  m_fDistCoeff = fCoeff;
}

inline double CBarnesHut::get_crit_radius() const
{
  return m_fCritRadius;
}

inline UINT CBarnesHut::get_max_rec_depth() const
{
  return m_nMaxRecDepth;
}

inline void CBarnesHut::set_max_rec_depth(UINT nDepth)
{
  m_nMaxRecDepth = nDepth;
}

inline int CBarnesHut::get_sym_type() const
{
  return m_nSymmType;
}

inline void CBarnesHut::set_sym_type(int nSymType)
{
  m_nSymmType = nSymType;
}

inline bool CBarnesHut::get_enable_quad_terms() const
{
  return m_bEnableQuadTerms;
}

inline void CBarnesHut::set_enable_quad_terms(bool bEnable)
{
  m_bEnableQuadTerms = bEnable;
}

inline Vector3D CBarnesHut::get_center() const
{
  return m_pTree->box_center;
}

inline OctoTreeCell* CBarnesHut::get_main_cell() const
{
  return m_pTree;
}

inline const CellsColl& CBarnesHut::get_cells() const
{
  return m_vCells;
}

};  // namespace EvaporatingParticle
