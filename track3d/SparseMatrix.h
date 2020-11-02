#pragma once

#ifndef _SparseMatrix_
#define _SparseMatrix_

#include "CObject.h"
#include <intsafe.h>

namespace EvaporatingParticle
{

struct CSparseMtxElem
{
  CSparseMtxElem(double val = 1, UINT col = 0)
    : fVal(val), nCol(col)
  {
  }

  double fVal;
  UINT   nCol;

  bool operator < (const CSparseMtxElem& elem) { return nCol < elem.nCol; }
};

typedef std::vector<double> CVectorND;
typedef std::vector<CSparseMtxElem> CVectorMtxElem;
typedef std::vector<CVectorMtxElem> CSparseMtxRows;
//---------------------------------------------------------------------------------------
// A sparse quadratic matrix of very large size. The elements are stored using the CSR 
// (compact sparse row) method.
//---------------------------------------------------------------------------------------
class CSparseMatrix : public CObject
{
public:
  CSparseMatrix(UINT nRowCount);

  UINT                        get_dimension() const;
  void                        set_dimension(UINT nRowCount);

  double                      get_elem(UINT nRow, UINT nCol) const;
  void                        set_elem(UINT nRow, UINT nCol, double fVal);

// A special case of setting elements for very huge matrices - by rows.
//  clear();
//  for(UINT nRow = 0; nRow < nRowCount; nRow++)
//    set_row(nRow, vRowElems);
  void                        clear();
// Call this function ONLY in a cycle from nRow = 0 to m_nRowCount - 1. Before the first call the clear() function must be called.
  void                        set_row(UINT nRow, const CVectorMtxElem& vElems);

  void                        get_column(UINT nCol, CVectorND& vRes) const;

// This matrix is multiplied from right with a vector of proper size: vRes = this * vSrc.
  bool                        multiply(const CVectorND& vSrc, CVectorND& vRes) const;

// This matrix is multiplied from right with another matrix of the same size: mRes = this * mSrc.
  bool                        multiply(const CSparseMatrix& mSrc, CSparseMatrix& mRes) const;

// Incomplete LU decomposition: this = L * U + R; L is a low-triangle matrix with unit diagonal; U is an upper-triangle matrix; R is a residual matrix.
  bool                        decomp_ILU(CSparseMatrix& mL, CSparseMatrix& mU);

// This function inverts only a sparse triangle U-matrix.
  static bool                 inverse_U(const CSparseMatrix& mU, CSparseMatrix& mInv);

// This function inverts only a sparse triangle L-matrix with unit main diagonal.
  static bool                 inverse_L(const CSparseMatrix& mL, CSparseMatrix& mInv);

protected:
  UINT                        find_elem(UINT nRow, UINT nCol) const;  // note: no range check is inside.

  void                        shrink_to_fit();

  void                        init();   // empty matrix.

private:
  UINT                        m_nRowCount;  // the number of equations in A*X = b.

  CSparseMtxRows              m_vElems;     // array of all the non-zero elements.

public:
// DEBUG static test functions:
  static void                 run_test(const std::string& cFileName);
  static void                 output_mtx(FILE* pStream, const CSparseMatrix& mtx, const char* pHeader);
  static void                 set_mtx_direct(CSparseMatrix& mtx);
  static void                 set_mtx_by_rows(CSparseMatrix& mtx);
// END DEBUG
};

//---------------------------------------------------------------------------------------
// Inline implementation
//---------------------------------------------------------------------------------------
inline UINT CSparseMatrix::get_dimension() const
{
  return m_nRowCount;
}

inline void CSparseMatrix::set_dimension(UINT nRowCount)
{
  m_nRowCount = nRowCount;
  init();
}

inline double CSparseMatrix::get_elem(UINT nRow, UINT nCol) const
{
  UINT nElem = find_elem(nRow, nCol);
  return nElem == UINT_MAX ? 0.0f : m_vElems[nRow].at(nElem).fVal;
}


};  // namespace EvaporatingParticle

#endif
