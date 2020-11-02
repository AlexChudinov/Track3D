
#include "stdafx.h"
#include "SparseMatrix.h"
#include "constant.hpp"
#include <algorithm>

namespace EvaporatingParticle
{

CSparseMatrix::CSparseMatrix(UINT nRowCount)
  : m_nRowCount(nRowCount)
{
  init();
}

void CSparseMatrix::init()
{
  m_vElems.resize(m_nRowCount);
}

void CSparseMatrix::clear()
{
  size_t nRowCount = m_vElems.size();
  for(UINT nRow = 0; nRow < nRowCount; nRow++)
    m_vElems[nRow].clear();

  m_vElems.clear();
}

void CSparseMatrix::set_row(UINT nRow, const CVectorMtxElem& vElems)
{
  UINT nRowSize = vElems.size();
  m_vElems[nRow].reserve(nRowSize);
  for(UINT i = 0; i < nRowSize; i++)
    m_vElems[nRow].push_back(vElems.at(i));
}

UINT CSparseMatrix::find_elem(UINT nRow, UINT nCol) const
{
  const CVectorMtxElem& vRowElems = m_vElems[nRow];
  UINT nRowSize = vRowElems.size();
  for(UINT j = 0; j < nRowSize; j++)
  {
    if(vRowElems[j].nCol == nCol)
      return j;
  }

  return UINT_MAX;
}

void CSparseMatrix::set_elem(UINT nRow, UINT nCol, double fVal)
{
  CVectorMtxElem& vRowElems = m_vElems[nRow];

  UINT nElem = find_elem(nRow, nCol);
  if(nElem != UINT_MAX)
  {
    vRowElems[nElem].fVal = fVal;
    return;
  }
  
// New element to be inserted:
  CSparseMtxElem elem(fVal, nCol);
  vRowElems.push_back(elem);
}

void CSparseMatrix::get_column(UINT nCol, CVectorND& vRes) const
{
  vRes.resize(m_nRowCount, 0);

  for(UINT nRow = 0; nRow < m_nRowCount; nRow++)
    vRes[nRow] = get_elem(nRow, nCol);
}

bool CSparseMatrix::multiply(const CVectorND& vSrc, CVectorND& vRes) const
{
  if(vSrc.size() != m_nRowCount)
    return false;

  vRes.resize(m_nRowCount);

  for(UINT nRow = 0; nRow < m_nRowCount; nRow++)
  {
    vRes[nRow] = 0;
    const CVectorMtxElem& vRowElems = m_vElems[nRow];
    UINT nRowSize = vRowElems.size();
    for(UINT i = 0; i < nRowSize; i++)
    {
      vRes[nRow] += vRowElems[i].fVal * vSrc[vRowElems[i].nCol];
    }
  }

  return true;
}

bool CSparseMatrix::multiply(const CSparseMatrix& mSrc, CSparseMatrix& mRes) const
{
  if(mSrc.get_dimension() != m_nRowCount)
    return false;

  if(mRes.get_dimension() != m_nRowCount)
    mRes.set_dimension(m_nRowCount);

  double fAij;
  CVectorND vCol(m_nRowCount, 0), vRes(m_nRowCount, 0);
  for(UINT nCol = 0; nCol < m_nRowCount; nCol++)
  {
    mSrc.get_column(nCol, vCol);
    multiply(vCol, vRes);
    for(UINT nRow = 0; nRow < m_nRowCount; nRow++)
    {
      fAij = vRes[nRow];
      if(fabs(fAij) < Const_Almost_Zero)
        continue;

      mRes.set_elem(nRow, nCol, fAij);
    }
  }

  return true;
}

bool CSparseMatrix::decomp_ILU(CSparseMatrix& mL, CSparseMatrix& mU)
{
  mL.set_dimension(m_nRowCount);
  mU.set_dimension(m_nRowCount);

  set_job_name("System matrix ILU decomposition");

  UINT nElem, nElemL, nElemU;
  double fAkj, fLki, fUij, fSum, fUjj;
  for(UINT k = 0; k < m_nRowCount; k++)
  {
    if(k % 100 == 0)
      set_progress(100 * double(k) / m_nRowCount);

    if(get_terminate_flag())
      return false;

    const CVectorMtxElem& vRowElems = m_vElems[k];

    for(UINT j = 0; j < k; j++)
    {
      nElem = find_elem(k, j);
      if(nElem != UINT_MAX)
      {
        fAkj = vRowElems[nElem].fVal;
        fSum = 0;
        for(UINT i = 0; i < j; i++)
        {
          nElemL = mL.find_elem(k, i);
          if(nElemL == UINT_MAX)
            continue;

          fLki = mL.m_vElems[k].at(nElemL).fVal;

          nElemU = mU.find_elem(i, j);
          if(nElemU == UINT_MAX)
            continue;

          fUij = mU.m_vElems[i].at(nElemU).fVal;

          fSum += fLki * fUij;
        }

        fUjj = mU.get_elem(j, j);
        if(fabs(fUjj) < Const_Almost_Zero)
        {
          AfxMessageBox("decomp_ILU: bad input matrix.");
          return false;
        }

        mL.set_elem(k, j, (fAkj - fSum) / fUjj);
      }
    }

    mL.set_elem(k, k, 1);

    for(UINT j = k; j < m_nRowCount; j++)
    {
      nElem = find_elem(k, j);
      if(nElem != UINT_MAX)
      {
        fAkj = vRowElems[nElem].fVal;
        fSum = 0;
        for(UINT i = 0; i < k; i++)
        {
          nElemL = mL.find_elem(k, i);
          if(nElemL == UINT_MAX)
            continue;

          fLki = mL.m_vElems[k].at(nElemL).fVal;

          nElemU = mU.find_elem(i, j);
          if(nElemU == UINT_MAX)
            continue;

          fUij = mU.m_vElems[i].at(nElemU).fVal;

          fSum += fLki * fUij;
        }

        mU.set_elem(k, j, fAkj - fSum);
      }
    }
  }

  return true;
}

static void bad_diag_msg(UINT i)
{
  char buff[16];
  CString sMsg = CString("inverse_U: small value on the input matrix diagonal, i = ") + CString(itoa(i, buff, 10));
  AfxMessageBox((const char*)sMsg);
}

bool CSparseMatrix::inverse_L(const CSparseMatrix& mL, CSparseMatrix& mInv)
{
  UINT N = mL.get_dimension();
  if(mInv.get_dimension() != N)
    mInv.set_dimension(N);

  UINT nElem;
  double fAik, fBij, fBkj, fSum;
  for(UINT i = 0; i < N; i++)
  {
    mInv.set_elem(i, i, 1);

    const CVectorMtxElem& vRowElems = mL.m_vElems[i];

    for(UINT j = 0; j < i; j++)
    {
      fSum = 0;
      for(UINT k = j; k < i; k++)
      {
        nElem = mL.find_elem(i, k);
        if(nElem == UINT_MAX)
          continue;

        fAik = vRowElems[nElem].fVal;

        nElem = mInv.find_elem(k, j);
        if(nElem == UINT_MAX)
          continue;

        fBkj = mInv.m_vElems[k].at(nElem).fVal;

        fSum += fAik * fBkj;
      }

      fBij = -fSum;
      if(fabs(fBij) < Const_Almost_Zero)
        continue;

      mInv.set_elem(i, j, fBij);
    }
  }

  mInv.shrink_to_fit();
  return true;
}

bool CSparseMatrix::inverse_U(const CSparseMatrix& mU, CSparseMatrix& mInv)
{
  UINT N = mU.get_dimension();
  if(mInv.get_dimension() != N)
    mInv.set_dimension(N);

  UINT nElem;
  double fAii, fAik, fBij, fBkj, fSum;
  for(int j = 0; j < N; j++)
  {
    for(int i = j; i >= 0; i--)
    {
      if(i == j)
      {
        fAii = mU.get_elem(i, i);
        if(fabs(fAii) < Const_Almost_Zero)
        {
          bad_diag_msg(i);
          return false;
        }

        mInv.set_elem(i, i, 1. / fAii);
      }
      else
      {
        fSum = 0;
        for(UINT k = i + 1; k <= j; k++)
        {
          nElem = mU.find_elem(i, k);
          if(nElem == UINT_MAX)
            continue;

          fAik = mU.m_vElems[i].at(nElem).fVal;

          nElem = mInv.find_elem(k, j);
          if(nElem == UINT_MAX)
            continue;

          fBkj = mInv.m_vElems[k].at(nElem).fVal;

          fSum += fAik * fBkj;
        }

        fBij = -fSum / mU.get_elem(i, i);
        if(fabs(fBij) < Const_Almost_Zero)
          continue;

        mInv.set_elem(i, j, fBij);
      }
    }
  }

  mInv.shrink_to_fit();
  return true;
}

void CSparseMatrix::shrink_to_fit()
{
  for(UINT nRow = 0; nRow < m_nRowCount; nRow++)
    m_vElems[nRow].shrink_to_fit();
}

// DEBUG
void CSparseMatrix::run_test(const std::string& cFileName)
{
  CSparseMatrix mtx(7);
//  set_mtx_direct(mtx);
  set_mtx_by_rows(mtx);

  FILE* pStream;
  errno_t nErr = fopen_s(&pStream, cFileName.c_str(), (const char*)("w"));
  if((nErr != 0) || (pStream == NULL))
    return;

// Matrix itself:
  output_mtx(pStream, mtx, "Matrix itself:");

  int nRes;
  double fVal;
  UINT N = mtx.m_nRowCount;

// Matrix by vector multiplication:
  CVectorND vVec = { 1, 2, 3, 4, 5, 6, 7 }, vRes(N);
  mtx.multiply(vVec, vRes);
  fputs("Matrix by vector multiplication:\n", pStream);
  for(UINT i = 0; i < N; i++)
    nRes = i < N - 1 ? fprintf(pStream, "%10.5lf", vRes[i]) : fprintf(pStream, "%10.5lf\n\n", vRes[i]); 

  mtx.get_column(3, vRes);
  fputs("Column #3 (beginning from 0):\n", pStream);
  for(UINT i = 0; i < N; i++)
    nRes = i < N - 1 ? fprintf(pStream, "%10.5lf", vRes[i]) : fprintf(pStream, "%10.5lf\n\n", vRes[i]); 

// Incomplete LU decomposition:
  CSparseMatrix mL(N);
  CSparseMatrix mU(N);
  bool bRes = mtx.decomp_ILU(mL, mU);
  if(bRes)
  {
    output_mtx(pStream, mL, "Low triangle");
    fputs("\n\n", pStream);
    output_mtx(pStream, mU, "Upper triangle");
  }

  CSparseMatrix mLU(N);
  if(mL.multiply(mU, mLU))
  {
    fputs("\n\n", pStream);
    output_mtx(pStream, mLU, "L * U matrix:");
  }

  CSparseMatrix mE(N);
  CSparseMatrix mInvL(N);
  if(inverse_L(mL, mInvL) && mL.multiply(mInvL, mE))
  {
    fputs("\n\n", pStream);
    output_mtx(pStream, mE, "Test of inverse_L() - production of mL * mInvL:");
  }

  mE.clear();
  mE.set_dimension(N);
  CSparseMatrix mInvU(N);
  if(inverse_U(mU, mInvU) && mU.multiply(mInvU, mE))
  {
    fputs("\n\n", pStream);
    output_mtx(pStream, mE, "Test of inverse_U() - production of mU * mInvU:");
  }

  fclose(pStream);
}

void CSparseMatrix::output_mtx(FILE* pStream, const CSparseMatrix& mtx, const char* pHeader)
{
  UINT nElem;
  double fVal;
  int nRes, nZero = 0;
  fputs(pHeader, pStream);
  fputs("\n", pStream);
  UINT N = mtx.m_nRowCount;
  for(size_t i = 0; i < N; i++)
  {
    for(size_t j = 0; j < N; j++)
    {
      nElem = mtx.find_elem(i, j);
      if(nElem == UINT_MAX)
      {
        nRes = j < N - 1 ? fprintf(pStream, "%10d", nZero) : fprintf(pStream, "%10d\n\n", nZero);
        continue;
      }

      fVal = mtx.m_vElems[i].at(nElem).fVal;
      nRes = j < N - 1 ? fprintf(pStream, "%10.5lf", fVal) : fprintf(pStream, "%10.5lf\n\n", fVal);
    }
  }
}

void CSparseMatrix::set_mtx_direct(CSparseMatrix& mtx)
{
  mtx.set_elem(0, 0, 9.0f);
  mtx.set_elem(0, 3, 3.0f);
  mtx.set_elem(0, 4, 1.0f);
  mtx.set_elem(0, 6, 1.0f);

  mtx.set_elem(1, 1, 11.0f);
  mtx.set_elem(1, 2, 2.0f);
  mtx.set_elem(1, 3, 1.0f);
  mtx.set_elem(1, 6, 2.0f);

  mtx.set_elem(2, 1, 1.0f);
  mtx.set_elem(2, 2, 10.0f);
  mtx.set_elem(2, 3, 2.0f);

  mtx.set_elem(3, 0, 2.0f);
  mtx.set_elem(3, 1, 1.0f);
  mtx.set_elem(3, 2, 2.0f);
  mtx.set_elem(3, 3, 9.0f);
  mtx.set_elem(3, 4, 1.0f);

  mtx.set_elem(4, 0, 1.0f);
  mtx.set_elem(4, 3, 1.0f);
  mtx.set_elem(4, 4, 12.0f);
  mtx.set_elem(4, 6, 1.0f);

  mtx.set_elem(5, 5, 8.0f);

  mtx.set_elem(6, 0, 2.0f);
  mtx.set_elem(6, 1, 2.0f);
  mtx.set_elem(6, 4, 3.0f);
  mtx.set_elem(6, 6, 8.0f);
}

void CSparseMatrix::set_mtx_by_rows(CSparseMatrix& mtx)
{
  CVectorMtxElem vElems(7);

  vElems.clear();
  vElems.push_back(CSparseMtxElem(1.0, 6)); vElems.push_back(CSparseMtxElem(1.0, 4)); vElems.push_back(CSparseMtxElem(3.0, 3)); vElems.push_back(CSparseMtxElem(9.0, 0));
  mtx.set_row(0, vElems);

//                mtx.set_elem(1, 6, 2.0f);                 mtx.set_elem(1, 3, 1.0f);                 mtx.set_elem(1, 2, 2.0f);                 mtx.set_elem(1, 1, 11.0f);
  vElems.clear();
  vElems.push_back(CSparseMtxElem(2.0, 6)); vElems.push_back(CSparseMtxElem(1.0, 3)); vElems.push_back(CSparseMtxElem(2.0, 2)); vElems.push_back(CSparseMtxElem(11.0, 1));
  mtx.set_row(1, vElems);

//                mtx.set_elem(2, 3, 2.0f);                mtx.set_elem(2, 2, 10.0f);                 mtx.set_elem(2, 1, 1.0f);
  vElems.clear();
  vElems.push_back(CSparseMtxElem(2.0, 3)); vElems.push_back(CSparseMtxElem(10.0, 2)); vElems.push_back(CSparseMtxElem(1.0, 1));
  mtx.set_row(2, vElems);

//                mtx.set_elem(3, 4, 1.0f);                 mtx.set_elem(3, 3, 9.0f);                 mtx.set_elem(3, 2, 2.0f);                 mtx.set_elem(3, 1, 1.0f);                 mtx.set_elem(3, 0, 2.0f);
  vElems.clear();
  vElems.push_back(CSparseMtxElem(1.0, 4)); vElems.push_back(CSparseMtxElem(9.0, 3)); vElems.push_back(CSparseMtxElem(2.0, 2)); vElems.push_back(CSparseMtxElem(1.0, 1)); vElems.push_back(CSparseMtxElem(2.0, 0));
  mtx.set_row(3, vElems);

//                mtx.set_elem(4, 6, 1.0f);                 mtx.set_elem(4, 4, 12.0f);                 mtx.set_elem(4, 3, 1.0f);                 mtx.set_elem(4, 0, 1.0f);
  vElems.clear();
  vElems.push_back(CSparseMtxElem(1.0, 6)); vElems.push_back(CSparseMtxElem(12.0, 4)); vElems.push_back(CSparseMtxElem(1.0, 3)); vElems.push_back(CSparseMtxElem(1.0, 0));
  mtx.set_row(4, vElems);

//                mtx.set_elem(5, 5, 8.0f);
  vElems.clear();
  vElems.push_back(CSparseMtxElem(8.0, 5));
  mtx.set_row(5, vElems);

//                mtx.set_elem(6, 6, 8.0f);                 mtx.set_elem(6, 4, 3.0f);                 mtx.set_elem(6, 1, 2.0f);                 mtx.set_elem(6, 0, 2.0f);
  vElems.clear();
  vElems.push_back(CSparseMtxElem(8.0, 6)); vElems.push_back(CSparseMtxElem(3.0, 4)); vElems.push_back(CSparseMtxElem(2.0, 1)); vElems.push_back(CSparseMtxElem(2.0, 0));
  mtx.set_row(6, vElems);
}

// END DEBUG

};  // namespace EvaporatingParticle