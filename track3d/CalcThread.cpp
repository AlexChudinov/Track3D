
#include "stdafx.h"

#include "CalcThread.h"
#include "ParticleTracking.h"
#include "Tracker.hpp"

namespace EvaporatingParticle
{

//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
CalcThread::CalcThread(UINT nID, CExecProc pThreadFunc, void* pData)
  : m_nID(nID), m_pThreadFunc(pThreadFunc), m_pData(pData), m_hCalcThread(NULL), m_nCalcThreadID(0)
{
}

CalcThread::~CalcThread()
{
  stop_calc_thread();
}

void CalcThread::start_calc_thread()
{
  if(m_hCalcThread == NULL)
    m_hCalcThread = (HANDLE)_beginthreadex(0, 0, m_pThreadFunc, this, 0, &m_nCalcThreadID);

// Threads priority support:
  CTracker* pObj = (CTracker*)m_pData;
  if(pObj->get_use_multi_thread())
  {
    BOOL bRes = FALSE;
    int nPriority = pObj->get_calc_thread_priority();
    switch(nPriority)
    {
      case CTracker::thpNormal: bRes = ::SetThreadPriority(m_hCalcThread, THREAD_PRIORITY_NORMAL); break;
      case CTracker::thpBelowNorm: bRes = ::SetThreadPriority(m_hCalcThread, THREAD_PRIORITY_BELOW_NORMAL); break;
      case CTracker::thpLowest: bRes = ::SetThreadPriority(m_hCalcThread, THREAD_PRIORITY_LOWEST); break;
    }
  }
}

void CalcThread::stop_calc_thread()
{
  if(m_hCalcThread != NULL)
  {
    WaitForSingleObject(m_hCalcThread, INFINITE);
    m_hCalcThread = NULL;
  }
}


//---------------------------------------------------------------------------------------
//
//---------------------------------------------------------------------------------------
CalcThreadVector::CalcThreadVector()
{
  SYSTEM_INFO sysinfo;
  GetSystemInfo(&sysinfo);
  m_nCount = sysinfo.dwNumberOfProcessors;
}

CalcThreadVector::~CalcThreadVector()
{
  size_t nSize = size();
  if(nSize != 0)
  {
    for(size_t i = 0; i < nSize; i++)
    {
      CalcThread* pThread = at(i);
      pThread->stop_calc_thread();
      delete pThread;
    }
  }

  clear();
}

void CalcThreadVector::distribute_jobs(UINT nFirstJob, UINT nLastJob, CExecProc pThreadFunc, void* pData)
{
  m_nJobCount = 1 + nLastJob - nFirstJob;
  if(m_nJobCount <= m_nCount) // one job for each thread.
  {
    for(UINT i = 0; i < m_nJobCount; i++)
    {
      CalcThread* pThread = new CalcThread(i, pThreadFunc, pData);
      pThread->set_jobs(i, i);
      push_back(pThread);
    }
  }
  else
  {
    UINT nStartJob, nEndJob;
    for(UINT i = 0; i < m_nCount; i++)
    {
      nStartJob = i == 0 ? nFirstJob : nEndJob + 1;
      nEndJob = (i + 1) * (nLastJob - nFirstJob) / m_nCount;
      CalcThread* pThread = new CalcThread(i, pThreadFunc, pData);
      pThread->set_jobs(nStartJob, nEndJob);
      push_back(pThread);
    }
  }
}

void CalcThreadVector::start_execution()
{
  size_t nSize = size();
  for(size_t i = 0; i < nSize; i++)
  {
    CalcThread* pThread = at(i);
    pThread->start_calc_thread();
  }
}

void CalcThreadVector::stop_execution()
{
  CTracker* pObj = CParticleTrackingApp::Get()->GetTracker();
  pObj->terminate();
  wait();
}

bool CalcThreadVector::wait()
{
  size_t nSize = size();
  HANDLE* pHandles = new HANDLE[nSize];
  for(size_t i = 0; i < nSize; i++)
  {
    CalcThread* pThread = at(i);
    pHandles[i] = pThread->get_thread_handle();
  }

  DWORD nRes = WaitForMultipleObjects(nSize, pHandles, TRUE, INFINITE);
  delete[] pHandles;
  return nRes != WAIT_TIMEOUT;
}

int CalcThreadVector::get_progress() const
{
  int nJobsDone = 0;
  size_t nSize = size();
  for(size_t i = 0; i < nSize; i++)
    nJobsDone += at(i)->get_done_job();

  return int(0.5 + 100 * ((double)nJobsDone / m_nJobCount));
}

}; // namespace EvaporatingParticle