
#pragma once

#include "stdafx.h"

#include <vector>


namespace EvaporatingParticle
{

typedef unsigned int(__stdcall *CExecProc)(LPVOID);
//---------------------------------------------------------------------------------------
//  CalcThread
//---------------------------------------------------------------------------------------
class CalcThread
{
public:
	CalcThread(UINT nID, CExecProc pThreadFunc, void* pData = NULL);
	~CalcThread();

	void          start_calc_thread();
	void          stop_calc_thread();

  UINT          get_first_job() const;
  UINT          get_last_job() const;
  UINT          get_done_job() const;

  void          done_job();

  void          set_jobs(UINT nFirstJob, UINT nLastJob);

  HANDLE        get_thread_handle() const;

  void*         m_pData;

protected:
	HANDLE        m_hCalcThread;
	UINT          m_nCalcThreadID;  // system ID of the thread, assigned in start_calc_thread().

  UINT          m_nFirstJob,
                m_nLastJob,
                m_nDoneJob;

  UINT          m_nID;  // ordinal number of the thread.

private:
	CExecProc     m_pThreadFunc;
};

inline UINT CalcThread::get_first_job() const
{
  return m_nFirstJob;
}

inline UINT CalcThread::get_last_job() const
{
  return m_nLastJob;
}

inline void CalcThread::done_job()
{
  m_nDoneJob++;
}

inline UINT CalcThread::get_done_job() const
{
  return m_nDoneJob;
}

inline void CalcThread::set_jobs(UINT nFirstJob, UINT nLastJob)
{
  m_nFirstJob = nFirstJob;
  m_nLastJob = nLastJob;
  m_nDoneJob = 0;
}

inline HANDLE CalcThread::get_thread_handle() const
{
  return m_hCalcThread;
}


typedef std::vector<CalcThread*> _CalcThreadVector;
//---------------------------------------------------------------------------------------
//  CalcThreadVector
//---------------------------------------------------------------------------------------
class CalcThreadVector : public _CalcThreadVector
{
public:
  CalcThreadVector();
  ~CalcThreadVector();

  void        distribute_jobs(UINT nFirstJob, UINT nLastJob, CExecProc ThreadFunc, void* pData);

  void        start_execution();
  void        stop_execution();

  bool        wait();

  int         get_progress() const;   // integer percentage of all done jobs.

protected:
  UINT        m_nCount;     // count of the system processors. 
  UINT        m_nJobCount;  // count of jobs to do.
};

}; // namespace EvaporatingParticle