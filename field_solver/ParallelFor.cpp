#include <stdafx.h>
#include "ParallelFor.h"

size_t CParFor::s_nProcNum = 0;

void CParFor::joinAll()
{
	for (auto& t : m_threads) t.join();
}

CParFor::CParFor()
{
	if (s_nProcNum == 0)
		s_nProcNum = std::thread::hardware_concurrency();
	m_threads.resize(s_nProcNum);
}

void TaskQueue::addTask(const Task& task)
{
	Locker lock(m_queueAccess);
	m_tasks.push(task);
}

TaskQueue::Task TaskQueue::getTask()
{
	Locker lock(m_queueAccess);
	Task task = m_tasks.front();
	m_tasks.pop();
	return task;
}

void TaskQueue::clear()
{
	Locker lock(m_queueAccess);
	m_tasks.swap(Tasks());
}

bool TaskQueue::isEmpty()
{
	Locker lock(m_queueAccess);
	return m_tasks.empty();
}

ThreadPool::ThreadPool(size_t nThreadNumber)
	:
	m_nThreadNumber(nThreadNumber ? nThreadNumber : Thread::hardware_concurrency()),
	m_nThreadsActive(0),
	m_bStopFlag(false)
{
	start();
}

ThreadPool::~ThreadPool()
{
	stop();
	for (auto& t : m_threads) t.join();
}

void ThreadPool::threadNumber(size_t nThreadNumber)
{
	stop();
	waitForAll();
	m_nThreadNumber = nThreadNumber;
	start();
}

size_t ThreadPool::threadNumber()
{
	return m_nThreadNumber;
}

void ThreadPool::waitForOne()
{
	Locker lock(m_stopMutex);
	m_stopCondition.wait(lock);
}

void ThreadPool::waitForAll()
{
	Locker lock(m_stopMutex);
	m_stopCondition.wait(lock, [&]()->bool { return m_tasks.isEmpty() && threadsActive() == 0; });
}

void ThreadPool::addTask(const TaskQueue::Task & task)
{
	m_tasks.addTask(task);
	m_startCondition.notify_one();
}

void ThreadPool::splitInPar(size_t n, const std::function<void(size_t)>& atomicOp)
{
	size_t nn = n / threadNumber() + 1;
	if (nn == 1)
	{
		for (size_t i = 0; i < n; ++i)
			addTask([=]() { atomicOp(i); });
	}
	else
	{
		for(size_t i = 0; i < n; i += nn)
			addTask([=]() 
		{
			size_t end = i + nn < n ? i + nn : n;
			for (size_t j = i; j < end; ++j)
				atomicOp(j);
		});
	}
}

void ThreadPool::threadEvtLoop()
{
	while (true)
	{
		{
			Locker lock(m_startMutex);
			m_startCondition.wait(lock, [&]()->bool { return m_tasks.isEmpty(); });
		}

		if (stopFlag()) break;
		else if (m_tasks.isEmpty()) continue;
		else
		{
			TaskQueue::Task task = m_tasks.getTask();
			incThreadsActive();
			task();
			decThreadsActive();
			m_stopCondition.notify_all();
		}
	}
}

void ThreadPool::start()
{
	m_threads = Threads(m_nThreadNumber);
	for (size_t i = 0; i < m_nThreadNumber; ++i)
		m_threads[i] = Thread(&ThreadPool::threadEvtLoop, this);
}

void ThreadPool::stop()
{
	m_tasks.clear();
	stopFlag(true);
	m_startCondition.notify_all();
}

bool ThreadPool::stopFlag() const
{
	return m_bStopFlag;
}

void ThreadPool::stopFlag(bool bStopFlag)
{
	m_bStopFlag = bStopFlag;
}

size_t ThreadPool::threadsActive() const
{
	return m_nThreadsActive;
}

void ThreadPool::decThreadsActive()
{
	--m_nThreadsActive;
}

void ThreadPool::incThreadsActive()
{
	++m_nThreadsActive;
}
