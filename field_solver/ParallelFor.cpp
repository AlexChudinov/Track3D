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

TaskQueue::Future TaskQueue::addTask(const Fun& task)
{
	Locker lock(m_queueAccess);
	m_tasks.push(Task(task));
	return m_tasks.back().get_future();
}

TaskQueue::Task TaskQueue::getTask()
{
	Locker lock(m_queueAccess);
	Task task(std::move(m_tasks.front()));
	m_tasks.pop();
	return std::move(task);
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
	m_nThreadNumber = nThreadNumber;
	start();
}

size_t ThreadPool::threadNumber()
{
	return m_nThreadNumber;
}

TaskQueue::Future ThreadPool::addTask(TaskQueue::Fun&& task)
{
	TaskQueue::Future future = m_tasks.addTask(task);
	m_startCondition.notify_one();
	return std::move(future);
}

void ThreadPool::splitInPar(size_t n, const std::function<void(size_t)>& atomicOp)
{
	size_t nn = n / threadNumber() + 1;
	std::vector<TaskQueue::Future> vFutures;
	if (nn == 1)
	{
		for (size_t i = 0; i < n; ++i)
			addTask([=]() { atomicOp(i); });
	}
	else
	{
		for(size_t i = 0; i < n; i += nn)
			vFutures.push_back(
			addTask([=]() 
		{
			size_t end = i + nn < n ? i + nn : n;
			for (size_t j = i; j < end; ++j)
				atomicOp(j);
		}));
	}
	std::for_each(vFutures.begin(), vFutures.end(), 
		[](TaskQueue::Future& future) { future.wait(); });
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
			task();
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
