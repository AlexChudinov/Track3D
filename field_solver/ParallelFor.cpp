#include <stdafx.h>
#include "ParallelFor.h"

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

ThreadPool::ThreadPool()
	:
	m_bValid(false),
	m_nThreadNumber(Thread::hardware_concurrency()),
	m_bStopFlag(false)
{}

ThreadPool::~ThreadPool()
{
	stop();
}

ThreadPool& ThreadPool::getInstance()
{
	static ThreadPool g_pThreadPool;
	if (!g_pThreadPool.valid())
	{
		g_pThreadPool.init();
	}
	return g_pThreadPool;
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

void ThreadPool::splitInPar(size_t n, const std::function<void(size_t)>& atomicOp, Progress * progress)
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
		for (size_t i = 0; i < n; i += nn)
			vFutures.push_back(
				addTask([=]()
		{
			size_t end = i + nn < n ? i + nn : n;
			for (size_t j = i; j < end; ++j)
				atomicOp(j);
		}));
	}
	size_t nProgressVal = 0;
	progress->set_progress(nProgressVal);
	std::for_each(vFutures.begin(), vFutures.end(),
		[&](TaskQueue::Future& future) 
	{ 
		future.wait();
		progress->set_progress(++nProgressVal * 100 / threadNumber());
	});
}

const std::string & ThreadPool::error() const
{
	return m_sErrorDescription;
}

void ThreadPool::threadEvtLoop()
{
	try {
		while (!stopFlag())
		{
			Locker lock(m_startMutex);
			m_startCondition.wait(lock, [&]()->bool { return !m_tasks.isEmpty() || stopFlag(); });

			if (stopFlag()) break;
			else
			{
				TaskQueue::Task task = m_tasks.getTask();
				lock.unlock();
				task();
			}
		}
	}
	catch (const std::exception& ex)
	{
		Locker lock(m_startMutex);
		(m_sErrorDescription += "\n") += ex.what();
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
	joinAll();
}

bool ThreadPool::stopFlag() const
{
	return m_bStopFlag;
}

void ThreadPool::stopFlag(bool bStopFlag)
{
	m_bStopFlag = bStopFlag;
}

void ThreadPool::joinAll()
{
	std::for_each(m_threads.begin(), m_threads.end(), [](Thread& t) { t.join(); });
}

void ThreadPool::init()
{
	Locker lock(m_globalLock);
	if (!m_bValid)
	{
		m_bValid = true;
		start();
	}
}

bool ThreadPool::valid()
{
	Locker lock(m_globalLock);
	return m_bValid;
}

void ThreadPool::valid(bool bValid)
{
	Locker lock(m_globalLock);
	m_bValid = bValid;
}
