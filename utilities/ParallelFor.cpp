#include <stdafx.h>
#include "ParallelFor.h"

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
	g_pThreadPool.init();
	return g_pThreadPool;
}

ThreadPool::Future ThreadPool::addTask(Fun&& task)
{
	Locker lock(m_globalLock);
	m_tasks.push(Task(task));
	m_startCondition.notify_one();
	return m_tasks.back().get_future();
}

ThreadPool::Task ThreadPool::getTask()
{
	Locker lock(m_globalLock);
	Task task(std::move(m_tasks.front()));
	m_tasks.pop();
	return task;
}

void ThreadPool::splitInPar(size_t n, const std::function<void(size_t)>& atomicOp)
{
	size_t nn = n / m_nThreadNumber + 1;
	std::vector<Future> vFutures;
	if (nn == 1)
	{
		for (size_t i = 0; i < n; ++i)
			vFutures.push_back(addTask([=]() { atomicOp(i); }));
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
		[](Future& future) { future.wait(); });
}

void ThreadPool::splitInPar(size_t n, std::function<void(size_t)>&& atomicOp, Progress * progress, size_t nThreads)
{
	nThreads = nThreads == 0 ? m_nThreadNumber : nThreads;
	size_t
		nn = n / nThreads + 1,
		nProgressStep = n / 90,
		nProgressCounter = 0;
	Mutex mutex;
	std::vector<Future> vFutures;

	auto setProgress = [&]()->bool
	{
		Locker lock(mutex);
		progress->set_progress(++nProgressCounter);
		return progress->get_terminate_flag();
	};

	if (nn == 1)
	{
		for (size_t i = 0; i < n; ++i)
			vFutures.push_back(addTask([=]() { atomicOp(i); }));;
	}
	else
	{
		for (size_t i = 0; i < n; i += nn)
			vFutures.push_back(
				addTask([=]()
		{
			size_t end = i + nn < n ? i + nn : n;
			for (size_t j = i; j < end; ++j)
			{
				if ((j - i) % nProgressStep == 0 && setProgress()) break;
				atomicOp(j);
			}
		}));
	}
	std::for_each(vFutures.begin(), vFutures.end(), [&](Future& future) { future.wait(); });
}

std::string ThreadPool::error()
{
	Locker lock(m_globalLock);
	return m_sErrorDescription;
}

void ThreadPool::threadEvtLoop()
{
	try 
	{
		while (!m_bStopFlag)
		{
			Locker lock(m_startMutex);
			m_startCondition.wait(lock, [&]()->bool { return !empty() || m_bStopFlag; });

			if (m_bStopFlag) break;
			else
			{
				Task task = getTask();
				lock.unlock();
				task();
			}
		}
	}
	catch (const std::exception& ex)
	{
		Locker lock(m_globalLock);
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
	Locker lock(m_globalLock);
	std::swap(m_tasks, Tasks());
	m_bStopFlag = true;
	m_startCondition.notify_all();
	lock.unlock();
	joinAll();
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

bool ThreadPool::empty()
{
	Locker lock(m_globalLock);
	return m_tasks.empty();
}
