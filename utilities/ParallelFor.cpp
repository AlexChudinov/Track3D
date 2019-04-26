#include <stdafx.h>
#include "ParallelFor.h"

ThreadPool::Mutex ThreadPool::mMtx;
ThreadPool::Tasks ThreadPool::mTasks;

ThreadPool::ThreadPool()
{
	SYSTEM_INFO sysInfo;
	if (!(mPool = CreateThreadpool(NULL)))
		goto exception;

	GetSystemInfo(&sysInfo);
	nThreadCnt = sysInfo.dwNumberOfProcessors;
	SetThreadpoolThreadMinimum(mPool, nThreadCnt);
	SetThreadpoolThreadMaximum(mPool, nThreadCnt);

	InitializeThreadpoolEnvironment(&mCbe);
	SetThreadpoolCallbackPool(&mCbe, mPool);

	return; 

exception:
	CloseThreadpool(mPool);
	std::string s("Exception #");
	s += std::to_string(GetLastError());
	throw(std::runtime_error(s));
}

ThreadPool & ThreadPool::getInstance()
{
	static ThreadPool sThreadPool;
	return sThreadPool;
}

ThreadPool::~ThreadPool()
{
	CloseThreadpool(mPool);
}

ThreadPool::Future ThreadPool::addTask(Fun&& task)
{
	PTP_WORK threadPoolWork = CreateThreadpoolWork(workCallback, NULL, &mCbe);
	if (!threadPoolWork)
		throw(std::runtime_error(std::string("Error #") + std::to_string(GetLastError())));
	Task work = Task(task);
	Future fut = work.get_future();
	pushTask(std::move(work));
	SubmitThreadpoolWork(threadPoolWork);
	CloseThreadpoolWork(threadPoolWork);
	return fut;
}

void ThreadPool::pushTask(Task&& task)
{
	Locker lock(mMtx);
	mTasks.push(std::move(task));
}

ThreadPool::Task ThreadPool::popTask()
{
	Locker lock(mMtx);
	Task task;
	if (!mTasks.empty())
	{
		task = std::move(mTasks.front());
		mTasks.pop();
	}
	return std::move(task);
}

void ThreadPool::splitInPar(size_t n, const std::function<void(size_t)>& atomicOp)
{
	size_t nn = n / getInstance().threadNumber() + 1;
	std::vector<Future> vFutures;
	if (nn == 1)
	{
		for (size_t i = 0; i < n; ++i)
			vFutures.push_back(getInstance().addTask([=]() { atomicOp(i); }));
	}
	else
	{
		for(size_t i = 0; i < n; i += nn)
			vFutures.push_back(
				getInstance().addTask([=]()
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
	nThreads = nThreads == 0 ? getInstance().threadNumber() : nThreads;
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
			vFutures.push_back(getInstance().addTask([=]() { atomicOp(i); }));;
	}
	else
	{
		for (size_t i = 0; i < n; i += nn)
			vFutures.push_back(
				getInstance().addTask([=]()
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

void ThreadPool::workCallback(PTP_CALLBACK_INSTANCE, PVOID, PTP_WORK)
{
	Task task;
	while ((task = popTask()).valid())
	{
		task();
	}
}

DWORD ThreadPool::threadNumber()
{
	return nThreadCnt;
}

ThreadPool::Mutex::Mutex()
{
	InitializeCriticalSection(&mCritSec);
}

ThreadPool::Mutex::~Mutex()
{
	DeleteCriticalSection(&mCritSec);
}

void ThreadPool::Mutex::lock()
{
	EnterCriticalSection(&mCritSec);
}

void ThreadPool::Mutex::unlock()
{
	LeaveCriticalSection(&mCritSec);
}
