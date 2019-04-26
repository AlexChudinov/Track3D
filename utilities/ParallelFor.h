#pragma once
#ifndef _PAR_FOR_
#define _PAR_FOR_

#include <queue>
#include <mutex>
#include <thread>
#include <vector>
#include <future>
#include "../track3d/CObject.h"



//A pool of threads
class ThreadPool
{
public:
	
	class Mutex
	{
		CRITICAL_SECTION mCritSec;
	public:
		Mutex();
		~Mutex();
		void lock();
		void unlock();
	};

	using Fun = std::function<void()>;
	using Task = std::packaged_task<void()>;
	using Tasks = std::queue<Task>;
	using Future = std::future<void>;
	using Locker = std::unique_lock<Mutex>;
	using Progress = EvaporatingParticle::CObject;

private:
	
	ThreadPool();
	ThreadPool(const ThreadPool&) = delete;
	ThreadPool& operator=(const ThreadPool&) = delete;

	static Mutex mMtx;
	static Tasks mTasks;
public:
	//Returns thread pool global instance
	static ThreadPool& getInstance();
	~ThreadPool();

	//Adds task to task queue
	Future addTask(Fun&& task);

	static void pushTask(Task&& task);
	static Task popTask();

	//Splits array into subarray and does parallel operation on it
	static void splitInPar(size_t n, const std::function<void(size_t)>& atomicOp);
	static void splitInPar(size_t n, std::function<void(size_t)>&& atomicOp, 
		Progress* progress, 
		size_t nThreads = 0);

	//Oarallel computation of max element
	template<class Iterator> static Iterator max_element(Iterator _First, Iterator _Last, 
		std::random_access_iterator_tag);

	template<class Iterator> static Iterator max_element(Iterator _First, Iterator _Last);

private:
	static void CALLBACK workCallback
	(
		PTP_CALLBACK_INSTANCE,
		PVOID,
		PTP_WORK
	);

	PTP_POOL mPool;
	TP_CALLBACK_ENVIRON mCbe;
	DWORD nThreadCnt;

	static Mutex mMutex;

	//Returns a number of current threads
	DWORD threadNumber();
};

template<class Iterator>
inline Iterator ThreadPool::max_element(Iterator _First, Iterator _Last, 
	std::random_access_iterator_tag)
{
	size_t dist = static_cast<size_t>(std::distance(_First, _Last));
	if (dist > getInstance().threadNumber())
	{
		size_t n = dist / getInstance().threadNumber() + 1;
		std::vector<Iterator> results(getInstance().threadNumber());
		std::vector<ThreadPool::Future> futures(getInstance().threadNumber());

		for (size_t i = 0; i < futures.size(); ++i)
		{
			n = n > std::distance(_First, _Last) ? std::distance(_First, _Last) : n;
			Iterator _Next = _First + n;
			futures[i] = getInstance().addTask([&results, i, _First, _Next]()
			{
				results[i] = std::max_element(_First, _Next);
			});
			_First = _Next;
		}

		for (auto& f : futures) f.wait();

		std::vector<Iterator>::iterator pRes = std::max_element
		(
			results.begin(),
			results.end(),
			[](Iterator a, Iterator b)->bool
		{
			return *a < *b;
		}
		);

		return *pRes;
	}
	else
	{
		return std::max_element(_First, _Last);
	}
}

template<class Iterator>
inline Iterator ThreadPool::max_element(Iterator _First, Iterator _Last)
{
	typedef std::iterator_traits<Iterator>::iterator_category category;
	return max_element(_First, _Last, category());
}

#endif // !_PAR_FOR_