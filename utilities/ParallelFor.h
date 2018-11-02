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
	
	using Fun = std::function<void()>;
	using Task = std::packaged_task<void()>;
	using Tasks = std::queue<Task>;
	using Future = std::shared_ptr<std::future<void>>;
	using Thread = std::thread;
	using Threads = std::vector<Thread>;
	using Mutex = std::mutex;
	using Locker = std::unique_lock<Mutex>;
	using ConditionVar = std::condition_variable;
	using String = std::string;
	using Progress = EvaporatingParticle::CObject;

private:
	//True if thread pool is valid
	bool m_bValid;

	Mutex m_startMutex;
	ConditionVar m_startCondition;

	size_t m_nThreadNumber;
	bool m_bStopFlag;
	Threads m_threads;

	Tasks m_tasks;

	String m_sErrorDescription;

	Mutex m_globalLock;
	Mutex m_initLock;

	ThreadPool();
	ThreadPool(const ThreadPool&) = delete;
	ThreadPool& operator=(const ThreadPool&) = delete;

	//Returns thread pool global instance
	static ThreadPool& getInstance();
public:
	~ThreadPool();

	//Adds task to task queue
	static Future addTask(Fun&& task);

	//Gets next task from queue
	static Task getTask();

	//Splits array into subarray and does parallel operation on it
	static void splitInPar(size_t n, const std::function<void(size_t)>& atomicOp);
	static void splitInPar(size_t n, std::function<void(size_t)>&& atomicOp, 
		Progress* progress, 
		size_t nThreads = 0);

	//Returns error string from a thread pool
	static String error();

	//Oarallel computation of max element
	template<class Iterator> static Iterator max_element(Iterator _First, Iterator _Last, 
		std::random_access_iterator_tag);

	template<class Iterator> static Iterator max_element(Iterator _First, Iterator _Last);

private:
	void threadEvtLoop();

	//Returns a number of current threads
	static size_t threadNumber();

	//Starts thread event loops
	void start();

	//Stops thread event loops
	void stop();

	//Joins to all current threads
	void joinAll();

	//Initialises thread pool
	void init();
};

#endif // !_PAR_FOR_

template<class Iterator>
inline Iterator ThreadPool::max_element(Iterator _First, Iterator _Last, 
	std::random_access_iterator_tag)
{
	size_t dist = static_cast<size_t>(std::distance(_First, _Last));
	size_t n = dist / ThreadPool::threadNumber() + 1;
	std::vector<Iterator> results(ThreadPool::threadNumber());
	std::vector<ThreadPool::Future> futures(ThreadPool::threadNumber());

	for (size_t i = 0; i < futures.size(); ++i)
	{
		Iterator _Next = _First + n;
		_Next = _Next < _Last ? _Next : _Last;
		futures[i] = ThreadPool::addTask([&results, i, _First, _Next]()
		{
			results[i] = std::max_element(_First, _Next);
		});
		_First = _Next;
	}

	for (auto f : futures) f->wait();

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

template<class Iterator>
inline Iterator ThreadPool::max_element(Iterator _First, Iterator _Last)
{
	typedef std::iterator_traits<Iterator>::iterator_category category;
	return max_element(_First, _Last, category());
}
