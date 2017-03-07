#pragma once
#ifndef _PAR_FOR_
#define _PAR_FOR_

#include <queue>
#include <mutex>
#include <atomic>
#include <thread>
#include <vector>
#include <future>
#include <algorithm>

//Thread pool tast queue
class TaskQueue
{
public:
	using Fun = std::function<void()>;
	using Task = std::packaged_task<void()>;
	using Tasks = std::queue<Task>;
	using Mutex = std::mutex;
	using Locker = std::unique_lock<Mutex>;
	using Future = std::future<void>;

private:
	Tasks m_tasks;
	Mutex m_queueAccess;

public:
	Future addTask(const Fun& task);
	Task getTask();
	void clear();

	//Returns false if there are some tasks
	bool isEmpty();
};

//A pool of threads
class ThreadPool
{
public:
	using Thread = std::thread;
	using Threads = std::vector<Thread>;
	using Mutex = std::mutex;
	using Locker = std::unique_lock<Mutex>;
	using ConditionVar = std::condition_variable;
	using AtomicBool = std::atomic_bool;

private:
	Mutex m_startMutex;
	ConditionVar m_startCondition;

	size_t m_nThreadNumber;
	AtomicBool m_bStopFlag;
	Threads m_threads;

	TaskQueue m_tasks;
public:

	ThreadPool(size_t nThreadNumber = 0);
	~ThreadPool();

	void threadNumber(size_t nThreadNumber);
	size_t threadNumber();

	//Adds task to task queue
	TaskQueue::Future addTask(TaskQueue::Fun&& task);

	//Splits array into subarray and does parallel operation on it
	void splitInPar(size_t n, const std::function<void(size_t)>& atomicOp);

private:
	void threadEvtLoop();

	//Starts thread event loops
	void start();

	//Stops thread event loops
	void stop();

	//Signals to threads to stop
	bool stopFlag() const;
	void stopFlag(bool bStopFlag);
};

//Splites for each loops into parallel
class CParFor
{
public:
	using Thread = std::thread;
	using Threads = std::vector<Thread>;
	using Mutex = std::mutex;
	using Locker = std::unique_lock<Mutex>;
private:
	static size_t s_nProcNum;
	Threads m_threads;

	void joinAll();
public:
	CParFor();

	template<typename It, typename Task>
	void parallelForEach(It start, It end, Task task);

	template<typename Task>
	void parallelForEach(size_t start, size_t end, Task task);
};

template<typename It, typename Task>
inline void CParFor::parallelForEach(It start, It end, Task task)
{
	size_t n =
		std::distance(start, end),
		nn = n / s_nProcNum + 1;
	size_t 
		first = start,
		last = (first + nn) < end ? first + nn : end;
	for (Thread& t : m_threads)
	{
		t = Thread([=]() {std::for_each(first, last, task)});
		first = last;
		last = (first + nn) < end ? first + nn : end;
	}

	joinAll();
}

template<typename Task>
inline void CParFor::parallelForEach(size_t start, size_t end, Task task)
{
	size_t
		n = end - start,
		nn = n / s_nProcNum + 1,
		first = start,
		last = first + nn < end ? first + nn : end;

	for (size_t i = 0; i < m_threads.size(); ++i)
	{
		m_threads[i] = Thread([=]()
		{
			task(first, last);
		});
		first = last;
		last = first + nn < end ? first + nn : end;
	}

	joinAll();
}

#endif // !_PAR_FOR_
