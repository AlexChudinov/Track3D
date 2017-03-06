#pragma once
#ifndef _PAR_FOR_
#define _PAR_FOR_

#include <queue>
#include <mutex>
#include <atomic>
#include <thread>
#include <vector>
#include <algorithm>

//Thread pool tast queue
class TaskQueue
{
public:
	using Task = std::function<void()>;
	using Tasks = std::queue<Task>;
	using Mutex = std::mutex;
	using Locker = std::unique_lock<Mutex>;
private:
	Tasks m_tasks;
	Mutex m_queueAccess;
public:
	void addTask(const Task& task);
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
	using AtomicUint = std::atomic_uint32_t;
private:
	Mutex m_startMutex;
	ConditionVar m_startCondition;
	Mutex m_stopMutex;
	ConditionVar m_stopCondition;
	Mutex m_mutex; //just mutex

	size_t m_nThreadNumber;
	AtomicUint m_nThreadsActive;
	AtomicBool m_bStopFlag;
	Threads m_threads;

	TaskQueue m_tasks;
public:

	ThreadPool(size_t nThreadNumber = 0);
	~ThreadPool();

	void threadNumber(size_t nThreadNumber);
	size_t threadNumber();

	//Waits when one thread is ending the work
	void waitForOne();

	//Waits when all thread is finished the work
	void waitForAll();

	//Adds task to task queue
	void addTask(const TaskQueue::Task& task);

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

	//Counts active threads
	size_t threadsActive() const;
	void decThreadsActive();
	void incThreadsActive();
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
