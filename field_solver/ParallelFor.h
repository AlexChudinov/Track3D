#pragma once
#ifndef _PAR_FOR_
#define _PAR_FOR_

#include <queue>
#include <mutex>
#include <thread>
#include <vector>
#include <future>
#include <algorithm>
#include "../track3d/CObject.h"

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

	TaskQueue m_tasks;

	String m_sErrorDescription;

	Mutex m_globalLock;

	ThreadPool();
	ThreadPool(const ThreadPool&) = delete;
	ThreadPool& operator=(const ThreadPool&) = delete;

public:
	~ThreadPool();

	//Returns thread pool global instance
	static ThreadPool& getInstance();

	void threadNumber(size_t nThreadNumber);
	size_t threadNumber();

	//Adds task to task queue
	TaskQueue::Future addTask(TaskQueue::Fun&& task);

	//Splits array into subarray and does parallel operation on it
	void splitInPar(size_t n, const std::function<void(size_t)>& atomicOp);
	void splitInPar(size_t n, const std::function<void(size_t)>& atomicOp, Progress* progress);

	//Returns error string from a thread pool
	const String& error() const;

private:
	void threadEvtLoop();

	//Starts thread event loops
	void start();

	//Stops thread event loops
	void stop();

	//Signals to threads to stop
	bool stopFlag() const;
	void stopFlag(bool bStopFlag);

	//Joins to all current threads
	void joinAll();

	//Initialises thread pool
	void init();

	//Returns true if the thread pool is valid
	bool valid();
	void valid(bool bValid);
};

#endif // !_PAR_FOR_
