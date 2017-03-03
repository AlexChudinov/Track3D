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

bool TaskQueue::isEmpty() const
{
	return m_tasks.empty();
}

void ThreadPool::threadEvtLoop()
{
	Locker lock(m_startMutex);
	while (true)
	{
		m_startCondition.wait(lock);
	}
}
