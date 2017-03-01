#pragma once
#ifndef _PAR_FOR_
#define _PAR_FOR_

#include <thread>
#include <vector>
#include <algorithm>

//Splites for each loops into paralel
class CParFor
{
public:
	using Thread = std::thread;
	using Threads = std::vector<Thread>;
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
		t = Thread([&]() {std::for_each(first, last, task)});
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
