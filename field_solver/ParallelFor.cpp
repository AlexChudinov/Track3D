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
