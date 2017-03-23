#include "stdafx.h"
#include "MemoryPool.h"
#include <algorithm>

const size_t CMemoryPool::s_nPageSize = 65536;

CMemoryPool::CMemoryPool()
	:
	m_pages()
{
}

void CMemoryPool::clear()
{
	std::for_each(m_pages.begin(), m_pages.end(),
		[=](PVoid page)
	{
		VirtualFree(page, NULL, MEM_RELEASE);
	});
	m_pages.clear();
}

CMemoryPool & CMemoryPool::getInstance()
{
	static CMemoryPool s_memoryPool;
	return s_memoryPool;
}

void * CMemoryPool::getPage()
{
	PVoid pNextPage = VirtualAlloc(NULL, s_nPageSize, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
	if (pNextPage)
	{
		m_pages.push_front(pNextPage);
		return pNextPage;
	}
	throw std::runtime_error("CMemoryPool::getPage: Memory allocation error.\n");
}

CMemoryPool::~CMemoryPool()
{
	clear();
}

size_t CMemoryPool::pageSize()
{
	return s_nPageSize;
}
