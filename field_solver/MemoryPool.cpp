#include "stdafx.h"
#include "MemoryPool.h"

PagePool::PagePool()
	:
	m_nPageSize(4096)
{
	SYSTEM_INFO info;
	GetSystemInfo(&info);
	*const_cast<size_t*>(&m_nPageSize) = info.dwPageSize;
}

PagePool & PagePool::getInstance()
{
	static PagePool s_pagePool;
	return s_pagePool;
}

void * PagePool::allocatePage()
{
	return VirtualAlloc(NULL, m_nPageSize * NUM_OF_PAGES_TO_COMMIT, MEM_COMMIT, PAGE_READWRITE);
}

size_t PagePool::pageSize() const
{
	return m_nPageSize;
}
