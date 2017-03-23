#pragma once
#ifndef _MEMORY_POOL_
#define _MEMORY_POOL_

#include <forward_list>
#include <mutex>

class CMemoryPool
{
public:
	using PVoid = void*;
	using PagesList = std::forward_list<PVoid>;
private:
	const static size_t s_nPageSize;
	//Use it as a singletone
	CMemoryPool();

	CMemoryPool(const CMemoryPool&) = delete;
	CMemoryPool(CMemoryPool&&) = delete;
	CMemoryPool& operator=(const CMemoryPool&) = delete;
	CMemoryPool& operator=(CMemoryPool&&) = delete;
public:
	static CMemoryPool& getInstance();

	//Allocate new page
	void* getPage();

	~CMemoryPool();

	//Returns current size of a one page
	static size_t pageSize();

	//Cleares pool's memory
	void clear();
private:
	PagesList m_pages;
};

//Keeps memory blocks of a certain type
template<typename T>
class BlockPool
{
public:
	using PVoid = void*;
	using Mutex = std::mutex;
	using Locker = std::unique_lock<Mutex>;
private:
	void formatNewPage();

	//Use it as a singletone
	BlockPool();

	BlockPool(const BlockPool&) = delete;
	BlockPool(BlockPool&&) = delete;
	BlockPool& operator=(const BlockPool&) = delete;
	BlockPool& operator=(BlockPool&&) = delete;
public:
	static BlockPool& getInstance();

	//Allocates a one block of a memory
	T* allocBlock();

	//Frees a one block of a memory
	void freeBlock(T* pT);

private:
	const size_t m_nBlockSize;
	const size_t m_nBlocksNum;
	PVoid m_head;
	Mutex m_mutexAccess;
};

template<typename T>
inline void BlockPool<T>::formatNewPage()
{
	PVoid tmp = CMemoryPool::getInstance().getPage();
	m_head = tmp;
	for (size_t i = 0; i < m_nBlocksNum-1; ++i)
	{
		PVoid next = static_cast<char*>(tmp) + m_nBlockSize;
		*static_cast<PVoid*>(tmp) = next;
		tmp = next;
	}
	*static_cast<PVoid*>(tmp) = NULL;
}

template<typename T>
inline BlockPool<T>::BlockPool()
	:
	m_nBlockSize(max(sizeof(T), sizeof(PVoid))),
	m_nBlocksNum(CMemoryPool::pageSize() / m_nBlockSize),
	m_head(NULL)
{
}

template<typename T>
inline BlockPool<T> & BlockPool<T>::getInstance()
{
	static BlockPool<T> s_blockPool;
	return s_blockPool;
}

template<typename T>
inline T * BlockPool<T>::allocBlock()
{
	Locker lock(m_mutexAccess);
	if (!m_head) formatNewPage();
	PVoid res = m_head;
	m_head = *static_cast<PVoid*>(m_head);
	return res;
}

template<typename T>
inline void BlockPool<T>::freeBlock(T * pT)
{
	Locker lock(m_mutexAccess);
	*static_cast<PVoid*>(pT) = m_head;
	m_head = pT;
}

//My allocator class
template<typename T>
class Allocator: std::allocator<T>
{
public:
	//Allocates nCount of objects of type T
	T* allocate(size_t nCount);

	//Deallocates pointer to an array of T objects of size nCount
	void deallocate(T* pT, size_t nCount);
};

template<typename T>
inline T * Allocator<T>::allocate(size_t nCount)
{
	switch (nCount)
	{
	case 0: return nullptr;
	case 1: return BlockPool<T>::getInstance().allocBlock();
	default:
	{
		size_t nBytes = nCount * sizeof(T);
		void* Ptr = VirtualAlloc(NULL, nBytes, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
		if (Ptr && VirtualLock(Ptr, nBytes))
		{
			return static_cast<T*>(Ptr);
		}
		throw std::runtime_error(std::string("Allocator<") + typeid(T).name()
			+ ">::allocate: Allocation error!");
	}
	}
}

template<typename T>
inline void Allocator<T>::deallocate(T * pT, size_t nCount)
{
	switch (nCount)
	{
	case 0: break;
	case 1: if (pT) BlockPool<T>::getInstance().freeBlock(pT); break;
	default: 
		if (pT)
		{
			VirtualUnlock(pT, nCount * sizeof(T));
			VirtualFree(pT, NULL, MEM_RELEASE);
		}
	}
}

#endif // !_MEMORY_POOL_


