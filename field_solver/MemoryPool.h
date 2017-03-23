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
	return static_cast<T*>(res);
}

template<typename T>
inline void BlockPool<T>::freeBlock(T * pT)
{
	Locker lock(m_mutexAccess);
	*(PVoid*)pT = m_head;
	m_head = pT;
}

//My allocator class
template<typename T>
class Allocator
{
public:
	static_assert(!std::is_const<T>::value,
		"The C++ Standard forbids containers of const elements "
		"because allocator<const T> is ill-formed.");

	typedef void _Not_user_specialized;

	typedef T value_type;

	typedef value_type *pointer;
	typedef const value_type *const_pointer;

	typedef value_type& reference;
	typedef const value_type& const_reference;

	typedef size_t size_type;
	typedef ptrdiff_t difference_type;

	typedef std::true_type propagate_on_container_move_assignment;
	typedef std::true_type is_always_equal;

	template<class _Other>
	struct rebind
	{	// convert this type to allocator<_Other>
		typedef Allocator<_Other> other;
	};

	pointer address(reference _Val) const _NOEXCEPT
	{	// return address of mutable _Val
		return (_STD addressof(_Val));
	}

	const_pointer address(const_reference _Val) const _NOEXCEPT
	{	// return address of nonmutable _Val
		return (_STD addressof(_Val));
	}

	Allocator() _THROW0()
	{	// construct default allocator (do nothing)
	}

	Allocator(const Allocator<T>&) _THROW0()
	{	// construct by copying (do nothing)
	}

	template<class _Other>
	Allocator(const Allocator<_Other>&) _THROW0()
	{	// construct from a related allocator (do nothing)
	}

	template<class _Other>
	Allocator<T>& operator=(const Allocator<_Other>&)
	{	// assign from a related allocator (do nothing)
		return (*this);
	}

	//Allocates nCount of objects of type T
	pointer allocate(size_type nCount);

	//Deallocates pointer to an array of T objects of size nCount
	void deallocate(pointer pT, size_type nCount);

	pointer allocate(size_type _Count, const void *)
	{	// allocate array of _Count elements, ignore hint
		return (allocate(_Count));
	}

	template<class _Objty,
		class... _Types>
		void construct(_Objty *_Ptr, _Types&&... _Args)
	{	// construct _Objty(_Types...) at _Ptr
		::new ((void *)_Ptr) _Objty(_STD forward<_Types>(_Args)...);
	}


	template<class _Uty>
	void destroy(_Uty *_Ptr)
	{	// destroy object at _Ptr
		_Ptr->~_Uty();
	}

	size_t max_size() const _NOEXCEPT
	{	// estimate maximum array size
		return ((size_t)(-1) / sizeof(T));
	}
};

template<typename T>
inline T * Allocator<T>::allocate(size_type nCount)
{
	switch (nCount)
	{
	case 0: return nullptr;
	case 1: return BlockPool<T>::getInstance().allocBlock();
	default:
	{
		size_t nBytes = nCount * sizeof(T);
		void* Ptr = VirtualAlloc(NULL, nBytes, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
		if (Ptr)
		{
			return static_cast<T*>(Ptr);
		}
		throw std::runtime_error(std::string("Allocator<") + typeid(T).name()
			+ ">::allocate: Allocation error!");
	}
	}
}

template<typename T>
inline void Allocator<T>::deallocate(pointer pT, size_type nCount)
{
	switch (nCount)
	{
	case 0: break;
	case 1: if (pT) BlockPool<T>::getInstance().freeBlock(pT); break;
	default: if (pT) VirtualFree(pT, NULL, MEM_RELEASE);
	}
}

template<typename T1, typename T2>
bool operator==(const Allocator<T1>&, const Allocator<T2>&)
{
	return true;
}

#endif // !_MEMORY_POOL_