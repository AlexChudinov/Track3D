#pragma once
#ifndef _MEMORY_POOL_
#define _MEMORY_POOL_

#include <mutex>

//Keeps memory blocks of a certain type
template<typename T>
class BlockPool
{
public:
	using Mutex = std::mutex;
	using Locker = std::unique_lock<Mutex>;
private:
	struct BlocksList
	{
		BlocksList * pNext;
	};

	void expandPoolSize();

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

	//Clenups memory pool
	void clean();

	~BlockPool();
private:
	BlocksList * m_head;
	Mutex m_mutexAccess;
	const size_t m_nPageSize; //Size of a one machine page in bytes
};

template<typename T>
inline void BlockPool<T>::expandPoolSize()
{
	constexpr size_t nBlockSize = max(sizeof(T), sizeof(BlocksList*));
	const size_t nSize = nBlockSize * m_nPageSize;
	m_head = (BlocksList*)VirtualAlloc(NULL, nSize, MEM_COMMIT, PAGE_READWRITE);
	BlocksList* head = m_head;

	char
		*pFirst = (char*)head,
		*pLast = pFirst + (nSize - nBlockSize);

	for (; pFirst < pLast; head = head->pNext)
		head->pNext = (BlocksList*)(pFirst += nBlockSize);

	head->pNext = nullptr;
}

template<typename T>
inline BlockPool<T>::BlockPool()
	:
	m_head(nullptr),
	m_nPageSize(0)
{
	SYSTEM_INFO sysInfo;
	GetSystemInfo(&sysInfo);
	*const_cast<size_t*>(&m_nPageSize) = sysInfo.dwPageSize;
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
	if (!m_head) expandPoolSize();
	T* res = (T*) m_head;
	m_head = m_head->pNext;
	return res;
}

template<typename T>
inline void BlockPool<T>::freeBlock(T * pT)
{
	Locker lock(m_mutexAccess);
	((BlocksList*) pT)->pNext = m_head;
	m_head = (BlocksList*) pT;
}

template<typename T>
inline void BlockPool<T>::clean()
{//No clean possibility
}

template<typename T>
inline BlockPool<T>::~BlockPool()
{
	clean();
}

//Class for a one block allocation from pool
template<typename T>
class BlockAllocator
{
public:
	static void* operator new(size_t n)
	{
		if (n != sizeof(T)) return ::operator new(n);
		else return (void*)BlockPool<T>::getInstance().allocBlock();
	}

	static void* operator new(size_t, void* m) { return m; }

	static void operator delete(void* m, size_t n)
	{
		if (n != sizeof(T)) return ::operator delete(m, n);
		else return BlockPool<T>::getInstance().freeBlock((T*)m);
	}

	static void operator delete(void*, void*){}
};

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
		return (T*)::operator new(sizeof(T)*nCount);
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
	default: ::operator delete((void*)pT, sizeof(T)*nCount);
	}
}

template<typename T1, typename T2>
bool operator==(const Allocator<T1>&, const Allocator<T2>&)
{
	return true;
}

#endif // !_MEMORY_POOL_