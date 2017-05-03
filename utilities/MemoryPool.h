#pragma once
#ifndef _MEMORY_POOL_
#define _MEMORY_POOL_

#include <mutex>
#include <forward_list>
#include <algorithm>
#include <vector>

//Interface to an every instance of template of BlockSizePool
class BlockPoolInterface {
	//Pointers to all created instances
	static std::vector<BlockPoolInterface*> s_blockPools;

public:

	//Every instance can be cleaned up
	virtual void cleanUp() = 0;

	//Inserts new one BlockAllocator into container
	static void insert(BlockPoolInterface * p);

	//Cleans up all created template instances
	static void cleanUpEveryPool();
};

//Keeps memory blocks of a certain size
template<size_t nBlockSize>
class BlockSizePool : public BlockPoolInterface
{
public:
	using Mutex = std::mutex;
	using Locker = std::unique_lock<Mutex>;
	using BlocksList = std::forward_list<char*>;
	using BlocksLayout = BlocksList;

private:
	//Cleans all blocks from the pool
	void cleanUpEverything() {
		for (char* p : m_blocks) {
			VirtualFree((void*)p, NULL, MEM_RELEASE);
		}
	}

	//Expands poolsize whenever it is needed
	void expandPoolSize()
	{
		const size_t nSize = nBlockSize * m_nPageSize;
		m_blocks.push_front((char*)VirtualAlloc(NULL, nSize, MEM_COMMIT, PAGE_READWRITE));
		m_layout.push_front(m_blocks.front());

		char
			*pFirst = m_layout.front(),
			*pLast = pFirst + (nSize - nBlockSize);

		for (pFirst += nBlockSize; pFirst < pLast; pFirst += nBlockSize)
			m_layout.push_front(pFirst);
	}

	//Use it as a singletone
	BlockSizePool()
		:
		m_nPageSize(0)
	{
		SYSTEM_INFO sysInfo;
		GetSystemInfo(&sysInfo);
		*const_cast<size_t*>(&m_nPageSize) = sysInfo.dwPageSize;
		BlockPoolInterface::insert(this);
	}

	BlockSizePool(const BlockSizePool&) = delete;
	BlockSizePool(BlockSizePool&&) = delete;
	BlockSizePool& operator=(const BlockSizePool&) = delete;
	BlockSizePool& operator=(BlockSizePool&&) = delete;
public:
	static BlockSizePool& getInstance()
	{
		static BlockSizePool g_blockSizePool;
		return g_blockSizePool;
	}

	//Allocates a one block of a memory
	void * allocBlock()
	{
		Locker lock(m_mtxDataAccess);
		if (m_layout.empty()) expandPoolSize();
		void* res = (void*)m_layout.front();
		m_layout.pop_front();
		return res;
	}

	//Frees a one block of a memory
	void freeBlock(void * pT)
	{
		Locker lock(m_mtxDataAccess);
		m_layout.push_front((char*)pT);
	}

	//Cleans up unused blocks
	void cleanUp()
	{
		m_layout.sort();
		auto it = m_blocks.begin();
		while (it != m_blocks.end()) {
			BlocksList::iterator first = std::lower_bound(m_layout.begin(),
				m_layout.end(), *it);
			if (*first == *it) {
				BlocksList::iterator cur = first, next = first;
				++next;
				size_t nCnt = 0;
				for (; nCnt < m_nPageSize - 1; ++nCnt) {
					//Checks if subsequent blocks were put together
					if (next == m_layout.end() || *next++ - *cur++ != nBlockSize) break;
				}
				if (nCnt == m_nPageSize - 1) {
					//Free all blocks on the page
					m_layout.erase_after(first, next);
					VirtualFree((void*)*it, NULL, MEM_RELEASE);
					first = it;
					it = m_blocks.erase_after(first, ++it);
				}
				else ++it;
			}
			else ++it;
		}
	}

	//Releses all blocks from memory
	~BlockSizePool() { /*cleanUpEverything();*/ }

private:
	BlocksList m_blocks;
	BlocksLayout m_layout;
	const size_t m_nPageSize;
	Mutex m_mtxDataAccess;
};

//Keeps memory blocks of a certain type
template<typename T>
class BlockPool
{
	using UsedBlockSizePool = BlockSizePool<sizeof(T)>;

	//Use it as a singletone
	BlockPool()
		:
		m_bspGlobalRef(UsedBlockSizePool::getInstance())
	{
	}

	BlockPool(const BlockPool&) = delete;
	BlockPool(BlockPool&&) = delete;
	BlockPool& operator=(const BlockPool&) = delete;
	BlockPool& operator=(BlockPool&&) = delete;
public:
	static BlockPool& getInstance()
	{
		static BlockPool g_blockPool;
		return g_blockPool;
	}

	//Allocates a one block of a memory
	T* allocBlock()
	{
		return (T*)m_bspGlobalRef.allocBlock();
	}

	//Frees a one block of a memory
	void freeBlock(T* pT)
	{
		m_bspGlobalRef.freeBlock((void*)pT);
	}

	//Cleans up the global pool of blocks
	void cleanUp() {
		m_bspGlobalRef.cleanUp();
	}

private:
	UsedBlockSizePool& m_bspGlobalRef;
};

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

private:
	std::allocator<T> m_stdAllocator;
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
		return m_stdAllocator.allocate(nCount);
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
	default: m_stdAllocator.deallocate(pT, nCount);
	}
}

template<typename T1, typename T2>
bool operator==(const Allocator<T1>&, const Allocator<T2>&)
{
	return true;
}

#endif // !_MEMORY_POOL_