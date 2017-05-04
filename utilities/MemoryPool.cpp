#include "stdafx.h"
#include "MemoryPool.h"
#include <vector>
#include "ParallelFor.h"

std::vector<BlockPoolInterface*> BlockPoolInterface::s_blockPools;

void BlockPoolInterface::insert(BlockPoolInterface * p){
	s_blockPools.push_back(p);
}

void BlockPoolInterface::cleanUpEveryPool(){
	std::vector<ThreadPool::Future> futures;
	for (BlockPoolInterface* p : s_blockPools) {
		futures.push_back(
			ThreadPool::getInstance().addTask([&]() { p->cleanUp(); })
		);
	}

	for (auto& f : futures) f.get();
}
