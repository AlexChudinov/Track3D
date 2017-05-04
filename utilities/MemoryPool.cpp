#include "stdafx.h"
#include "MemoryPool.h"
#include <vector>
#include "ParallelFor.h"

std::vector<BlockPoolInterface*> BlockPoolInterface::s_blockPools;

void BlockPoolInterface::insert(BlockPoolInterface * p){
	s_blockPools.push_back(p);
}

void BlockPoolInterface::cleanUpEveryPool(){
	ThreadPool::getInstance().addTask([&]() {
		for (BlockPoolInterface* p : s_blockPools) {
			p->cleanUp();
		}
	});
}
