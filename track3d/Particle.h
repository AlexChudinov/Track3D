#pragma once
#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "../utilities/Integrator.h"
#include "../utilities/MemoryPool.h"

class Particle : public IntegratorState
{
public:
	static void operator delete(void* ptr, size_t n);
	virtual void deleteObj() = 0;
};

#endif // !_PARTICLE_H_
