#pragma once
#ifndef _DCMESH_ADAPTER_
#define _DCMESH_ADAPTER_

#include "MeshData.h"
#include "DirichletTesselation.h"

/**
 * @class DCMeshAdapter implementation of CMeshAdapter class for the dirichlet tesselation ussage
 */
class DCMeshAdapter : public CMeshAdapter
{
	using DirTess = EvaporatingParticle::CDirichletTesselation;
public:
	DCMeshAdapter(const Elements& es, const Nodes& ns, const DirTess& tess);

	//Creates operator for field values on a given mesh
	//Scalar operators first
	enum ScalarOperatorType
	{
		LaplacianSolver,
		LaplacianSolver1,
		GradX,
		GradY,
		GradZ
	};

	PScalFieldOp createOperator(ScalarOperatorType type, const BaseOperatorParams*);

	
	//Calculate gradient at a given node
	InterpCoefs gradX(Label idx) const;
	InterpCoefs gradY(Label idx) const;
	InterpCoefs gradZ(Label idx) const;

	//Returns operators calculating gradient components
	ScalarFieldOperator gradX() const;
	ScalarFieldOperator gradY() const;
	ScalarFieldOperator gradZ() const;
	
	ScalarFieldOperator laplacian() const;
	ScalarFieldOperator laplacian1() const;
private:
	const DirTess& m_tess;
};

#endif