#pragma once
#ifndef _MESH_DATA_
#define _MESH_DATA_

#include <set>
#include <map>
#include <string>
#include <cassert>

#include <linearAlgebra/matrixTemplate.h>

#include "BoundaryMesh.h"

#include "../track3d/CObject.h"
#include "../track3d/Elements.h"
#include "../track3d/vector3d.hpp"
#include "../utilities/ParallelFor.h"

//Mesh edges structure
class CMeshConnectivity
{
public:
	using Label = uint32_t;
	using NodeConnections = std::set<uint32_t, std::less<uint32_t>, Allocator<uint32_t>>;
	using Graph = std::vector<NodeConnections>;
	using Elem = EvaporatingParticle::CElem3D;
	using Node = EvaporatingParticle::CNode3D;
	using Nodes = std::vector<Node*>;
	using Elems = std::vector<Elem*>;

private:
	Graph m_graph;

	//Creates tetrahedral connectivity
	void addTet(Label n0, Label n1, Label n2, Label n3);
	//Creates pyramid connectivity
	void addPyr(Label n0, Label n1, Label n2, Label n3, Label n4);
	//Creates wedge connectivity
	void addWedge(Label n0, Label n1, Label n2, Label n3, Label n4, Label n5);
	//Creates hexahedral connectivity
	void addHexa(Label n0, Label n1, Label n2, Label n3, Label n4, Label n5, Label n6, Label n7);

public:
	//Returns a size of a graph
	Label size() const;
	//Adds new connection
	void addEdge(Label i, Label j);
	//Adds element to a mesh connectivity
	void addElem(const Elem* e);
	//Returns neigbour elements
	const NodeConnections& neighbor(Label i) const;
};

//Linear field transformation interface
class COperator
{
public:
	using Field = std::vector<double>;

	virtual Field applyToField(const Field& f) const = 0;

	//It's important!!!
	virtual ~COperator() {}
};

class CFieldOperator : public COperator
{
public:
	using MatrixCoef = std::pair<uint32_t, double>;
	using MatrixRow = std::map<uint32_t, double,
		std::less<uint32_t>, Allocator<std::pair<uint32_t, double>>>;
	using Matrix = std::vector<MatrixRow>;

	friend class CMeshAdapter;
	friend class DCMeshAdapter;
	friend CFieldOperator& operator+=(CFieldOperator& op1, const CFieldOperator& op2);
	friend CFieldOperator& operator*=(CFieldOperator& op1, const CFieldOperator& op2);
private:
	Matrix m_matrix;

public:

	//Applies operator to a field
	Field applyToField(const Field& f) const;
};

// Addapts mesh to an external usage
class CMeshAdapter
{
public:
	using Elements = EvaporatingParticle::CElementsCollection;
	using Nodes = EvaporatingParticle::CNodesCollection;
	using PBoundary = std::unique_ptr<BoundaryMesh>;
	using InterpCoef = std::pair<uint32_t, double>;
	using InterpCoefs = std::map<uint32_t, double,
		std::less<uint32_t>, Allocator<std::pair<uint32_t, double>>>;
	using Label = uint32_t;
	using Labels = std::vector<Label>;
	using Vector3D = BoundaryMesh::Vector3D;
	using Vector3DOp = std::array<InterpCoefs, 3>;
	using Matrix3D = math::matrix_c<double, 3, 3>;
	using Matrix2D = math::matrix_c<double, 2, 2>;
	using Node = EvaporatingParticle::CNode3D;
	using Element = EvaporatingParticle::CElem3D;
	using ScalarFieldOperator = CFieldOperator;
	using PScalFieldOp = std::unique_ptr<COperator>;
	using ProgressBar = EvaporatingParticle::CObject;
	using PProgressBar = std::unique_ptr<ProgressBar>;
	using Graph = CMeshConnectivity;

	//Keeps base set of parameters for operator instantiation
	struct BaseOperatorParams{};

protected:
	const Elements& m_elems;
	const Nodes& m_nodes;
private:
	PBoundary m_pBoundary;

	//Lazy because it will be created only once by demand
	Graph m_meshGraph;

	//Progress bar interface
	PProgressBar m_pProgressBar;

	//Factor for small step calculation depending on position inside mesh
	double m_fSmallStepFactor;
	static double s_fEpsilon;

	//Looks for an element containing point with coordinates v
	const Element* element(const Vector3D& v, Label& nCurNode, const Label& nPrevNode = Label(0)) const;

public:
	//Math operations with interpolation coeffs
	static InterpCoefs& add(InterpCoefs& ic1, const InterpCoefs& ic2);
	static InterpCoefs& mul(double h, InterpCoefs& ic);
	//Removes zero elements making memory consumption by InterpCoefs smaller
	static InterpCoefs& removeZeros(InterpCoefs& ic);

protected:
	//Creates covariance matrix for dirrections around label 
	Matrix3D covarianceOfDirections(Label l) const;
	//Returns vector of finite difference total projections on a directions x,y,z
	Vector3DOp finDiffDirCov(Label l) const;
	//Returns vector of coefficients for grad in point calculation using directed derivatives averaging
	InterpCoefs gradX(Label l) const;
	InterpCoefs gradY(Label l) const;
	InterpCoefs gradZ(Label l) const;

	//Creates mesh graph connections
	void createGraph();

	//Node type
	enum NodeType : uint8_t
	{
		FirstTypeBoundaryNode,
		SecondTypeBoundaryNode,
		InnerNode
	};
	using NodeTypes = std::vector<NodeType>;
	//Returns array of node types
	NodeTypes nodeTypes() const;
	//Returns true if this node and its neighbour boundary nodes lie in an one plane
	bool isFlatBoundary(Label nNodeIdx, const Vector3D& norm, const NodeTypes& types) const;

	//Create rough and fast LaplacianField solver DU = 0 for a zero approximation
	ScalarFieldOperator laplacianSolver0() const;
	//Simple Laplacian solver with step myltiplication by a factor
	ScalarFieldOperator laplacianSolver1() const;
	//Creates solver which uses graph and operators arithmetics
	ScalarFieldOperator laplacianSolver2() const;
	//Creates solver which uses equal steps
	ScalarFieldOperator laplacianSolver3() const;
	//Calculates laplacian operator
	ScalarFieldOperator laplacian() const;
	//Directed derivative calculation
	ScalarFieldOperator directedDerivative(const Vector3D& dir);

	//Obtains interpolating coefficients if containing element is known
	InterpCoefs interpCoefs(const Vector3D& pos, const Element* e) const;

	//Looks for space position pos in closest neighbor elements
	const Element* lookInClosestElements(const Vector3D& pos, Label l) const;
public:
	CMeshAdapter(const Elements& es, const Nodes& ns, double fSmallStepFactor = 0.3);

	//Returns mesh graph
	const CMeshConnectivity* meshGraph() const { return &m_meshGraph; }

	//Access to a progress bar interface
	ProgressBar* progressBar() const;

	//Gets and sets small step factor value
	double smallStepFactor() const;
	void smallStepFactor(double fVal);
	//Gets and sets epsilon
	static double eps();
	static void eps(size_t nFactor);

	//Returns minimum characteristic size of surrounding elements
	double minElemSize(Label l) const;
	//Returns minimum length of connected edges
	double minEdgeLength(Label l) const;
	//Returns maximum length of connected edges
	double maxEdgeLength(Label l) const;
	//Returns optimal space step that still will be inside neighbor elements
	double optimalStep(const Vector3D& dir, Label l, Label deep = 0) const;

	//Looks for the space coordinate point inside the neighbor elements
	const Element* lookInNeighbor(const Vector3D& pos, Label l, Label deep = 0) const;

	//Returns pointer to a boundary mesh interface
	BoundaryMesh * boundaryMesh() const;
	//Obtains interpolating coefficients for a field value in v-location
	InterpCoefs interpCoefs(const Vector3D& v, Label& nCurNode, const Label& nPrevNode = Label(0)) const;

	//Creates operator for field values on a given mesh
	//Scalar operators first
	enum ScalarOperatorType
	{
		LaplacianSolver0,
		LaplacianSolver1,
		LaplacianSolver2,
		LaplacianSolver3,
		GradX,
		GradY,
		GradZ
	};
	PScalFieldOp createOperator(ScalarOperatorType type = LaplacianSolver1, 
		const BaseOperatorParams* params = nullptr);

	//Returns operators calculating gradient components
	ScalarFieldOperator gradX() const;
	ScalarFieldOperator gradY() const;
	ScalarFieldOperator gradZ() const;
};

#endif // !_MESH_DATA_
