#pragma once
#ifndef _MESH_DATA_
#define _MESH_DATA_

#include <set>
#include <map>
#include <string>
#include <vector>
#include <memory>
#include <numeric>
#include <cassert>
#include <algorithm>
#include <unordered_map>

#include <linearAlgebra\matrixTemplate.h>
#include "ParallelFor.h"
#include "../track3d/Elements.h"
#include "../track3d/vector3d.hpp"
#include "../track3d/CObject.h" // [MS] 10-02-2017 progress bar support.

//Boundary conditions for the mesh
class BoundaryMesh
{
	using StrRef = std::reference_wrapper<const std::string>;

	//Declare reference wrappers less operator
	class StrRefLess : public std::binary_function<const StrRef&, const StrRef&, bool>
	{
	public: bool operator()(const StrRef& s1, const StrRef& s2) const{ return s1.get() < s2.get(); }
	};

	template<typename FieldType> using BoundaryValsEntry = std::map<uint32_t, FieldType>;
	template<typename FieldType> using BoundaryValsTab = std::map<std::string, BoundaryValsEntry<FieldType>>;

	class BoundaryValsBase
	{
	public: 
		virtual void removeBoundary(const std::string&) = 0;
		virtual const std::type_info& fieldType() const = 0;
	};
	//Keeps boundary values of field
	template<typename FieldType>
	class BoundaryVals : public BoundaryValsBase
	{
		BoundaryValsTab<FieldType> m_boundariesVal;
	public:
		void addBoundaryVals(const std::string& strName, const BoundaryValsEntry<FieldType>& mpEntry)
		{
			m_boundariesVal[strName] = mpEntry;
		}

		void removeBoundary(const std::string& strName) { m_boundariesVal.erase(strName); }

		FieldType& getVal(const std::string& strName, uint32_t nLabel)
		{ 
			return m_boundariesVal.at(strName).at(nLabel);
		}

		const FieldType& getVal(const std::string& strName, uint32_t nLabel) const
		{
			return m_boundariesVal.at(strName).at(nLabel);
		}

		const std::type_info& fieldType() const { return typeid(FieldType); }
	};

public:
	enum BoundaryType { ZERO_GRAD, FIXED_VAL };

	using Vector3D = EvaporatingParticle::Vector3<double>;
	using SetLabels = std::set<uint32_t>;
	using BoundaryDescription = std::pair<BoundaryType, SetLabels>;
	using BoundariesMap = std::map<std::string, BoundaryDescription>;
	using NamesList = std::set<StrRef, StrRefLess>;
	using ReversedBoundariesMap = std::map<uint32_t, std::pair<Vector3D, NamesList>>;
	using BoundaryNormals = std::map<uint32_t, Vector3D>;
	using PBoundaryVals = std::unique_ptr<BoundaryValsBase>;

	using iterator = ReversedBoundariesMap::iterator;
	using const_iterator = ReversedBoundariesMap::const_iterator;

private:
	BoundariesMap m_mapBoundariesList;
	ReversedBoundariesMap m_mapReversedBoundariesList;
	PBoundaryVals m_pBoundaryVals;

public:
	//Checks if boundary is empty
	bool empty() const; 
	//Returns biggest number of boundary
	uint32_t maxLabel() const;

	//Adds new boundary patch
	void addBoundary(
		const std::string& strName,
		const std::vector<uint32_t>& vLabels,
		const std::vector<Vector3D>& vNormals,
		BoundaryType type = FIXED_VAL);

	//Removes existing boundary patch
	void removeBoundary(const std::string& strName);

	//Sets type of a boundary with a name strName
	void boundaryType(const std::string& strName, BoundaryType type);

	//Returns type of a boundary with a name strName
	BoundaryType boundaryType(const std::string& strName) const;

	//Returns a set of boundaries connected to a given label
	const NamesList& boundaryNames(uint32_t l) const;

	//Returns numbers of nodes which are belong to a given boundary
	const SetLabels& boundaryLabels(const std::string& strName) const;

	//Get iterators for boundary
	const_iterator begin() const;
	iterator begin();
	const_iterator end() const;
	iterator end();

	//Checks if the name is in the boundaries list
	bool isBoundary(const std::string& sName) const;
	//Checks if the label belongs to the boundary
	bool isBoundary(uint32_t l) const;
	
	//Checks if it is a first-type (Dirichlet) boundary condition
	bool isFirstType(const std::string& sName) const;	
	bool isFirstType(uint32_t l) const;

	//Gets the normal for given node label
	Vector3D normal(uint32_t l) const;

	//Returns size of a patch strName
	size_t patchSize(const std::string& strName) const;

	//Sets boundary vals
	template<typename FieldType>
	void boundaryVals(const std::string& strName, const std::vector<FieldType>& vVals)
	{
		if (!isBoundary(strName)) 
			throw std::runtime_error(std::string("BoundaryMesh::addBoundaryVals<")
				+ typeid(FieldType).name() + ">"
				+ ": " + strName + " boundary was not found.");
		const SetLabels& vBoundaryLabels = m_mapBoundariesList[strName].second;
		if (vVals.size() != vBoundaryLabels.size())
			throw std::runtime_error(std::string("BoundaryMesh::addBoundaryVals<")
				+ typeid(FieldType).name() + ">"
				+ ": Input vector of values and vector of boundary labels have different sizes.");
		BoundaryValsEntry<FieldType> vals;
		SetLabels::const_iterator itLabels = vBoundaryLabels.begin();
		for (size_t i = 0; i < vVals.size(); ++i)
			vals[*(itLabels++)] = vVals[i];
		if (!m_pBoundaryVals) m_pBoundaryVals.reset(new BoundaryVals<FieldType>);
		else
		{
			if (typeid(FieldType) != m_pBoundaryVals->fieldType())
				throw std::runtime_error(std::string("BoundaryMesh::addBoundaryVals<")
					+ typeid(FieldType).name() + ">: Try to call the function with pointer to"
					" BoundaryVals<" + m_pBoundaryVals->fieldType().name() + "> object.");
		}
		dynamic_cast<BoundaryVals<FieldType>*>(m_pBoundaryVals.get())->addBoundaryVals(strName, vals);
	}

	//Applies boundary values to a field
	template<typename FieldType>
	void applyBoundaryVals(std::vector<FieldType>& f) const
	{
		if (!m_pBoundaryVals) return;
		if (m_pBoundaryVals->fieldType() != typeid(FieldType))
			throw std::runtime_error(std::string("BoundaryMesh::applyBoundaryVals: ")
				+ "Try to apply boundary vals of type " + m_pBoundaryVals->fieldType().name()
				+ " to field of type" + typeid(FieldType).name());
		if (maxLabel() >= f.size())
			throw std::runtime_error("BoundaryMesh::applyBoundaryVals: "
				"Field size is too small.");

		for (const auto& boundaryPatch : *this)
		{
			int primaryCondition = 0;
			FieldType primaryCondAcc(0.0);
			const NamesList& listNames = boundaryPatch.second.second;
			for (const auto& name : listNames)
			{
				switch (boundaryType(name))
				{
				case BoundaryMesh::ZERO_GRAD:
					break;
				case BoundaryMesh::FIXED_VAL:
					primaryCondAcc += 
						dynamic_cast<const BoundaryVals<FieldType>*>(m_pBoundaryVals.get())
						->getVal(name, boundaryPatch.first);
					primaryCondition++;
					break;
				default:
					throw std::runtime_error(std::string("BoundaryMesh::applyBoundaryVals<")
						+ typeid(FieldType).name() + ">: Unexpected boundary condition type.");
				}
			}
			if (primaryCondition != 0)
				f[boundaryPatch.first] = primaryCondAcc / primaryCondition;
			else
				f[boundaryPatch.first] = 0;
		}
	}
};

//Mesh edges structure
class CMeshConnectivity
{
public:
	using Label = uint32_t;
	using NodeConnections = std::set<uint32_t>;
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

//Field implementation
template<typename FieldType>
class CField
{
public:
	using FieldData = std::vector<FieldType>;
private:
	FieldData m_vFieldVals;
public:
	FieldData& fieldVals() { return m_vFieldVals; }
	const FieldData& fieldVals() const { return m_vFieldVals; }

	//Returns standart square deviation of two fields
	static double std(const FieldData& f1, const FieldData& f2)
	{
		assert(f1.size() == f2.size());
		return std::inner_product(f1.begin(), f1.end(), f2.begin, 0.0, std::plus<double>(),
			[](const FieldType& fv1, const FieldType& fv2)->double
		{
			FieldType df = fv2 - fv1;
			return df*df;
		});
	}
	static double std(const CField& f1, const CField& f2) 
	{ 
		return std(f1.fieldVals(), f2.fieldVals()); 
	}
	//Returns maximum square difference of two fields
	static double maxSquareDiff(const FieldData& f1, const FieldData& f2)
	{
		assert(f1.size() == f2.size());
		return std::inner_product(f1.begin(), f1.end(), f2.begin(), 0.0, 
			[](const double& d1, const double& d2)->double { return max(d1, d2); },
			[](const FieldType& fv1, const FieldType& fv2)->double
		{
			FieldType df = fv2 - fv1;
			return df*df;
		});
	}
	static double maxSquareDiff(const CField& f1, const CField& f2) 
	{ 
		return maxSquareDiff(f1.fieldVals(), f2.fieldVals()); 
	}
};

//Linear field transformation
template<typename FieldType>
class CFieldOperator
{
public:
	using MatrixCoef = std::pair<uint32_t, double>;
	using MatrixRow = std::unordered_map<uint32_t, double>;
	using Matrix = std::vector<MatrixRow>;
	using Field = std::vector<FieldType>;

	friend class CMeshAdapter;

private:
	Matrix m_matrix;

public:

	//Applies operator to a field
	Field applyToField(const Field& f)
	{
		if (m_matrix.size() != f.size())
			throw std::runtime_error("CFieldOperator::applyToField:"
				" Matrix and field sizes are different!");
		Field result(f.size(), 0.0);

		ThreadPool::getInstance().splitInPar(f.size(), [&](size_t n)
		{
			for (const MatrixCoef& c : m_matrix[n])
				result[n] += f[c.first] * c.second;
		});

		return result;
	}
};

// Addapts mesh to an external usage
class CMeshAdapter
{
public:
	using Elements = EvaporatingParticle::CElementsCollection;
	using Nodes = EvaporatingParticle::CNodesCollection;
	using PBoundary = std::unique_ptr<BoundaryMesh>;
	using InterpCoef = std::pair<uint32_t, double>;
	using InterpCoefs = std::unordered_map<uint32_t, double>;
	using Label = UINT;
	using Labels = std::vector<Label>;
	using Vector3D = BoundaryMesh::Vector3D;
	using Vector3DOp = std::array<InterpCoefs, 3>;
	using Matrix3D = math::matrix_c<double, 3, 3>;
	using Matrix2D = math::matrix_c<double, 2, 2>;
	using Node = EvaporatingParticle::CNode3D;
	using Element = EvaporatingParticle::CElem3D;
	using ScalarFieldOperator = CFieldOperator<double>;
	using PScalFieldOp = std::unique_ptr<ScalarFieldOperator>;
	using ProgressBar = EvaporatingParticle::CObject;
	using PProgressBar = std::unique_ptr<ProgressBar>;
	using Graph = CMeshConnectivity;
	using PGraph = std::unique_ptr<Graph>;

private:
	const Elements& m_elems;
	const Nodes& m_nodes;
	PBoundary m_pBoundary;

	//Lazy because it will be created only once by demand
	PGraph m_lazyGraph;

	//Progress bar interface
	PProgressBar m_pProgressBar;

	//Factor for small step calculation depending on position inside mesh
	double m_fSmallStepFactor;
	double m_fEpsilon;

	//Looks for an element containing point with coordinates v
	const Element* element(const Vector3D& v, Label& nCurNode, const Label& nPrevNode = Label(0)) const;

	//Math operations with interpolation coeffs
	static InterpCoefs& add(InterpCoefs& ic1, const InterpCoefs& ic2);
	static InterpCoefs& mul(double h, InterpCoefs& ic);
	
	//Creates covariance matrix for dirrections around label 
	Matrix3D covarianceOfDirections(Label l) const;
	//Returns vector of finite difference total projections on a directions x,y,z
	Vector3DOp finDiffDirCov(Label l) const;
	//Returns vector of coefficients for grad in point calculation using directed derivatives averaging
	InterpCoefs gradX(Label l) const;
	InterpCoefs gradY(Label l) const;
	InterpCoefs gradZ(Label l) const;

	//Creates mesh graph connections
	void lazyGraphCreation();

	//Create rough and fast LaplacianField solver DU = 0 for a zero approximation
	ScalarFieldOperator laplacianSolver0();
	//Simple Laplacian solver with step myltiplication by a factor
	ScalarFieldOperator laplacianSolver1();
	//Directed derivative calculation
	ScalarFieldOperator directedDerivative(const Vector3D& dir);

	//Obtains interpolating coefficients if containing element is known
	InterpCoefs interpCoefs(const Vector3D& pos, const Element* e) const;

	//Looks for space position pos in closest neighbor elements
	const Element* lookInClosestElements(const Vector3D& pos, Label l) const;
public:
	CMeshAdapter(const Elements& es, const Nodes& ns, double fSmallStepFactor = 0.3);

	//Access to a progress bar interface
	ProgressBar* progressBar() const;

	//Release mesh graph and free memory
	void releaseGraph();

	//Gets and sets small step factor value
	double smallStepFactor() const;
	void smallStepFactor(double fVal);
	//Gets and sets epsilon
	double eps() const;
	void eps(size_t nFactor);

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
		GradX,
		GradY,
		GradZ
	};
	PScalFieldOp createOperator(ScalarOperatorType type = LaplacianSolver1);

	//Returns operators calculating gradient components
	ScalarFieldOperator gradX();
	ScalarFieldOperator gradY();
	ScalarFieldOperator gradZ();
};

#endif // !_MESH_DATA_
