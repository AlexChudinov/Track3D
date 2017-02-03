#ifndef MESH_GEOMETRY_H
#define MESH_GEOMETRY_H

#include <map>

#include <linearAlgebra\vectorTemplate.h>
#include <linearAlgebra\linearInterpolation.h>
#include <data_structs\graph.h>

/**
 * Mesh connectivity and node space positions
 */
template<typename Float, typename label>
class mesh_geometry
{
public:
	using string         = std::string;
	using node_labels    = std::vector<label>;
	using vector3f       = math::vector_c<Float, 3>;
	using node_positions = std::vector<vector3f>;
	using graph          = data_structs::graph<label>;
    using box3D          = std::pair<vector3f, vector3f>;
    using label_list	 = std::set<label>;

	//Interpolation coefs
	using InterpCoef  = std::pair<label, Float>;
	using InterpCoefs = std::map<label, Float>;

	//Boundary conditions for the mesh
	class BoundaryMesh
	{
		const mesh_geometry& m_mesh;

		using StrRef = std::reference_wrapper<const std::string>;

		//Declare reference wrappers less operator
		class StrRefLess : public std::binary_function<const StrRef&, const StrRef&, bool>
		{
		public:
			bool operator()(const StrRef& s1, const StrRef& s2) const
			{
				return s1.get() < s2.get();
			}
		};

	public:
		enum BoundaryType { ZERO_GRAD, FIXED_VAL };

		using BoundaryDescription = std::pair<BoundaryType, label_list>;
		using BoundariesMap = std::map<std::string, BoundaryDescription>;
		using NamesList = std::set<StrRef, StrRefLess> ;
		using ReversedBoundariesMap = std::map<label, std::pair<vector3f, NamesList>>;
		using BoundaryNormals = std::map<label, vector3f>;

		using iterator = typename ReversedBoundariesMap::iterator;
		using const_iterator = typename ReversedBoundariesMap::const_iterator;

	private:
		BoundariesMap m_mapBoundariesList;
		ReversedBoundariesMap m_mapReversedBoundariesList;

	public:
		//Creates empty boundary mesh
		BoundaryMesh(const mesh_geometry& mesh):m_mesh(mesh){}

		//Adds new boundary patch
		void addBoundary(
			const std::string& strName, 
			const std::vector<label>& vLabels,
			const std::vector<vector3f>& vNormals,
			BoundaryType type = FIXED_VAL)
		{
			if (vLabels.size() != vNormals.size())
				throw std::runtime_error("BoundaryMesh::addBoundary: Sizes of normals and labels vectors are different.");
			if (*std::max_element(vLabels.begin(), vLabels.end()) >= m_mesh.size())
				throw std::runtime_error("BoundaryMesh::addBoundary: Too big label for used mesh.");
			if (isBoundary(strName)) removeBoundary(strName);
			m_mapBoundariesList[strName] = std::make_pair(type, label_list(vLabels.begin(), vLabels.end()));
			typename BoundariesMap::const_iterator it = m_mapBoundariesList.lower_bound(strName);
			for (size_t i = 0; i < vLabels.size(); ++i)
			{
				std::pair<vector3f, NamesList>& entry = m_mapReversedBoundariesList[vLabels[i]];
				entry.first = vNormals[i];
				entry.second.insert(std::cref(it->first));
			}
		}

		//Removes existing boundary patch
		void removeBoundary(const std::string& strName)
		{
			for (label l : m_mapBoundariesList.at(strName).second)
			{
				m_mapReversedBoundariesList[l].second.erase(strName);
				if (m_mapReversedBoundariesList[l].second.empty()) m_mapReversedBoundariesList.erase(l);
			}
			m_mapBoundariesList.erase(strName);
		}

		//Sets type of a boundary with a name strName
		void boundaryType(const std::string& strName, BoundaryType type)
		{
			m_mapBoundariesList.at(strName).first = type;
		}

		//Returns type of a boundary with a name strName
		BoundaryType boundaryType(const std::string& strName) const
		{
			return m_mapBoundariesList.at(strName).first;
		}

		//Returns a set of boundaries connected to a given label
		const NamesList& boundaryNames(label l) const
		{
			return m_mapReversedBoundariesList.at(l).second;
		}

		//Returns numbers of nodes which are belong to a given boundary
		const label_list& boundaryLabels(const std::string& strName) const
		{
			return m_mapBoundariesList.at(strName).second;
		}

		//Get iterators for boundary
		const_iterator begin() const { return m_mapReversedBoundariesList.begin(); }
 		iterator begin() { return m_mapReversedBoundariesList.begin(); }
		const_iterator end() const { return m_mapReversedBoundariesList.end(); }
		iterator end() { return m_mapReversedBoundariesList.end(); }

		//Checks if the name is in the boundaries list
		bool isBoundary(const std::string& sName) const 
		{ 
			return m_mapBoundariesList.find(sName) != m_mapBoundariesList.end();
		}
		//Checks if the label belongs to the boundary
		bool isBoundary(label l) const
		{
			return m_mapReversedBoundariesList.find(l) != m_mapReversedBoundariesList.end();
		}
		//Checks if it is a first-type (Dirichlet) boundary condition
		bool isFirstType(const std::string& sName) const
		{
			return m_mapBoundariesList.at(sName).first == FIXED_VAL;
		}
		bool isFirstType(label l) const
		{
			return std::accumulate(m_mapReversedBoundariesList.at(l).second.begin(),
				m_mapReversedBoundariesList.at(l).second.end(), false, 
				[=](bool val, const std::string& sName)->bool
			{
				return val |= isFirstType(sName);
			});
		}

		//Gets the normal for given node label
		vector3f normal(label l) const
		{
			return m_mapReversedBoundariesList.at(l).first;
		}
	};

private:
    graph mesh_connectivity_;
    node_positions node_positions_;

	//Numeric limit for floating point precision
	Float m_fEpsilon;
public:
	mesh_geometry(const graph& g, const node_positions& np)
        : mesh_connectivity_(g), node_positions_(np), m_fEpsilon(std::numeric_limits<Float>::epsilon()*100.0)
    { 
		if(g.size() != np.size()) 
			throw(std::runtime_error("Sizes of graph and node positions array mismatch!"));
	}

	//Sets the precision limit
	void eps(size_t fFactor)
	{
		m_fEpsilon = std::numeric_limits<Float>::epsilon()*static_cast<Float>(fFactor);
	}
	//Gets the precision limit
	Float eps() const
	{
		return m_fEpsilon;
	}

	//Creates boundary conditions for given mesh
	BoundaryMesh createBoundary() const
	{
		return BoundaryMesh(*this);
	}

    /**
     * Returns all mesh points in the order where points at odd indexes coincides with points at even indexes,
	 * Therefore, if we put the return to a glDrawArray function with GL_LINES it will draw the mesh
     */
    node_positions meshConnections() const
    {
        node_positions res(mesh_connectivity_.connectionsNum()*2); //Each connection contains two vertices
        size_t counter = 0;
        auto observer = [&](size_t i, size_t j)
        {
            res[counter++] = node_positions_[i];
            res[counter++] = node_positions_[j];
        };
        mesh_connectivity_.iterateOverUniqueConnections(observer);
        return res;
    }

    /**
     * Returns a number of nodes in a mesh
     */
    inline size_t size() const { return node_positions_.size(); }

    /**
     * Returns minimal box that contains whole mesh
     */
    box3D box() const
    {
        Float   min_x = node_positions_[0][0],
                min_y = node_positions_[0][1],
                min_z = node_positions_[0][2],
                max_x = node_positions_[0][0],
                max_y = node_positions_[0][1],
                max_z = node_positions_[0][2];

        typename node_positions::const_iterator it = node_positions_.begin() + 1;
        for(; it != node_positions_.end(); ++it)
        {
            min_x = std::min(min_x, (*it)[0]);
            min_y = std::min(min_y, (*it)[1]);
            min_z = std::min(min_z, (*it)[2]);
            max_x = std::max(max_x, (*it)[0]);
            max_y = std::max(max_y, (*it)[1]);
            max_z = std::max(max_z, (*it)[2]);
        }
        return box3D({min_x, min_y, min_z}, {max_x, max_y, max_z});
    }

	/**
	 * Returns point3D by a label
	 */
	inline const vector3f& spacePositionOf(label id) const { return node_positions_[id]; }

	/**
	 * Returns shortest edge length incident to the given label id
	 */
	inline double shortestEdgeLength(label id) const
	{
		return math::abs(spacePositionOf(
			*std::min_element(
			mesh_connectivity_.getNeighbour(id).begin(),
			mesh_connectivity_.getNeighbour(id).end(),
			[&](label l1, label l2)->bool
		{
			return math::sqr(spacePositionOf(l1) - spacePositionOf(id))
				< math::sqr(spacePositionOf(l2) - spacePositionOf(id));
		})) - spacePositionOf(id));
	}

	/**
	 * Gets nodes incident to a node with the label id
	 */
	inline const label_list& neighbour(label id) const { return mesh_connectivity_.getNeighbour(id); }

	/**
	 * Visit neigbour points
	 */
	template<typename visitor>
	void visit_neigbour(label l, visitor V) const
	{
		for (label ll : mesh_connectivity_.getNeighbour(l)) V(ll);
	}

	/**
	 * Search for a clossest point from the start point
	 * Returns the label of that point
	 */
	label find_closest(Float x, Float y, Float z, label start = 0) const
	{
		const vector3f pos{ x, y, z };
		label result = start;
		double minSqrDist = math::sqr(node_positions_[start] - pos);

		if (minSqrDist == 0.0) return start;

		mesh_connectivity_.bfs_iterative(start, 
			[&](label l)->bool 
		{
			double testSqrDist = math::sqr(node_positions_[l] - pos);
			if (testSqrDist <= minSqrDist)
			{
				minSqrDist = testSqrDist;
				result = l;
				return true;
			}
			return false;
		});

		return result;
	}

	//Find closest line if one point is known
	//Assume, that it was found using function find_closest
	label find_line(Float x, Float y, Float z, label start) const
	{
		const vector3f pos = vector3f{ x,y,z } - node_positions_[start];
		const label_list& neighbor = mesh_connectivity_.getNeighbour(start);

		return *std::max_element(neighbor.begin(), neighbor.end(),
			[&](label l1, label l2)->bool
		{
			vector3f 
				v1 = node_positions_[l1] - node_positions_[start],
				v2 = node_positions_[l2] - node_positions_[start];
			return pos*v1 / math::abs(v1) < pos*v2 / math::abs(v2);
		});
	}

	//Find closest plane
	//It is assumed that next and start already correspond to a closest line
	label find_plane(Float x, Float y, Float z, label start, label next) const
	{
		vector3f
			pos = vector3f{ x,y,z } -node_positions_[start],
			e0 = node_positions_[next] - node_positions_[start];
		pos -= (pos*e0)*e0 / math::sqr(e0);

		const label_list& neighbor = mesh_connectivity_.getNeighbour(start);

		return *std::max_element(neighbor.begin(), neighbor.end(),
			[&](label l1, label l2)->bool
		{
			if (l1 == next) return true;
			if (l2 == next) return false;
			vector3f
				v1 = node_positions_[l1] - node_positions_[start],
				v2 = node_positions_[l2] - node_positions_[start];
			return pos*v1 / math::abs(v1) < pos*v2 / math::abs(v2);
		});
	}

	//Find best tetrahedral for the x,y,z point
	//First three nodes of the tet were previously found by the find_closest, 
	//find_line, and find_plane functions
	label find_tet(Float x, Float y, Float z, label start, label next1, label next2) const
	{
		vector3f 
			pos = vector3f{ x,y,z } - node_positions_[start],
			e0 = node_positions_[next1] - node_positions_[start],
			e1 = node_positions_[next2] - node_positions_[start];

		const label_list& neighbor = mesh_connectivity_.getNeighbour(start);

		return *std::max_element(neighbor.begin(), neighbor.end(),
			[&](label l1, label l2)->bool
		{
			vector3f
				v1 = node_positions_[l1] - node_positions_[start],
				v2 = node_positions_[l2] - node_positions_[start];
			double 
				norm1 = ::fabs(math::det(math::matrix_c<double, 3, 3>{v1, e0, e1})),
				norm2 = ::fabs(math::det(math::matrix_c<double, 3, 3>{v2, e0, e1}));
			if (norm1 < eps() && norm2 >= eps()) return true;
			if (norm2 < eps()) return false;
			return pos*v1 / math::abs(v1) < pos*v2 / math::abs(v2);
		});
	}

	/**
	 * Returns coeffs for field interpolation
	 */
	InterpCoefs interpCoefs(Float x, Float y, Float z, label start = 0) const
	{
		label l0, l1, l2, l3;
		vector3f pos{ x,y,z };
		l0 = find_closest(x, y, z, start);
		vector3f dp0 = pos - node_positions_[l0];
		if (math::sqr(dp0) < eps())
		{
			return InterpCoefs{ InterpCoef(l0, 1.0) };
		}
		else
		{
			l1 = find_line(x, y, z, l0);
			vector3f e0 = node_positions_[l1] - node_positions_[l0],
				dp1 = dp0 - (dp0 * e0) * e0 / math::sqr(e0);
			if (math::sqr(dp1) < eps())
			{
				std::tuple<Float, Float> coefs
					= math::lineInterpolation(pos, node_positions_[l0], node_positions_[l1]);
				return InterpCoefs
				{ 
					InterpCoef(l0, std::get<0>(coefs)), 
					InterpCoef(l1, std::get<1>(coefs)) 
				};
			}
			else
			{
				l2 = find_plane(x, y, z, l0, l1);
				vector3f 
					e1 = node_positions_[l2] - node_positions_[l0],
					dp2 = dp1 - (dp1 * e1) * e1 / math::sqr(e1);
				if (math::sqr(dp2) < eps())
				{
					std::tuple<Float, Float, Float> coefs
						= math::triInterpolation(pos, node_positions_[l0], node_positions_[l1], node_positions_[l2]);
					return InterpCoefs
					{
						InterpCoef(l0, std::get<0>(coefs)),
						InterpCoef(l1, std::get<1>(coefs)),
						InterpCoef(l2, std::get<2>(coefs))
					};
				}
				else
				{
					l3 = find_tet(x, y, z, l0, l1, l2);
					std::tuple<Float, Float, Float, Float> coefs
						= math::tetInterpolation(pos, node_positions_[l0], node_positions_[l1], node_positions_[l2], node_positions_[l3]);
					return InterpCoefs
					{
						InterpCoef(l0, std::get<0>(coefs)),
						InterpCoef(l1, std::get<1>(coefs)),
						InterpCoef(l2, std::get<2>(coefs)),
						InterpCoef(l3, std::get<3>(coefs))
					};
				}
			}
		}
	}
};

#endif // MESH_GEOMETRY_H
