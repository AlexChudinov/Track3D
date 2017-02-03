#pragma once
#ifndef _MESH_IMPLEMENTATION_H_
#define _MESH_IMPLEMENTATION_H_ 1

#include <memory>

#include "..\LSExport.h"
#include "..\mesh_math\mesh_geometry.h"

using graph = data_structs::graph<UINT>;
using vector3f = math::vector_c<double, 3>;
using node_positions = std::vector<vector3f>;
using mesh_geom = mesh_geometry<double, UINT>;

class MeshImplementation : public Mesh
{
	using basic_mesh_geometry = mesh_geom;
	std::shared_ptr<mesh_geom> _geometry;
public:
	MeshImplementation(const graph& g, const node_positions& np);

	std::shared_ptr<mesh_geom> geometryPtr();

	std::pair<V3D, V3D> getBox() const;
};

#endif // !_MESH_IMPLEMENTATION_H_
