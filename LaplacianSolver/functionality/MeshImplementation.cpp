#include "MeshImplementation.h"

MeshImplementation::MeshImplementation(const graph& g, const node_positions& np)
	: Mesh(), _geometry(new mesh_geom(g, np))
{}

std::shared_ptr<mesh_geom> MeshImplementation::geometryPtr()
{
	return _geometry;
}

std::pair<V3D, V3D> MeshImplementation::getBox() const
{
	mesh_geom::box3D box_ = _geometry->box();
	V3D min, max;
	min.x = box_.first[0]; max.x = box_.second[0];
	min.y = box_.first[1]; max.y = box_.second[1];
	min.z = box_.first[2]; max.z = box_.second[2];
	return std::make_pair(min, max);
}
