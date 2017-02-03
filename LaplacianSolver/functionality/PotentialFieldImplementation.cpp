#include "PotentialFieldImplementation.h"
#include "MeshImplementation.h"

PotentialFieldImplementation::PotentialFieldImplementation(Mesh* meshGeom)
	: 
	basic_field(dynamic_cast<MeshImplementation*>(meshGeom)->geometryPtr())
{}

const std::vector<double>& PotentialFieldImplementation::getPotentialVals() const
{
	return basic_field::data();
}

void PotentialFieldImplementation::setBoundaryVal(const std::string & name, double val)
{
	basic_field::set_boundary_uniform_val(name, val);
}

void PotentialFieldImplementation::setBoundaryVal(const std::string & name, const std::vector<double>& vals)
{
	basic_field::set_boundary_vals(name, vals);
}

void PotentialFieldImplementation::addBoundary(
	const std::string& sName,
	const std::vector<UINT>& vLabels,
	const std::vector<V3D>& vNormals)
{
	std::vector<vector3f> vNormalsInner(vNormals.size());
	std::transform(vNormals.begin(), vNormals.end(), vNormalsInner.begin(),
		[](V3D v)->vector3f { return vector3f{ v.x, v.y, v.z }; });
	basic_field::add_boundary(sName, vLabels, vNormalsInner);
}

void PotentialFieldImplementation::setBoundaryType(const std::string & name, BOUNDARY_TYPE type)
{
	switch (type)
	{
	case FIXED_VAL: return basic_field::set_boundary_type(name, basic_field::BoundaryMesh::FIXED_VAL);
	case ZERO_GRAD: return basic_field::set_boundary_type(name, basic_field::BoundaryMesh::ZERO_GRAD);
	default: throw std::runtime_error(
		"PotentialFieldImplementation::setBoundaryType :"
		"Unsupported boundary field type");
	}
}

void PotentialFieldImplementation::applyBoundaryConditions()
{
	basic_field::applyBoundaryConditions();
}

void PotentialFieldImplementation::diffuse()
{
	basic_field next = basic_field::diffuse();
	data() = next.data();
}

double PotentialFieldImplementation::interpolate(double x, double y, double z, UINT * track_label) const
{
	return basic_field::interpolate(x, y, z, track_label);
}

