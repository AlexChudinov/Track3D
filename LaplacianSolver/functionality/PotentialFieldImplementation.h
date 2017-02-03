#pragma once
#ifndef _POTENTIAL_FIELD_IMPLEMENTATION_H_
#define _POTENTIAL_FIELD_IMPLEMENTATION_H_ 1

#include <memory>

#include "..\LSExport.h"
#include "..\mesh_math\Field.h"

class PotentialFieldImplementation : public PotentialField, public field<double>
{
	using mesh_geom = mesh_geometry<double, UINT>;
	using basic_field = field<double>;

public:
	PotentialFieldImplementation(Mesh* meshGeom);

	const std::vector<double>& getPotentialVals() const;

	void setBoundaryVal(const std::string& name, double val);

	void setBoundaryVal(const std::string& name, const std::vector<double>& vals);

	void addBoundary(const std::string& sName, const std::vector<UINT>& vLabels, const std::vector<V3D>& vNormals);

	void setBoundaryType(const std::string& name, BOUNDARY_TYPE type);

	void applyBoundaryConditions();

	void diffuse();

	double interpolate(double x, double y, double z, UINT* track_label) const;
};

#endif // !_POTENTIAL_FIELD_IMPLEMENTATION_H_
