#pragma once
#ifndef _FIELD_OPERATOR_IMPLEMENTATION_
#define _FIELD_OPERATOR_IMPLEMENTATION_

#include "../LSExport.h"
#include "../mesh_math/fieldOperator.h"

class FieldOperatorImplementation : public ScalarFieldOperator, public FieldLinearOp<double>
{
	using basic_operator = FieldLinearOp<double>;
public:
	FieldOperatorImplementation(const field<double>& field, ScalarFieldOperator::OperatorType type);

	void applyToField(PotentialField* field) const;
};

#endif //_FIELD_OPERATOR_IMPLEMENTATION_