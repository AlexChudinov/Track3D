#include "fieldOperatorImplementation.h"

FieldOperatorImplementation::FieldOperatorImplementation(const field<double>& field,
	ScalarFieldOperator::OperatorType type)
	:
	basic_operator(field)
{
	switch (type)
	{
	case ScalarFieldOperator::Identity: basic_operator::setToIdentity(); break;
	case ScalarFieldOperator::LaplacianSolver: basic_operator::laplacianSolver(); break;
	default: throw std::runtime_error("FieldOperatorImplementation::FieldOperatorImplementation:"
										 " Unsupported operator type.");
	}

}

void FieldOperatorImplementation::applyToField(PotentialField * field) const
{
	basic_operator::applyToField(*dynamic_cast<basic_operator::Field*>(field));
}
