//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementIntegralVariablePostprocessorOP.h"

registerMooseObject("MooseApp", ElementIntegralVariablePostprocessorOP);

defineLegacyParams(ElementIntegralVariablePostprocessorOP);

InputParameters
ElementIntegralVariablePostprocessorOP::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  params.addClassDescription("Computes a volume integral of the specified variable");
  params.addRequiredCoupledVar(
      "second_variable",
      "The name of the second variable in the inner product (variable, second_variable)");
  return params;
}

ElementIntegralVariablePostprocessorOP::ElementIntegralVariablePostprocessorOP(
    const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable"))
    // _v(coupledValue("second_variable"))
{
  addMooseVariableDependency(&mooseVariableField());
}

Real
ElementIntegralVariablePostprocessorOP::computeQpIntegral()
{
  std::cout << "num = " <<  _u[_qp] << std::endl;
  return _u[_qp];
}
