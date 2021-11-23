//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Test1TensorMechanicsHardeningLinear.h"
#include <vector>

registerMooseObject("TensorMechanicsApp", Test1TensorMechanicsHardeningLinear);

InputParameters
Test1TensorMechanicsHardeningLinear::validParams()
{
  InputParameters params = TensorMechanicsHardeningModel::validParams();
  // params.addRequiredParam<Real>("value_0", "The yield strength when internal variable = 0");
  // params.addRequiredParam<Real>("HardFactor", "The hardening factor");
  params.addClassDescription("Hardening defined by linear hardenig rule");
  // params.addRequiredCoupledVarWithAutoBuild(
  //     "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

Test1TensorMechanicsHardeningLinear::Test1TensorMechanicsHardeningLinear(
    const InputParameters & parameters)
  : TensorMechanicsHardeningModel(parameters),
    _value_0(getParam<Real>("value_0")),
    // _sigma_y_init(getMaterialProperty<Real>("sigma_y_init")),
    // _hard_factor(getParam<Real>("HardFactor"))
    _hard_factor(getMaterialProperty<Real>("HardFactor"))
    // _op_num(coupledComponents("v")),
    // _vals(coupledValues("v"))
{
}

Real
Test1TensorMechanicsHardeningLinear::value(Real intnl) const
{
  return _value_0 + _hard_factor[_qp] * intnl;
  // return _value_0 * std::pow(intnl / _epsilon0 + 1, _exponent);
}

Real
Test1TensorMechanicsHardeningLinear::derivative(Real intnl) const // d*/d_intnl
{
  return _hard_factor;
}

std::string
Test1TensorMechanicsHardeningLinear::modelName() const
{
  return "LinearHardening";
}
