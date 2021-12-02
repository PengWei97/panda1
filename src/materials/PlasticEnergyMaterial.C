//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PlasticEnergyMaterial.h"

registerMooseObject("PhaseFieldApp", PlasticEnergyMaterial);

InputParameters
PlasticEnergyMaterial::validParams()
{
  InputParameters params = DerivativeFunctionMaterialBase::validParams();
  params.addClassDescription("Free energy material for the elastic energy contributions.");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredCoupledVar("args", "Arguments of F() - use vector coupling");
  return params;
}

PlasticEnergyMaterial::PlasticEnergyMaterial(const InputParameters & parameters)
  : DerivativeFunctionMaterialBase(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _eqv_plastic_strain(getMaterialPropertyByName<Real>(_base_name + "eqv_plastic_strain")),
    _hard_factor(getMaterialPropertyByName<Real>(_base_name + "hard_factor"))
    // _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "elastic_strain"))
{
  _D_hard_factor.resize(_nargs);
  // _d2elasticity_tensor.resize(_nargs);

  // fetch stress and elasticity tensor derivatives (in simple eigenstrain models this is is only
  // w.r.t. 'c')
  for (unsigned int i = 0; i < _nargs; ++i)
  {
    _D_hard_factor[i] = &getMaterialPropertyDerivativeByName<Real>(_base_name + "hard_factor", _arg_names[i]);
  }
}

void
PlasticEnergyMaterial::initialSetup()
{
  validateCoupling<Real>(_base_name + "eqv_plastic_strain");
  validateCoupling<Real>(_base_name + "hard_factor");
}

Real
PlasticEnergyMaterial::computeF()
{
  return 0.5 * _hard_factor[_qp]*_eqv_plastic_strain[_qp]*_eqv_plastic_strain[_qp];
}

Real
PlasticEnergyMaterial::computeDF(unsigned int i_var)
{
  unsigned int i = argIndex(i_var);

  // product rule d/di computeF (doubleContraction commutes)
  return 0.5 * (*_D_hard_factor[i])[_qp]*_eqv_plastic_strain[_qp]*_eqv_plastic_strain[_qp];
}
