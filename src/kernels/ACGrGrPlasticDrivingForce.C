//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACGrGrPlasticDrivingForce.h"

#include "Material.h"
// #include "RankFourTensor.h"
// #include "RankTwoTensor.h"

registerMooseObject("PhaseFieldApp", ACGrGrPlasticDrivingForce);

InputParameters
ACGrGrPlasticDrivingForce::validParams()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Adds elastic energy contribution to the Allen-Cahn equation");
  params.addRequiredParam<MaterialPropertyName>(
      "D_hard_factor_name", "The elastic tensor derivative for the specific order parameter");
  return params;
}

ACGrGrPlasticDrivingForce::ACGrGrPlasticDrivingForce(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    _D_hard_factor(getMaterialProperty<Real>("D_hard_factor_name")),
    _eqv_plastic_strain(getMaterialPropertyByName<Real>("eqv_plastic_strain"))
{
}

Real
ACGrGrPlasticDrivingForce::computeDFDOP(PFFunctionType type)
{
  // Access the heterogeneous strain calculated by the Solid Mechanics kernels
  // RankTwoTensor strain(_elastic_strain[_qp]);

  // Compute the partial derivative of the stress wrt the order parameter
  // RankTwoTensor D_stress = _D_elastic_tensor[_qp] * strain;

  switch (type)
  {
    case Residual:
      return 0.5 *
             _D_hard_factor[_qp] * _eqv_plastic_strain[_qp] * _eqv_plastic_strain[_qp]; // Compute the deformation energy driving force

    case Jacobian:
      return 0.0;
  }

  mooseError("Invalid type passed in");
}
