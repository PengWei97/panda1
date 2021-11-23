//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TensorMechanicsHardeningModel.h"
#include "GrainDataTracker.h"

/**
 * Power-Rule Hardening defined by:
 * assuming p = internal_parameter, then value = value_0 * (p / epsilon0 + 1)^{exponent})
 * Notice that if epsilon0 = 0, it will return not a number.
 */
class Test1TensorMechanicsHardeningLinear : public TensorMechanicsHardeningModel
{
public:
  static InputParameters validParams();

  Test1TensorMechanicsHardeningLinear(const InputParameters & parameters);

  virtual Real value(Real intnl) const override;

  virtual Real derivative(Real intnl) const override;

  virtual std::string modelName() const override;

private:
  /// The value = value_0 + _h_factor*p
  const Real _value_0;

  /// The value = value_0 + _h_factor*p
  // const Real _hard_factor;
  
  // /// Number of order parameters
  // const unsigned int _op_num;

  // /// Order parameters
  // const std::vector<const VariableValue *> _vals;

  // const MaterialProperty<Real> & _sigma_y_init;

  const MaterialProperty<Real> & _hard_factor;
};
