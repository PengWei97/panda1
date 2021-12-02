//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DerivativeFunctionMaterialBase.h"

/**
 * Material class to compute the elastic free energy and its derivatives
 */
class PlasticEnergyMaterial : public DerivativeFunctionMaterialBase
{
public:
  static InputParameters validParams();

  PlasticEnergyMaterial(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual Real computeF() override;
  virtual Real computeDF(unsigned int i_var) override;

  const std::string _base_name;

  /// The equivalent plastic strain 
  const MaterialProperty<Real> & _eqv_plastic_strain;

  ///@{ hard factor 
  const MaterialProperty<Real> & _hard_factor;
  std::vector<const MaterialProperty<Real> *> _D_hard_factor;
};
