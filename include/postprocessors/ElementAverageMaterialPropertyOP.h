//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralMaterialProperty.h"
#include "MooseVariableInterface.h"

#include "Coupleable.h"

/**
 * Computes the average of a material property over a volume.
 */
template <bool is_ad>
class ElementAverageMaterialPropertyOPTempl : public ElementIntegralMaterialPropertyTempl<is_ad>
{
public:
  static InputParameters validParams();

  ElementAverageMaterialPropertyOPTempl(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;
  virtual void threadJoin(const UserObject & y) override;

protected:
  /// Domain volume
  Real _volume;

  /// Holds the values of second_variable at current quadrature points  
  const VariableValue & _v;
};

typedef ElementAverageMaterialPropertyOPTempl<false> ElementAverageMaterialPropertyOP;
typedef ElementAverageMaterialPropertyOPTempl<true> ADElementAverageMaterialPropertyOP;
