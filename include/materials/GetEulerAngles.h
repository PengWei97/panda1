//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "EulerAngleProvider.h"
#include "GrainTracker.h"

// Forward Declarations
class GrainTracker;
class EulerAngleProvider;

/**
 * Output euler angles from user object to an AuxVariable.
 */
class GetEulerAngles : public Material
{
public:
  GetEulerAngles(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;
  // virtual void precalculateValue();


  /// Object providing the Euler angles
  const EulerAngleProvider & _euler;

  /// Grain tracker object
  const GrainTracker & _grain_tracker;

  /// Number of grains
  MooseEnum _output_euler_angle;

  /// precalculated element value
  Real _value;

  MaterialProperty<Real> & _euler_phi1;
};
