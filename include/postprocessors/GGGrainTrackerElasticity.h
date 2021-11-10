//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GrainDataTracker.h"
#include "RankFourTensor.h"

class GGEulerAngleProvider;

/**
 * Manage a list of elasticity tensors for the grains
 */
class GGGrainTrackerElasticity : public GrainDataTracker<RankFourTensor>
{
public:
  static InputParameters validParams();

  GGGrainTrackerElasticity(const InputParameters & parameters);

protected:
  RankFourTensor newGrain(unsigned int new_grain_id);

  Real hardFactor(unsigned int new_grain_id);

  /// generate random rotations when the Euler Angle provider runs out of data (otherwise error out)
  const bool _random_rotations;

  /// unrotated elasticity tensor
  RankFourTensor _C_ijkl;

  /// object providing the Euler angles
  const GGEulerAngleProvider & _euler;
};
