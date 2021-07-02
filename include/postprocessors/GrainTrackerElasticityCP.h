//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GrainDataTrackerCP.h"
#include "RankFourTensor.h"

class EulerAngleProvider;

/**
 * Manage a list of elasticity tensors for the grains
 */
class GrainTrackerElasticityCP : public GrainDataTrackerCP<RankFourTensor>
{
public:
  static InputParameters validParams();

  GrainTrackerElasticityCP(const InputParameters & parameters);

protected:
  RankFourTensor newGrain(unsigned int new_grain_id);
  // 

  /// generate random rotations when the Euler Angle provider runs out of data (otherwise error out)
  const bool _random_rotations;
  // 是否提供随机旋转矩阵

  /// unrotated elasticity tensor
  RankFourTensor _C_ijkl;
  // 没有旋转的弹性模量

  /// object providing the Euler angles
  const EulerAngleProvider & _euler;
  // 所提供的欧拉角
};
