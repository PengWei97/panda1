//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GrainTrackerElasticityCP.h"
#include "EulerAngleProvider.h"
#include "RotationTensor.h"

registerMooseObject("pandaApp", GrainTrackerElasticityCP);

InputParameters
GrainTrackerElasticityCP::validParams()
{
  InputParameters params = GrainTracker::validParams();
  params.addParam<bool>("random_rotations",
                        true,
                        "Generate random rotations when the Euler Angle "
                        "provider runs out of data (otherwise error "
                        "out)");
  params.addRequiredParam<std::vector<Real>>("C_ijkl", "Unrotated stiffness tensor");
  params.addParam<MooseEnum>(
      "fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  return params;
}

GrainTrackerElasticityCP::GrainTrackerElasticityCP(const InputParameters & parameters)
  : GrainDataTrackerCP<RankFourTensor>(parameters),
    _random_rotations(getParam<bool>("random_rotations")),
    _C_ijkl(getParam<std::vector<Real>>("C_ijkl"),
            getParam<MooseEnum>("fill_method").getEnum<RankFourTensor::FillMethod>()),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider"))
{
}

RankFourTensor
GrainTrackerElasticityCP::newGrain(unsigned int new_grain_id)
{
  // new_grain_id:输入的参数，
  EulerAngles angles;
  // 定义类型为EulerAngles的对象angles (Real phi1, Phi, phi2)

  if (new_grain_id < _euler.getGrainNum())
    angles = _euler.getEulerAngles(new_grain_id);
    // _euler.getGrainNum():文件中所提供欧拉角的列数
  else
  {
    if (_random_rotations)
      angles.random();
    else
      mooseError("GrainTrackerElasticityCP has run out of grain rotation data.");
  }

  RankFourTensor C_ijkl = _C_ijkl;
  C_ijkl.rotate(RotationTensor(RealVectorValue(angles)));

  return C_ijkl;
}
