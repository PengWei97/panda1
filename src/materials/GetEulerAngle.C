//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GetEulerAngles.h"
#include "GrainTracker.h"
#include "EulerAngleProvider.h"
#include "GrainTrackerInterface.h"

registerMooseObject("PhaseFieldApp", GetEulerAngles);

InputParameters
GetEulerAngles::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Get Euler angles from user object to an Material property.");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  params.addRequiredParam<UserObjectName>("grain_tracker",
                                          "The GrainTracker UserObject to get values from.");
  MooseEnum euler_angles("phi1 Phi phi2");
  params.addRequiredParam<MooseEnum>("output_euler_angle", euler_angles, "Euler angle to get");
  return params;
}

GetEulerAngles::GetEulerAngles(const InputParameters & parameters)
  : Material(parameters),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider")),
    _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
    _output_euler_angle(getParam<MooseEnum>("output_euler_angle")),
    _euler_phi1(declareProperty<Real>("euler_phi1"))
{
}

// void
// GetEulerAngles::precalculateValue()
// {
//   // ID of unique grain at current point
//   const auto grain_id =
//       _grain_tracker.getEntityValue(_current_elem->id(),
//                                     FeatureFloodCount::FieldType::UNIQUE_REGION,
//                                     0);
//   std::cout << "Grian_ID is: " << grain_id << std::endl;

//   // Recover euler angles for current grain
//   RealVectorValue angles;
//   if (grain_id >= 0)
//     angles = _euler.getEulerAngles(grain_id);

//   // Return specific euler angle
//   _value = angles(_output_euler_angle);
//   std::cout << "value is: " << _value << std::endl;
// }

// Real
// OutputEulerAngles::computeValue()
// {
//   return _value;
// }

void
GetEulerAngles::computeQpProperties()
{
    // ID of unique grain at current point
  const auto grain_id =
      _grain_tracker.getEntityValue(_current_elem->id(),
                                    FeatureFloodCount::FieldType::UNIQUE_REGION,
                                    0);
  std::cout << "Grian_ID is: " << grain_id << std::endl;

  // Recover euler angles for current grain
  RealVectorValue angles;
  if (grain_id >= 0)
    angles = _euler.getEulerAngles(grain_id);

  // Return specific euler angle
  _value = angles(_output_euler_angle);
  std::cout << "value is: " << _value << std::endl;

  _euler_phi1[_qp] = _value;
}
