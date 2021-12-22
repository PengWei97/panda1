//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementAverageMaterialPropertyOP.h"


registerMooseObject("MooseApp", ElementAverageMaterialPropertyOP);
registerMooseObject("MooseApp", ADElementAverageMaterialPropertyOP);

template <bool is_ad>
InputParameters
ElementAverageMaterialPropertyOPTempl<is_ad>::validParams()
{
  InputParameters params = ElementIntegralMaterialPropertyTempl<is_ad>::validParams();
  params.addClassDescription("Computes the average of a material property over a volume.");
  params.addRequiredCoupledVar(
      "second_variable",
      "The name of the second variable in the inner product (variable, second_variable)");
  return params;
}

template <bool is_ad>
ElementAverageMaterialPropertyOPTempl<is_ad>::ElementAverageMaterialPropertyOPTempl(
    const InputParameters & parameters)
  : ElementIntegralMaterialPropertyTempl<is_ad>(parameters), _volume(0.0), _v(coupledValue("second_variable"))
{
}

// _volume(0.0) 为整个域的体积值;

template <bool is_ad>
void
ElementAverageMaterialPropertyOPTempl<is_ad>::initialize()
{
  ElementIntegralMaterialPropertyTempl<is_ad>::initialize();
  // _integral_value = 0;

  _volume = 0.0;
}

template <bool is_ad>
void
ElementAverageMaterialPropertyOPTempl<is_ad>::execute()
{
  ElementIntegralMaterialPropertyTempl<is_ad>::execute();
  // _integral_value += computeIntegral(); // 积分所计算的值

      
}

template <bool is_ad>
Real
ElementAverageMaterialPropertyOPTempl<is_ad>::getValue()
{
  const Real integral = ElementIntegralMaterialPropertyTempl<is_ad>::getValue();
  // gatherSum(_integral_value);

  ElementIntegralMaterialPropertyTempl<is_ad>::gatherSum(_volume);
  // gatherSum:: a gather operation is required to collect the values computed on all processes to the root process;

  std::cout << "average = " << integral / _volume << std::endl;
  
  return integral / _volume;
}

template <bool is_ad>
void
ElementAverageMaterialPropertyOPTempl<is_ad>::threadJoin(const UserObject & y)
{
  ElementIntegralMaterialPropertyTempl<is_ad>::threadJoin(y);

  const ElementAverageMaterialPropertyOPTempl<is_ad> & pps =
      static_cast<const ElementAverageMaterialPropertyOPTempl<is_ad> &>(y);
  _volume += pps._volume;
}

//  This is called after the execution of the postprocessor and is intended to perform aggregation for shared memory parallelism.

template class ElementAverageMaterialPropertyOPTempl<false>;
template class ElementAverageMaterialPropertyOPTempl<true>;
