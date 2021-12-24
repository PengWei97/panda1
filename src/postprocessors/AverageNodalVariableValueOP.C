//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AverageNodalVariableValueOP.h"
#include "MooseMesh.h"
#include "SubProblem.h"

registerMooseObject("MooseApp", AverageNodalVariableValueOP);

defineLegacyParams(AverageNodalVariableValueOP);

InputParameters
AverageNodalVariableValueOP::validParams()
{
  InputParameters params = NodalVariablePostprocessor::validParams();

  params.addClassDescription("Computes the average value of a field by sampling all nodal "
                             "solutions on the domain or within a subdomain");
  params.addRequiredCoupledVar(
      "second_variable",
      "The name of the second variable in the inner product (variable, second_variable)");
  return params;
}

AverageNodalVariableValueOP::AverageNodalVariableValueOP(const InputParameters & parameters)
  : NodalVariablePostprocessor(parameters), _sum(0), _n(0), _v(coupledValue("second_variable"))
{
}

// doco-init-start
void
AverageNodalVariableValueOP::initialize()
{
  _sum = 0;
  _n = 0;
}
// doco-init-end

// doco-execute-get-start
void
AverageNodalVariableValueOP::execute()
{
  // if (_v[_qp] > 0.5)
  // {
    _sum += _u[_qp];
    _n++;
  // }
}

Real
AverageNodalVariableValueOP::getValue()
{
  return _sum / _n;
}
// doco-execute-get-end

// doco-final-start
void
AverageNodalVariableValueOP::finalize()
{
  gatherSum(_sum);
  gatherSum(_n);
}
// doco-final-end

// doco-thread-start
void
AverageNodalVariableValueOP::threadJoin(const UserObject & y)
{
  const AverageNodalVariableValueOP & pps = static_cast<const AverageNodalVariableValueOP &>(y);
  _sum += pps._sum;
  _n += pps._n;
}
// doco-thread-end
