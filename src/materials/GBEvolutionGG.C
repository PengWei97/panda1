//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBEvolutionGG.h"

registerMooseObject("PhaseFieldApp", GBEvolutionGG);
registerMooseObject("PhaseFieldApp", ADGBEvolutionGG);

template <bool is_ad>
InputParameters
GBEvolutionGGTempl<is_ad>::validParams()
{
  InputParameters params = GBEvolutionGGBaseTempl<is_ad>::validParams();
  params.addRequiredParam<Real>("GBenergy", "Grain boundary energy in J/m^2");
  return params;
}

template <bool is_ad>
GBEvolutionGGTempl<is_ad>::GBEvolutionGGTempl(const InputParameters & parameters)
  : GBEvolutionGGBaseTempl<is_ad>(parameters), _GBEnergy(this->template getParam<Real>("GBenergy"))
{
}

template <bool is_ad>
void
GBEvolutionGGTempl<is_ad>::computeQpProperties()
{
  // eV/nm^2
  this->_sigma[this->_qp] = _GBEnergy * this->_JtoeV * (this->_length_scale * this->_length_scale);

  GBEvolutionGGBaseTempl<is_ad>::computeQpProperties();
}

template class GBEvolutionGGTempl<false>;
template class GBEvolutionGGTempl<true>;
