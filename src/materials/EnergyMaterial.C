#include "EnergyMaterial.h"

registerMooseObject("pandaApp", EnergyMaterial);

template <>
InputParameters
validParams<EnergyMaterial>()
{
  InputParameters params = validParams<Material>();

  // Allow the user to specify which independent variable's gradient to use for calculating the
  // convection velocity property:
  return params;
}

EnergyMaterial::EnergyMaterial(const InputParameters & parameters)
  : Material(parameters)  
  // _D_hard_factor(getMaterialProperty<Real>("D_hard_factor_name")),
  // _eqv_plastic_strain(getMaterialPropertyByName<Real>("eqv_plastic_strain"))
{
}

void
EnergyMaterial::computeQpPlasticEnergyDensity()
{

}