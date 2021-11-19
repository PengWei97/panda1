#pragma once

#include "Material.h"

class EnergyMaterial;

template <>
InputParameters validParams<EnergyMaterial>();

class EnergyMaterial : public Material
{
public:
  EnergyMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpPlasticEnergyDensity(); //override;

private:
//   const MaterialProperty<Real> & _D_hard_factor;
//   const MaterialProperty<Real> & _eqv_plastic_strain;

//   MaterialProperty<Real> & _d_plastic_strain_energy;

};