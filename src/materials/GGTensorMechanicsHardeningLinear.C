#include "GGTensorMechanicsHardeningLinear.h"

registerMooseObject("pandaApp", GGTensorMechanicsHardeningLinear);

InputParameters
GGTensorMechanicsHardeningLinear::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<Real>("yield_strength_init", 2.0e3, "the initial vaule yield strength of the problem, in MPa");
  params.addClassDescription("Computes the permeability of a porous medium made up of packed "
                             "steel spheres of a specified diameter in accordance with Pamuk and "
                             "Ozdemir (2012). This also provides a specified dynamic viscosity of "
                             "the fluid in the medium.");
    return params;
}

GGTensorMechanicsHardeningLinear::GGTensorMechanicsHardeningLinear(const InputParameters & parameters)
  : Material(parameters),
    _yield_strength_init(getParam<Real>("yield_strength_init")),
    _sigma_y_init(declareProperty<Real>("sigma_y_init"))
{

}  

void
GGTensorMechanicsHardeningLinear::computeQpProperties()
{
    _sigma_y_init[_qp] = _yield_strength_init;
}
