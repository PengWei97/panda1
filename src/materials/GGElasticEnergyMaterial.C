#include "GGElasticEnergyMaterial.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

registerMooseObject("pandaApp", GGElasticEnergyMaterial);

InputParameters
GGElasticEnergyMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Free energy material for the elastic energy contributions.");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredCoupledVar("args", "Arguments of F() - use vector coupling");
  return params;
}

GGElasticEnergyMaterial::GGElasticEnergyMaterial(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
    _strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "elastic_strain")),
    _GG_elastic_energy(declareProperty<Real>(_base_name + "GG_elastic_energy"))
{

}

// void
// GGElasticEnergyMaterial::initialSetup()
// {
//   validateCoupling<RankTwoTensor>(_base_name + "elastic_strain");
//   validateCoupling<RankFourTensor>(_base_name + "elasticity_tensor");
// }

void
GGElasticEnergyMaterial::initQpStatefulProperties()
{
  // init the diffusivity property (this will become
  // _diffusivity_old in the first call of computeProperties)
  _GG_elastic_energy[_qp] = 0;
}

void 
GGElasticEnergyMaterial::computeQpProperties()
{
  _GG_elastic_energy[_qp] = 0.5 * _stress[_qp].doubleContraction(_strain[_qp]);
}


