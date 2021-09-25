//  Creation time: 2021.07.12

#pragma once

#include "Material.h"
#include "RankTwoTensorForward.h"
#include "RankFourTensorForward.h"

class GGElasticEnergyMaterial : public Material
{
public:
  static InputParameters validParams();

  GGElasticEnergyMaterial(const InputParameters & parameters);

  // virtual void initialSetup() override;
protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
   
   RankFourTensor newGrain(unsigned int new_grain_id);

  //  void computeGGElasticEnergy();

   const std::string _base_name;

   const MaterialProperty<RankTwoTensor> & _stress;
   const MaterialProperty<RankFourTensor> & _elasticity_tensor;
   const MaterialProperty<RankTwoTensor> & _strain;

   MaterialProperty<Real> & _GG_elastic_energy;
};

// MaterialProperty<Real> & _elastic_energy;