#pragma once

#include "Material.h"

class GGTensorMechanicsHardeningLinear : public Material
{
public:
    static InputParameters validParams();

    GGTensorMechanicsHardeningLinear(const InputParameters & parameters);
    
protected:
    virtual void computeQpProperties() override;

    const Real & _yield_strength_init;

    MaterialProperty<Real> & _sigma_y_init;
};