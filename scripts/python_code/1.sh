#!/bin/bash
# mkdir ../include/materials
# mkdir ../src/materials

cp ~/projects/moose/modules/phase_field/include/materials/DerivativeMultiPhaseMaterial.h ~/projects/panda/include/materials/DerivativeGrainGrowthMaterial.h

cp ~/projects/moose/modules/phase_field/src/materials/DerivativeMultiPhaseMaterial.C ~/projects/panda/src/materials/DerivativeGrainGrowthMaterial.C

code ~/projects/panda/include/materials/DerivativeGrainGrowthMaterial.h
code ~/projects/panda/src/materials/DerivativeGrainGrowthMaterial.C
    