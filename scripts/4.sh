#!/bin/bash
# mkdir ../include/materials
# mkdir ../src/materials

cp ~/projects/moose/modules/tensor_mechanics/include/materials/FiniteStrainPlasticMaterial.h ~/projects/panda/include/materials/Test5FiniteStrainPlasticMaterial.h

cp ~/projects/moose/modules/tensor_mechanics/src/materials/FiniteStrainPlasticMaterial.C ~/projects/panda/src/materials/Test5FiniteStrainPlasticMaterial.C

code ~/projects/panda/include/materials/Test5FiniteStrainPlasticMaterial.h
code ~/projects/panda/src/materials/Test5FiniteStrainPlasticMaterial.C
    