#!/bin/bash
# mkdir ../include/materials
# mkdir ../src/materials

cp ~/projects/moose/modules/tensor_mechanics/include/materials/FiniteStrainPlasticMaterial.h ~/projects/panda/include/materials/Test1FiniteStrainPlasticMaterial.h

cp ~/projects/moose/modules/tensor_mechanics/src/materials/FiniteStrainPlasticMaterial.C ~/projects/panda/src/materials/Test1FiniteStrainPlasticMaterial.C

code ~/projects/panda/include/materials/ACGrGrElasticDrivingForce.h
code ~/projects/panda/src/materials/ACGrGrElasticDrivingForce.C
    