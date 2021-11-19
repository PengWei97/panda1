#!/bin/bash
# mkdir ../include/materials
# mkdir ../src/materials

cp /home/pengwei/projects/panda/include/materials/Test4FiniteStrainPlasticMaterial.h ~/projects/panda/include/materials/Test5FiniteStrainPlasticMaterial.h

cp /home/pengwei/projects/panda/src/materials/Test4FiniteStrainPlasticMaterial.C ~/projects/panda/src/materials/Test5FiniteStrainPlasticMaterial.C

code ~/projects/panda/include/materials/EnergyMaterial.h
code ~/projects/panda/src/materials/EnergyMaterial.C
    