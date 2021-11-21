#!/bin/bash
# mkdir ../include/materials
# mkdir ../src/materials

cp ~/projects/moose/modules/phase_field/include/materials/ElasticEnergyMaterial.h ~/projects/panda/include/materials/GGElasticEnergyMaterial.h

cp ~/projects/moose/modules/phase_field/src/materials/ElasticEnergyMaterial.C ~/projects/panda/src/materials/GGElasticEnergyMaterial.C

code ~/projects/panda/include/materials/GGElasticEnergyMaterial.h
code ~/projects/panda/src/materials/GGElasticEnergyMaterial.C
    