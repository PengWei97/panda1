#!/bin/bash
# mkdir ../include/materials
# mkdir ../src/materials

cp ~/projects/moose/modules/phase_field/include/materials/ElasticEnergyMaterial.h ~/projects/panda/include/materials/PlasticEnergyMaterial.h

cp ~/projects/moose/modules/phase_field/src/materials/ElasticEnergyMaterial.C ~/projects/panda/src/materials/PlasticEnergyMaterial.C

code ~/projects/panda/include/materials/PlasticEnergyMaterial.h
code ~/projects/panda/src/materials/PlasticEnergyMaterial.C
    