#!/bin/bash
# mkdir ../include/materials
# mkdir ../src/materials

cp ~/projects/moose/modules/tensor_mechanics/include/materials/ComputeLinearElasticStress.h ~/projects/panda/include/materials/ComputeIncrementalLinearElasticStress.h

cp ~/projects/moose/modules/tensor_mechanics/src/materials/ComputeLinearElasticStress.C ~/projects/panda/src/materials/ComputeIncrementalLinearElasticStress.C

code ~/projects/panda/include/materials/ComputeIncrementalLinearElasticStress.h
code ~/projects/panda/src/materials/ComputeIncrementalLinearElasticStress.C
    