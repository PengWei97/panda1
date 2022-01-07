#!/bin/bash
# mkdir ../include/materials
# mkdir ../src/materials

cp ~/projects/moose/modules/phase_field/include/materials/ComputePolycrystalElasticityTensor.h ~/projects/panda/include/materials/ComputePolycrystalElasticityTensorGG.h

cp ~/projects/moose/modules/phase_field/src/materials/ComputePolycrystalElasticityTensor.C ~/projects/panda/src/materials/ComputePolycrystalElasticityTensorGG.C

code ~/projects/panda/include/materials/ComputePolycrystalElasticityTensorGG.h
code ~/projects/panda/src/materials/ComputePolycrystalElasticityTensorGG.C
    