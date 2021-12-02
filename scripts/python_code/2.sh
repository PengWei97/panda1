#!/bin/bash
# mkdir ../include/kernels
# mkdir ../src/kernels

cp ~/projects/moose/modules/phase_field/include/kernels/ACGrGrElasticDrivingForce.h ~/projects/panda/include/kernels/ACGrGrPlasticDrivingForce.h

cp ~/projects/moose/modules/phase_field/src/kernels/ACGrGrElasticDrivingForce.C ~/projects/panda/src/kernels/ACGrGrPlasticDrivingForce.C

code ~/projects/panda/include/kernels/ACGrGrPlasticDrivingForce.h
code ~/projects/panda/src/kernels/ACGrGrPlasticDrivingForce.C
    