#!/bin/bash
# mkdir ../include/materials
# mkdir ../src/materials

cp /home/pengwei/projects/panda/src/userobjects/TensorMechanicsHardeningLinear.C /home/pengwei/projects/panda/src/userobjects/Test1TensorMechanicsHardeningLinear.C

cp /home/pengwei/projects/panda/include/userobjects/TensorMechanicsHardeningLinear.h /home/pengwei/projects/panda/include/userobjects/Test1TensorMechanicsHardeningLinear.h

code /home/pengwei/projects/panda/src/userobjects/Test1TensorMechanicsHardeningLinear.C
code /home/pengwei/projects/panda/include/userobjects/Test1TensorMechanicsHardeningLinear.h
    