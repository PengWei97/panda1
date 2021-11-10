#!/bin/bash
# mkdir ../include/userobjects
# mkdir ../src/userobjects

cp ~/projects/moose/modules/tensor_mechanics/include/userobjects/EulerAngleProvider.h ~/projects/panda/include/userobjects/GGEulerAngleProvider.h

cp ~/projects/moose/modules/tensor_mechanics/src/userobjects/EulerAngleProvider.C ~/projects/panda/src/userobjects/GGEulerAngleProvider.C

code ~/projects/panda/include/userobjects/GGEulerAngleProvider.h
code ~/projects/panda/src/userobjects/GGEulerAngleProvider.C
    