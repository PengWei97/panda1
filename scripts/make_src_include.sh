name = action
mkdir ../include/name
mkdir ../src/name

cp ~/projects/moose/modules/tensor_mechanics/include/action/PolycrystalElasticDrivingForceAction.h ~/projects/panda/include/action/GGPolycrystalElasticDrivingForceAction.h

cp ~/projects/moose/modules/tensor_mechanics/src/action/PolycrystalElasticDrivingForceAction.C ~/projects/panda/src/action/GGPolycrystalElasticDrivingForceAction.C

code ~/projects/panda/include/action/GGPolycrystalElasticDrivingForceAction.h
code ~/projects/panda/src/action/GGPolycrystalElasticDrivingForceAction.C




