import numpy as np

def writesh(dirname, filename,editFileName, threads, file, outputshfilename):
    shfile = """#!/bin/bash
# mkdir ../include/%s
# mkdir ../src/%s

cp /home/pengwei/projects/panda/include/%s/%s.h ~/projects/panda/include/%s/%s.h

cp /home/pengwei/projects/panda/src/%s/%s.C ~/projects/panda/src/%s/%s.C

code ~/projects/panda/include/%s/%s.h
code ~/projects/panda/src/%s/%s.C
    """ % (dirname,dirname,dirname,filename,dirname,editFileName,dirname,filename,dirname,editFileName,dirname,editFileName,dirname,editFileName)
    with open(outputshfilename, "w") as f:
        f.write(shfile)
    f.close()

writesh( "materials", "TestFiniteStrainPlasticMaterial","Test2FiniteStrainPlasticMaterial", 1, 1, "2.sh")

# /home/pengwei/projects/moose/modules/tensor_mechanics/include/userobjects/EulerAngleProvider.h
# /home/pengwei/projects/moose/modules//src/materials/FiniteStrainPlasticMaterial.C

# /home/pengwei/projects/panda/src/materials/TestFiniteStrainPlasticMaterial.C