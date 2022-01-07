import numpy as np

def writesh(dirname, filename,editFileName, threads, file, outputshfilename):
    shfile = """#!/bin/bash
# mkdir ../include/%s
# mkdir ../src/%s

cp ~/projects/moose/modules/phase_field/include/%s/%s.h ~/projects/panda/include/%s/%s.h

cp ~/projects/moose/modules/phase_field/src/%s/%s.C ~/projects/panda/src/%s/%s.C

code ~/projects/panda/include/%s/%s.h
code ~/projects/panda/src/%s/%s.C
    """ % (dirname,dirname,dirname,filename,dirname,editFileName,dirname,filename,dirname,editFileName,dirname,editFileName,dirname,editFileName)
    with open(outputshfilename, "w") as f:
        f.write(shfile)
    f.close()

writesh( "materials", "ComputePolycrystalElasticityTensor","ComputePolycrystalElasticityTensorGG", 1, 1, "1.sh")

# E:\Github\moose\moose\modules\phase_field\src\materials\ComputePolycrystalElasticityTensor.C