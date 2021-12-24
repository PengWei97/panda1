#!/bin/bash
# mkdir ../include/postprocessors
# mkdir ../src/postprocessors

cp ~/projects/moose/framework/include/postprocessors/ElementAverageMaterialProperty.h ~/projects/panda/include/postprocessors/ElementAverageMaterialPropertyOP.h

cp ~/projects/moose/framework/src/postprocessors/ElementAverageMaterialProperty.C ~/projects/panda/src/postprocessors/ElementAverageMaterialPropertyOP.C

code ~/projects/panda/include/postprocessors/ElementAverageMaterialPropertyOP.h
code ~/projects/panda/src/postprocessors/ElementAverageMaterialPropertyOP.C
    


Dear MOOSE experts,

Recently, I am trying to find the element average or node average of a `MaterialProperty ` in the material module. I focused on referring to the `ElementAverageMaterialProperty` class in the processing module. `ElementAverageMaterialPropertyTest` was created to implement such operations, the code is as follows,
```c++
// ElementAverageMaterialPropertyTest.h
protected:


```

```c++
// ElementAverageMaterialPropertyTest.C

```
But an error occurred during the calculation process, the error is as follows,


Any suggestions or recommendations to fix the problem would be greatly appreciated.

Thank you
Wei