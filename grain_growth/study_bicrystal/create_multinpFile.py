#!/usr/bin/env python3
#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

import os
import unittest
import pyhit

class TestExamples(unittest.TestCase):
    def test(self):

        # MOOSEDOCS:example-begin
        # Load the packages
        import pyhit
        import moosetree
        # Read the file
        for num in range(3,6):
            root = pyhit.load('grain_growth_2D_bicrystal_addTime.i')

            # Locate and modify "xmax" parameter for the mesh
            Materials = moosetree.find(root, func=lambda n: n.fullpath == '/Materials/CuGrGr')
            # num = 4
            agrs = 
            Materials["wGB"] = num

            # Set the comment on altered parameter
            # mesh.setComment("xmax", "Changed from 3 to 4")

            str_wGB = str(num)
            globalName["my_filename"] = "addTime_wGB" + str_wGB
            input_name = "grain_growth_2D_bicrystal_addTime_wGB" + str_wGB +'.i'
            # Write the modified file
            pyhit.write(input_name, root)

        # MOOSEDOCS:example-end
        # self.assertEqual(mesh["xmax"], num)
        # self.assertEqual(mesh.comment("xmax"), "Changed from 3 to 4")

        # out = mesh.render()
        # self.assertIn("xmax = 4", out)
        # self.assertIn("Changed from 3 to 4", out)

if __name__ == '__main__':
    unittest.main(module=__name__, verbosity=2)

