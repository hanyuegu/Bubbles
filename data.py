from dolfin import *     
import os
from debug import *
from tests.test_geometry_three_d import *
from tests.test_geometry_two_d import *
from VolumeFractionSampler import *


for i in range(50):
    CircleSampling(i+1)

for i in range(50):
    os.system("gmsh -2 -format msh2 Circle_{}.geo -o Circle_{}.msh".format(i + 1, i + 1))
    os.system("dolfin-convert Circle_{}.msh Circle_{}.xml".format(i + 1, i + 1))