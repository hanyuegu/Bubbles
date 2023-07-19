from dolfin import *
import os
from VolumeFractionSampling import *

for i in range(1000):
    StellarSampling(1.4, i+1)

# for i in range(1,1001):
#     os.system("gmsh -2 -format msh2 data_stellar/{}.geo_unrolled -o data_stellar/{}.msh".format(i, i))
#     os.system("dolfin-convert data_stellar/{}.msh data_stellar/{}.xml".format(i, i))