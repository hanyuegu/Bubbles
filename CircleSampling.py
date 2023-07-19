from dolfin import *
from shapely.geometry import Polygon
import os
import random as rd
from pathlib import Path

import matplotlib.pyplot as plt
import math
import numpy as np

from bubbles.two_d.hole import Circle, Stellar
from bubbles.two_d.topology import Topology

# Safe .GEO files at DATA_DIR
DATA_DIR = Path(__file__).parent.joinpath("data_circle")
DATA_DIR.mkdir(exist_ok=True, parents=True)


def CirlceSampling(Inclusion_area_total, N):

    rd.seed(N)

    """Domain with four sub-rects and many stellar holes."""
    # Domain with many little stellar holes
    # Define the rectangles
    width = 2
    height = 2
    rect_out = [(0, 0), (2, 2)]

    # Set periodic_boundary, True or False for outer rectangle

    topo = Topology(
        rect_out, [], periodic_boundary=False
    )

    # Make some holes
    counter = 0
    i = 0

    stellar_area_total = 0

    Polygon_list = []
    refs = 60
    while True:

        while i < 2000:
            i += 1
            print(i)

            midpoint_s = rd.uniform(0, width), rd.uniform(0, height)
            radius = rd.uniform(0.01, 0.5)
            stellar_hole = Circle(midpoint_s, radius)
            discretized_hole = stellar_hole.discretize_hole(refs)

            if stellar_area_total + Polygon(discretized_hole).area > Inclusion_area_total:
                break

            Boolean_value, discretized_hole = topo.add_hole(stellar_hole, refs=refs, discretized_hole=discretized_hole)

            if Boolean_value:

                pgon = Polygon(discretized_hole)
                Polygon_list.append(discretized_hole)
                stellar_area_total += pgon.area
                counter += 1

                # Draw for debugging
                # x = [p[0] for p in discretized_hole]
                # y = [p[1] for p in discretized_hole]
                # plt.scatter(x, y)
                # plt.savefig("output.png")
                print(i)

        remain_area = Inclusion_area_total - stellar_area_total
        scaling_ratio = np.sqrt(remain_area / Polygon(discretized_hole).area)
        radius = radius * scaling_ratio
        stellar_hole.update_radius(radius)
        discretized_hole = stellar_hole.discretize_hole(refs)
        Boolean_value, discretized_hole = topo.add_hole(stellar_hole, refs=refs, discretized_hole=discretized_hole)

        if Boolean_value:

            pgon = Polygon(discretized_hole)
            Polygon_list.append(discretized_hole)
            stellar_area_total += pgon.area
            counter += 1
            # x = [p[0] for p in discretized_hole]
            # y = [p[1] for p in discretized_hole]
            # plt.scatter(x, y)
            # plt.savefig("output.png")
            print(i)
            break

    print("Needed ", i, "runs for ", counter, "holes")

    # Make geometry, set filled True or False if the inclusion is a hole or not
    topo.mesh(
        file_name=DATA_DIR.joinpath("{}".format(N)),
        write_geo=True,
        lc=1,
        lc_subrects=1,
        filled=True,
    )

    # check Inclusion_area_total
    check_sum = 0
    for P in Polygon_list:
        check_sum += Polygon(P).area
    print("check_sum", check_sum)


CirlceSampling(1.4, 1)  # volumefraction == 0.25
