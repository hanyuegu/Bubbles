import os
import random as rd
from dolfin import *
import matplotlib.pyplot as plt
from bubbles.geo_utils.utils import write_geo
from bubbles.two_d.hole import Circle, Stellar
from bubbles.two_d.topology import Topology
from shapely.geometry import Polygon
import random
import math
import numpy as np

# random.seed(20)

def CircleSampling(N):
    # Domain with many little stellar holes
    # Define the rectangles
    width = 2
    height = 2
    rect_out = [(0, 0), (2, 2)]

    # Set periodic_boundary, True or False for outer rectangle

    periodic_boundary = True
    topo = Topology(rect_out, [], periodic_boundary)

    Inclusion_area_total = 1.4  # volumefraction == 0.35

    random_list = [random.uniform(0, Inclusion_area_total) for _ in range(N-1)]
    random_list.append(0)
    random_list.sort()
    area_list = np.diff(random_list) # area of first N-1 circle
    
    counter = 0
    i = 0
    circal_area_total = 0
    area_index = 0
    holes = []
    
    if periodic_boundary:
        while counter < N-1 and i < 1000:
            i += 1
            print(i)

            midpoint_s = rd.uniform(0, width), rd.uniform(0, height)
            radius = sqrt(area_list[area_index]/math.pi)
            stellar_hole = Circle(midpoint_s, radius)
            Boolean_value, discretized_hole = topo.add_hole(stellar_hole, refs=60)
        
            if Boolean_value:
                circal_area_total += math.pi * radius ** 2
                hole = {}
                hole["midpoint"] = midpoint_s
                hole["radius"] = radius
                holes.append(hole)
                counter += 1
                area_index += 1

        # generate last circle
        while counter < N and i < 1000:
            i += 1
            print(i)
            radius_N = sqrt((Inclusion_area_total - circal_area_total)/math.pi)
            midpoint_N = rd.uniform(0, width-radius_N), rd.uniform(0, height-radius_N)

            Flag = False
            for h in holes:
                if math.dist(h["midpoint"], midpoint_N) < (h["radius"] + radius_N):
                    Flag = True
            if Flag:
                continue
  
            stellar_hole = Circle(midpoint_N, radius_N)
            Boolean_value, discretized_hole = topo.add_hole(stellar_hole, refs=60)

            if Boolean_value:
                hole = {}
                hole["midpoint"] = midpoint_N
                hole["radius"] = radius_N
                holes.append(hole)
                counter += 1

        print("Needed ", i, "runs for ", N, "holes")

        # check Inclusion_area_total
        # check_sum = 0
        # for h in holes:
        #     check_sum += math.pi * h["radius"] ** 2
        # print("check_sum", check_sum)

    else:
        # Make some holes
        Polygon_list = []
        while counter < N-1 and i < 1000:
            i += 1
            print(i)

            midpoint_s = rd.uniform(0, width), rd.uniform(0, height)
            radius = sqrt(area_list[area_index]/math.pi)
            stellar_hole = Circle(midpoint_s, radius)
            Boolean_value, discretized_hole = topo.add_hole(stellar_hole, refs=60)
        
            if Boolean_value:
                print(discretized_hole)
                pgon = Polygon(discretized_hole)
                Polygon_list.append(discretized_hole)
                circal_area_total += pgon.area
                hole = {}
                hole["midpoint"] = midpoint_s
                hole["radius"] = radius
                holes.append(hole)
                counter += 1
                area_index += 1

        # generate last circle
        while counter < N and i < 1000:
            i += 1
            print(i)
            radius_N = sqrt((Inclusion_area_total - circal_area_total)/math.pi)
            midpoint_N = rd.uniform(0, width-radius_N), rd.uniform(0, height-radius_N)
        
            Flag = False
            for h in holes:
                if math.dist(h["midpoint"], midpoint_N) < (h["radius"] + radius_N):
                    Flag = True
            if Flag:
                continue

            stellar_hole = Circle(midpoint_N, radius_N)
            Boolean_value, discretized_hole = topo.add_hole(stellar_hole, refs=60)
            if Boolean_value:
                Polygon_list.append(discretized_hole)
                counter += 1

        print("Needed ", i, "runs for ", N, "holes")

        # check Inclusion_area_total
        # check_sum = 0
        # for P in Polygon_list:
        #     check_sum += Polygon(P).area
        # print("check_sum", check_sum)

    # Make geometry, set filled True or False if the inclusion is a hole or not
    geo = topo.get_geometry(lc=1, lc_subrects=1, filled=True)
    write_geo(geo, "Circle_{}".format(N))

# CircleSampling(2)