# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 16:26:19 2019

@author: Dustin
"""

# JUST KNOX STUFF




# Some fairly standard modules
import os, csv, lzma
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import descartes
from itertools import product
from collections import Counter
import random
import statistics
import time

# The geopandas module does not come standard with anaconda,
# so you'll need to run the anaconda prompt as an administrator
# and install it via "conda install -c conda-forge geopandas".
# That installation will include pyproj and shapely automatically.
# These are useful modules for plotting geospatial data.
import geopandas as gpd
import pyproj
import shapely.geometry

# These modules are useful for tracking where modules are
# imported from, e.g., to check we're using our local edited
# versions of open_cp scripts.
import sys
import inspect
import importlib

# In order to use our local edited versions of open_cp
# scripts, we insert the parent directory of the current
# file ("..") at the start of our sys.path here.
sys.path.insert(0, os.path.abspath(".."))

# Elements from PredictCode's custom "open_cp" package
import open_cp
import open_cp.geometry
import open_cp.plot
import open_cp.sources.chicago as chicago
import open_cp.retrohotspot as retro
import open_cp.prohotspot as phs
import open_cp.knox





"""
Given a TimedPoints object that contains spatial coordinates paired
    with timestamps, return a similar object that only contains the
    data whose timestamps are within the given range
"""
def getTimedPointsInTimeRange(points, start_time, end_time):
    return points[(points.timestamps >= start_time)
                  & (points.timestamps <= end_time)]


"""
getRegionCells
Given a grid, return a list of all cells in the grid that are part of the
    defined region, i.e., that are not covered by the "mask" of the grid
"""
def getRegionCells(grid):
    # Make sure to do yextent then xextent, because cellcoords
    #  correspond to (row,col) in grid
    all_cells = product(range(grid.yextent), range(grid.xextent))
    return tuple(cc for cc in all_cells 
                 if not grid.mask[cc[0]][cc[1]])




def makeBins(size, num):
    return list( (i*size, (i+1)*size) for i in range(num))





def getKnoxResult(data, num_iter, sbins, tbins, tbin_unit="days"):
    
    knox = open_cp.knox.Knox()
    knox.set_time_bins(tbins, unit=tbin_unit)
    knox.space_bins = sbins
    knox.data = data
    result = knox.calculate(iterations=num_iter)
    return result
    




"""
knox_ratio
Computes the "knox ratio", defined as being:
    The "statistic" (# events in a given pair of distance and time bins)
    divided by
    the median #events in that bin pair among the randomized trials
Note: If the #trials is even, the "median" is actually the entry preceding
    the direct center of the list, rather than the average of the two entries
    on either side.
"""
def knox_ratio(statistic, distribution):
    """As in the paper, compute the ratio of the statistic to the median
    of the values in the distribution"""
    d = np.array(distribution)
    d.sort()
    return statistic / d[len(d)//2]





"""
all_ratios
Given a Knox result, this constructs a generator that returns Knox ratios for
    all (space,time) bin pairs
"""
def all_ratios(result):
    for i, space_bin in enumerate(result.space_bins):
        for j, time_bin in enumerate(result.time_bins):
            yield knox_ratio(result.statistic(i,j), result.distribution(i,j))
            






"""
#
# start: "2018-09-01"
# end: "2019-03-01"
# step: "3D", "2W", "6M", "1Y", etc
# output is always array containing type datetime64[D]
"""
def generateDateRange(start, end, step):
    step_num = int(step[:-1])
    step_unit = step[-1].upper()
    if step_unit not in "YMWD":
        print("generateDateRange only recognizes Y/M/W/D as units")
        sys.exit(1)
    if step_unit == "W":
        step_num *= 7
        step_unit = "D"
    start_datetime = np.datetime64(start)
    end_datetime = np.datetime64(end)
    date_range = np.arange(start_datetime, end_datetime, step=np.timedelta64(step_num, step_unit), dtype="datetime64[{}]".format(step_unit))
    return np.array(date_range, dtype="datetime64[D]")




def generateLaterDate(start, step):
    time_base = step[-1].upper()
    time_quantity = int(step[:-1])
    converted_start = np.datetime64(start, time_base)
    later_date = converted_start + np.timedelta64(time_quantity, time_base)
    return np.datetime64(later_date, "D")




"""
bin size
p value
diff btw time windows


"""



init_clock_time = time.time()


#Constants

# Heat-like color mapping
_cdict = {'red':   [(0.0,  1.0, 1.0),
                   (1.0,  1.0, 1.0)],
         'green': [(0.0,  1.0, 1.0),
                   (1.0,  0.0, 0.0)],
         'blue':  [(0.0,  0.2, 0.2),
                   (1.0,  0.2, 0.2)]}

yellow_to_red = matplotlib.colors.LinearSegmentedColormap("yellow_to_red", _cdict)


_day = np.timedelta64(1,"D")







clock_time = time.time()

print("Loading parameters...")

# Data location parameters

datadir = os.path.join("..", "..", "Data")
chicago_file_name = "chicago_all_old.csv"
chicago_side = "South"
#crime_type_set = {"THEFT"}
crime_type_set = {"BURGLARY"}
chicago_file_path = os.path.join(datadir, chicago_file_name)
chicago.set_data_directory(datadir)



# Test parameters

cell_width = 250
cell_height = cell_width



num_knox_iterations = 500

knox_sbin_size = 100 #meters
knox_sbin_num = 15
knox_tbin_size = 7 #days
knox_tbin_num = 10

knox_sbins = makeBins(knox_sbin_size, knox_sbin_num)
knox_tbins = makeBins(knox_tbin_size, knox_tbin_num)




earliest_time = "2017-05-01"
latest_time = "2017-08-02"
time_step = "6M"
time_len = "12M"
date_ref = "".join(earliest_time[2:].split("-"))


start_times = generateDateRange(earliest_time, latest_time, time_step)

num_exp = len(start_times)

print(start_times[0])
print(start_times[-1])
print("Num exp: {}".format(num_exp))




outfilebase = \
    "knox_ssX_burg_sbin{}-{}_tbin{}-{}_iter{}_{}-{}_{}-{}.txt".format(
                                            knox_sbin_num, 
                                            knox_sbin_size, 
                                            knox_tbin_num, 
                                            knox_tbin_size, 
                                            num_knox_iterations, 
                                            date_ref, 
                                            time_step, 
                                            num_exp, 
                                            time_len)
outfilename = os.path.join(datadir, outfilebase)

print("outfile: {}".format(outfilename))




print("...loaded parameters.\nTime taken: {}".format(time.time() - clock_time))

print("Loading data...")

clock_time = time.time()

#points_crime = chicago.load(chicago_file_path, crime_type_set)
points_crime = chicago.load(chicago_file_path, crime_type_set, type="all")


print("...loaded data.\nTime taken: {}".format(time.time() - clock_time))


### OBTAIN GRIDDED REGION

clock_time = time.time()
print("Loading region and data subset...")

# Obtain polygon shapely object for region of interest
region_polygon = chicago.get_side(chicago_side)

# Obtain data set
points_crime_region = open_cp.geometry.intersect_timed_points(points_crime, region_polygon)

# Obtain grid with cells only overlaid on relevant region
masked_grid_region = open_cp.geometry.mask_grid_by_intersection(
        region_polygon, 
        open_cp.data.Grid(
                xsize=cell_width, 
                ysize=cell_height, 
                xoffset=0, 
                yoffset=0))

# Get a list/tuple of all cellcoords in the region
cellcoordlist_region = getRegionCells(masked_grid_region)

# Obtain number of cells in the grid that contain relevant geometry
# (i.e., not the full rectangular grid, only relevant cells)
num_cells_region = len(cellcoordlist_region)

print("...loaded region and data subset.\nTime taken: {}".format(time.time() - clock_time))



all_exp_clock_time = time.time()
exp_clock_time_list = []

print("Starting experiments...")

for exp_index, start_time in enumerate(start_times):
    
    exp_clock_time = time.time()
    
    end_time = generateLaterDate(start_time, time_len)
    
    print("Exp {}/{}, start {} end {}".format(exp_index+1, 
                                              num_exp, start_time, end_time))
    
    
    
    
    ### SELECT TRAINING DATA
    
    
    # Get subset of data for training
    points_crime_region_train = getTimedPointsInTimeRange(points_crime_region, 
                                                      start_time, 
                                                      end_time)
    
    
    
    num_events = len(points_crime_region_train.timestamps)
    
    print("Number of events: {}".format(num_events))
    
    
    knox_result = getKnoxResult(points_crime_region_train, 
                                num_knox_iterations, 
                                knox_sbins, 
                                knox_tbins)
    
    
    
    
    
    
    
    
    # This section of code will output the Knox statistic info to a text file,
    #  containing the following for each experiment:
    #  - start date
    #  - end date
    #  - number of events
    #  - "Knox Statistics"
    #  - one line for each time bin,
    #  -  holding (space-separated) count for each space bin
    #  - "Monte Carlo Medians"
    #  - one line for each time bin,
    #  -  holding (space-separated) median Monte Carlo count for each space bin
    #  - "P Values"
    #  - one line for each time bin,
    #     holding (space-separated) p values for each space bin
    #  - newline
    
    
    with open(outfilename,"a") as fout:
        fout.write(str(start_time))
        fout.write("\n")
        fout.write(str(end_time))
        fout.write("\n")
        fout.write(str(num_events))
        fout.write("\n")
        fout.write("Knox Statistics\n")
        for i in range(knox_tbin_num):
            fout.write(" ".join([str(knox_result.statistic(j,i)) for j in range(knox_sbin_num)]))
            fout.write("\n")
        fout.write("Monte Carlo Medians\n")
        for i in range(knox_tbin_num):
            fout.write(" ".join([str(statistics.median(knox_result.distribution(j,i))) for j in range(knox_sbin_num)]))
            fout.write("\n")
        fout.write("P Values\n")
        for i in range(knox_tbin_num):
            fout.write(" ".join([str(knox_result.pvalue(j,i)) for j in range(knox_sbin_num)]))
            fout.write("\n")
        fout.write("\n")
    
    
    exp_clock_time_list.append(time.time()-exp_clock_time)
    

exp_clock_time_list.sort()

print("...Finished experiments.")
print("Min time: {}".format(exp_clock_time_list[0]))
print("Max time: {}".format(exp_clock_time_list[-1]))
print("Avg time: {}".format(sum(exp_clock_time_list)/len(exp_clock_time_list)))
print("All exp time: {}".format(time.time()-all_exp_clock_time))
print("Total script time: {}".format(time.time()-init_clock_time))