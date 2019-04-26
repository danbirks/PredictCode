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
    return np.datetime64(start) + np.timedelta64(int(step[:-1]), step[-1].upper())




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



num_knox_iterations = 10

knox_sbin_size = 100 #meters
knox_sbin_num = 12
knox_tbin_size = 7 #days
knox_tbin_num = 8

knox_sbins = makeBins(knox_sbin_size, knox_sbin_num)
knox_tbins = makeBins(knox_tbin_size, knox_tbin_num)




earliest_time = "2017-10-01"
latest_time = "2019-01-02"
time_step = "1M"
time_len = "2M"

start_times = generateDateRange(earliest_time, latest_time, time_step)

for t in start_times:
    print(t)

sys.exit(0)




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





for start_time in start_times:
    
    end_time = generateLaterDate(start_time, time_len)
    print(start_time, end_time)
    
    
    
    
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
    #  - newline
    #  - start date
    #  - end date
    #  - number of events
    #  - "p values"
    #  - one line for each space bin,
    #     holding (space-separated) p values for each time bin
    #  - "statistics"
    #  - one line for each space bin,
    #  -  holding (space-separated) count for each time bin
    #  - "distribution"
    #  - for each space&time bin, counts from n Monte Carlo iterations
    #     *(currently badly formatted, brackets and newlines to strip)
    #  - newline
    outfilebase = "knox_sschi_burg_sbin{}_tbin{}_iter{}.txt".format(
                                                knox_sbin_size, 
                                                knox_tbin_size, 
                                                num_knox_iterations)
    outfilename = os.path.join(datadir, outfilebase)
    
    with open(outfilename,"a") as fout:
        fout.write(str(start_time))
        fout.write("\n")
        fout.write(str(end_time))
        fout.write("\n")
        fout.write(str(num_events))
        fout.write("\n")
        fout.write("p values\n")
        for i in range(knox_sbin_num):
            fout.write(" ".join([str(knox_result.pvalue(i,j)) for j in range(knox_tbin_num)]))
            fout.write("\n")
        fout.write("statistics\n")
        for i in range(knox_sbin_num):
            fout.write(" ".join([str(knox_result.statistic(i,j)) for j in range(knox_tbin_num)]))
            fout.write("\n")
        fout.write("distribution\n")
        for i in range(knox_sbin_num):
            fout.write(" ".join([str(knox_result.distribution(i,j)) for j in range(knox_tbin_num)]))
            fout.write("\n")
        fout.write("\n")
    
    
    
    
    
    
    
    
    
    
    continue
    
    
    
    
    
    
    
    
    mappable = plt.cm.ScalarMappable(cmap=yellow_to_red)
    mappable.set_array(list(all_ratios(knox_result)))
    mappable.autoscale()
    
    
    
    
    
    
    p_thresh_list = [0.10, 0.05, 0.01]
    
    for p_thresh in p_thresh_list:
        
        
        # !!! Add display of pvalues in a grid. Only need 1 grid per dataset,
        #  i.e., 1 per trio of p-thresholds
        
        
        
        fig, ax = plt.subplots(figsize=(16,4))
        
        xmin = min(x for x,y in knox_result.space_bins)
        xmax = max(y for x,y in knox_result.space_bins)
        ax.set(xlim=(xmin, xmax), xlabel="Distance in metres")
        
        _day = np.timedelta64(1,"D")
        ymin = min(t / _day for t,s in knox_result.time_bins)
        ymax = max(s / _day for t,s in knox_result.time_bins)
        ax.set(ylim=(ymin, ymax), ylabel="Time in days")
        
        
        ax.set_title("Knox, {} events from {} to {}, p={}".format(len(points_crime_region_train.timestamps), start_time, end_time, p_thresh))
        
        
        
        for i, space_bin in enumerate(knox_result.space_bins):
            for j, time_bin in enumerate(knox_result.time_bins):
                if knox_result.pvalue(i,j) >= p_thresh:
                    continue
                ratio = knox_ratio(knox_result.statistic(i,j), knox_result.distribution(i,j))
                x, y = space_bin[0], time_bin[0] / _day
                width = space_bin[1] - space_bin[0]
                height = (time_bin[1] - time_bin[0]) / _day
                p = matplotlib.patches.Rectangle((x,y), width, height, fc=mappable.to_rgba(ratio))
                ax.add_patch(p)
                
        cbar = fig.colorbar(mappable, orientation="vertical")
        cbar.set_label("Knox ratio")
        
    
