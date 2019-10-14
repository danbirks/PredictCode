# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:34:27 2019

@author: lawdfo


Purpose:
    Visualise output of spatio-temporal crime risk evaluations. In contrast
    to previous code, the goal here is to not have any hard-coded parameters
    that are specific to a given dataset (e.g., the Chicago data we've been
    using for testing).


"""

import sys
import os
import datetime
import csv
import numpy as np
import time

from collections import Counter
from itertools import product

import geopandas as gpd
import matplotlib.pyplot as plt
from descartes import PolygonPatch
from matplotlib.collections import PatchCollection
import matplotlib


# Elements from PredictCode's custom "open_cp" package
sys.path.insert(0, os.path.abspath(".."))
#import open_cp
#import open_cp.geometry
#import open_cp.sources.chicago as chicago
from open_cp.data import TimedPoints
from open_cp.data import Grid as DataGrid
from open_cp.geometry import intersect_timed_points, \
                             mask_grid_by_intersection
from open_cp.plot import patches_from_grid
import open_cp.prohotspot as phs


# Define color map that ranges from yellow to red
_cdict = {'red':   [(0.0,  1.0, 1.0),
                    (1.0,  1.0, 1.0)],
          'green': [(0.0,  1.0, 1.0),
                    (1.0,  0.0, 0.0)],
          'blue':  [(0.0,  0.2, 0.2),
                    (1.0,  0.2, 0.2)]}
yellow_to_red = matplotlib.colors.LinearSegmentedColormap("yellow_to_red", _cdict)

# Define color map for 5 discrete areas (most significant to least)
discrete_colors = matplotlib.colors.ListedColormap(['red', 'yellow', 'green', 'blue', 'white'])



# Load custom functions that make dealing with datetime and timedelta easier
from crimeRiskTimeTools import generateLaterDate, \
                               generateEarlierDate, \
                               generateDateRange, \
                               getTimedPointsInTimeRange, \
                               getSixDigitDate




# Custom functions



"""
getRegionCells

Purpose:
Generate tuple of all cell coordinates in a region's grid.
Each element corresponds to (row, col) values of a cell --
 that is, the 0-up index of the cell, not geographic coordinates

Input: grid object (open_cp.data.Grid), primarily for its:
        yextent: number of cells in height of region
        xextent: number of cells in width of region
        mask: which cells of the full rectangle are within the area of interest
"""
def getRegionCells(grid):
    # Make sure to do yextent then xextent, because cellcoords
    #  correspond to (row,col) in grid
    all_cells = product(range(grid.yextent), range(grid.xextent))
    return tuple(cc for cc in all_cells 
                 if not grid.mask[cc[0]][cc[1]])



"""
countPointsPerCell

Purpose:
Given a TimedPoints object and a Grid (MaskedGrid?) object,
 return a Counter object that is a mapping from the grid cell
 coordinates to the number of recognised points within the cell.
 Note that "grid cell coordinates" refers to which row of cells
 and which column of cells it's located at, NOT spatial coords.
"""
def countPointsPerCell(points, grid):
    # Get xy coords from TimedPoints
    xcoords, ycoords = points.xcoords, points.ycoords
    # Convert coords to cellcoords
    xgridinds = np.floor((xcoords - grid.xoffset) / grid.xsize).astype(np.int)
    ygridinds = np.floor((ycoords - grid.yoffset) / grid.ysize).astype(np.int)
    # Count the number of crimes per cell
    # NOTE!: We do (y,x) instead of (x,y) because cells are (row,col)
    return Counter(zip(ygridinds, xgridinds))



"""
getHitRateList

Purpose:
Given a sorted list of cells, and a mapping from cells to number of events
 in those cells, return a list of numbers, of length equal to the given
 list of cells, representing the running total of number of events in all
 cells up to that point in the list.
This is a very useful function for measuring the hit rate of a given
 algorithm. An algorithm should output a ranked list of cells, so then we
 can use this function to see how many events the algorithm would have
 detected, given any coverage rate (where the coverage rate corresponds to
 how far along the list we are permitted to search).
"""
def getHitRateList(sorted_cells, cell_hit_map):
    running_total = 0
    hit_rate_list = []
    for cell in sorted_cells:
        running_total += cell_hit_map[cell]
        hit_rate_list.append(running_total)
    return hit_rate_list



"""
Purpose:
    Generate a map that displays the locations of a given set of events
Input:
    points =
    masked_grid = Grid object with a mask
    polygon = 
    title = desired title for output grid
    sizex = width of output grid
            default = 8
    sizey = height of output grid
            default = sizex
Output:
    Display a graph of the region (light blue) with overlaid cell squares
    (black border) and locations of events (red, plus signs), with xy axes as
    geographic coordinates.
"""
def plotPointsOnGrid(points, masked_grid, polygon, title=None, sizex=8, sizey=None):
    
    if sizey == None:
        sizey = sizex
    
    fig, ax = plt.subplots(figsize=(sizex,sizey))
    
    ax.add_patch(PolygonPatch(polygon, fc="none", ec="Black"))
    ax.add_patch(PolygonPatch(polygon, fc="Blue", ec="none", alpha=0.2))
    ax.scatter(points.xcoords,
               points.ycoords,
               marker="+", color="red")
    
    xmin, ymin, xmax, ymax = polygon.bounds
    # Set the axes to have a buffer of 500 around the polygon
    ax.set(xlim=[xmin-500,xmax+500], ylim=[ymin-500,ymax+500])
    
    pc = patches_from_grid(masked_grid)
    ax.add_collection(PatchCollection(pc, facecolor="None", edgecolor="black"))
    if title != None:
        ax.set_title(title)
    
    



"""
plotPointsOnColorGrid
"""

def plotPointsOnColorGrid(polygon, 
                          points, 
                          mesh_info, 
                          value_matrix, 
                          cmap_choice, 
                          title=None, 
                          sizex=8, 
                          sizey=None, 
                          edge_color = "black", 
                          point_color = "black", 
                          point_shape = "+"):
    
    if sizey == None:
        sizey = sizex
    
    fig, ax = plt.subplots(figsize=(sizex,sizey))
    ax.set_aspect(1)
    
    # Color the cells based on the value matrix
    ax.pcolormesh(*mesh_info, value_matrix, cmap=cmap_choice)
    
    # Add outline of region
    ax.add_patch(PolygonPatch(polygon, fc="none", ec="Black"))
    # Plot events
    ax.scatter(points.xcoords,
               points.ycoords,
               marker=point_shape, 
               color=point_color, 
               alpha=0.2)
    
    # Find bounds of the polygon
    xmin, ymin, xmax, ymax = polygon.bounds
    # Set the axes to have a buffer of 500 around the polygon
    ax.set(xlim=[xmin-500,xmax+500], ylim=[ymin-500,ymax+500])
    
    if title != None:
        ax.set_title(title)
    
    return



"""
sortCellsByRiskMatrix
"""

def sortCellsByRiskMatrix(cells, risk_matrix):
    # For each cellcoord, get its risk from the risk matrix
    cellcoord_risk_dict = dict()
    for cc in cells:
        cellcoord_risk_dict[cc] = risk_matrix[cc[0]][cc[1]]
    
    # Sort cellcoords by risk, highest risk first
    cells_risksort = sorted(cells, \
                            key=lambda x:cellcoord_risk_dict[x], \
                            reverse=True)
    return cells_risksort



"""
rankMatrixFromSortedCells

 sorted_cell_list = output from sortCellsByRiskMatrix, e.g.
 cutoff_list = proportions for tiered results, like [0.01,0.02,0.05,0.1]
 score_list = scores to assign to each tier
               length should be 1 more than cutoff list
               by default, evenly spaced from 0 to 1
"""
def rankMatrixFromSortedCells(masked_matrix, 
                              sorted_cell_list, 
                              cutoff_list, 
                              score_list=None):
    if score_list == None:
        #score_list = list(range(len(cutoff_list), -1, -1))
        score_list = np.linspace(0,1,len(cutoff_list)+1)
    if len(score_list) != len(cutoff_list)+1:
        print("Error! Score list is not 1 more than cutoff list: \
              Cutoff:{len(cutoff_list))} vs Score:{len(score_list))}")
        sys.exit(1)
    
    rank_matrix = np.zeros_like(masked_matrix.mask, dtype=float)
    """
    print("type, shape, rank_matrix:")
    print(type(rank_matrix))
    print(rank_matrix.shape)
    print(rank_matrix)
    print("type, shape, masked_matrix:")
    print(type(masked_matrix))
    print(masked_matrix.shape)
    print(masked_matrix)
    """
    
    rank_matrix = masked_matrix.mask_matrix(rank_matrix)
    
    num_cells = len(sorted_cell_list)
    curr_tier = 0
    for i, c in enumerate(sorted_cell_list):
        if masked_matrix.mask[c]:
            print("Error! Cell in sorted list is a masked cell!")
            print(c)
            sys.exit(1)
        while curr_tier < len(cutoff_list) \
              and i/num_cells >= cutoff_list[curr_tier]:
            curr_tier+=1
        rank_matrix[c] = score_list[curr_tier]
    
    return rank_matrix
    
    



"""
runPhsModel
"""
def runPhsModel(training_data, grid, cutoff_time, time_unit, dist_unit, time_bandwidth, dist_bandwidth, weight="linear"):
    
    
    # Obtain model and prediction on grid cells
    phs_predictor = phs.ProspectiveHotSpot(grid=grid)
    phs_predictor.data = training_data
    
    dist_band_in_units = dist_bandwidth/dist_unit
    time_band_in_units = time_bandwidth/time_unit
    
    if weight=="linear":
        phs_predictor.weight = phs.LinearWeightNormalised(space_bandwidth=dist_band_in_units, time_bandwidth=time_band_in_units)
    elif weight=="classic":
        phs_predictor.weight = phs.ClassicWeightNormalised(space_bandwidth=dist_band_in_units, time_bandwidth=time_band_in_units)
    
    phs_predictor.grid = dist_unit
    phs_predictor.time_unit = time_unit
    
    # Only include this method of establishing cutoff_time if we want a
    #  prediction for the day after the latest event in training data. If so,
    #  this will ignore any event-less period of time between training and
    #  test data, which means time decay may be less pronounced.
    #cutoff_time = sorted(training_data.timestamps)[-1] + _day
    
    phs_grid_risk = phs_predictor.predict(cutoff_time, cutoff_time)
    
    
    
    """
    phs_grid_risk_matrix = phs_grid_risk.intensity_matrix
    print("Type of risk matrix:")
    print(type(phs_grid_risk_matrix))
    print("Size of risk matrix:")
    print(phs_grid_risk_matrix.size)
    print("Shape of risk matrix:")
    print(phs_grid_risk_matrix.shape)
    
    md = phs_grid_risk.mesh_data()
    print("Type of mesh_data()")
    print(type(md))
    for i,x in enumerate(md):
        print(f"Type of md-{i}")
        print(type(x))
        print(f"Shape of md-{i}")
        print(x.shape)
        print(x[:5])
        print(x[-5:])
    """
    
    
    
    
    # We might need to mask the risk to the relevant region?
    # Not sure if there's a reason we wouldn't want it masked...
    # But then we would need to pass in masked_grid_region to this function?
    # Or maybe it's sufficient to do it later, if we even need it at all?
    #phs_grid_risk.mask_with(masked_grid_region)
    
    
    #return phs_grid_risk.intensity_matrix
    return grid.mask_matrix(phs_grid_risk.intensity_matrix)
    
    # Sort cellcoords by risk in intensity matrix, highest risk first
    #grid_cells = getRegionCells(grid=grid)
    #return sortCellsByRiskMatrix(grid_cells, phs_grid_risk_matrix)
    



"""
loadGenericData
"""
def loadGenericData(filepath, crime_type_set = {"BURGLARY"}, date_format_csv = "%m/%d/%Y %I:%M:%S %p", epsg = None, proj=None, longlat=True, infeet=False, has_header=True):
    
    # Note: Data is expected in a csv file with the following properties:
    # Row 0 = Header
    # Col 0 = Date/time
    # Col 1 = Longitude (or eastings)
    # Col 2 = Latitude (or northings)
    # Col 3 = Crime type
    # Col 4 = Location type (optional, currently not implemented)
    
    # EPSGs:
    # 3435 or 4326 or 3857 or...? = Chicago
    #   frame.crs = {"init": "epsg:4326"} # standard geocoords
    #   "We'll project to 'web mercator' and use tilemapbase to view the regions with an OpenStreetMap derived basemap"
    #   frame = frame.to_crs({"init":"epsg:3857"})
    # 27700 = UK???
    
    if longlat:
        try:
            import pyproj as _proj
        except ImportError:
            print("Package 'pyproj' not found: projection methods will not be supported.", file=sys.stderr)
            _proj = None
        if not _proj:
            print("_proj did not load!")
            sys.exit(1)
        if not proj:
            if not epsg:
                raise Exception("Need to provide one of 'proj' object or 'epsg' code")
            proj = _proj.Proj({"init": "epsg:"+str(epsg)})
    
    
    _FEET_IN_METERS = 3937 / 1200
    data = []
    with open(filepath) as f:
        csvreader = csv.reader(f)
        if has_header:
            _ = next(csvreader)
        for row in csvreader:
            # Confirm crime type is one we're interested in
            crime_type = row[3].strip()
            if crime_type not in crime_type_set:
                continue
            # Grab time, x, and y values (x & y may be long & lat)
            t = datetime.datetime.strptime(row[0], date_format_csv)
            x = float(row[1])
            y = float(row[2])
            if longlat:
                x, y = proj(x, y)
            else:
                if infeet:
                    x /= _FEET_IN_METERS
                    y /= _FEET_IN_METERS
            # Store data trio
            data.append((t, x, y))
    
    
    data.sort(key = lambda triple : triple[0])
    times = [triple[0] for triple in data]
    xcoords = np.empty(len(data))
    ycoords = np.empty(len(data))
    for i, triple in enumerate(data):
        xcoords[i], ycoords[i] = triple[1], triple[2]
    
    timedpoints = TimedPoints.from_coords(times, xcoords, ycoords)
    
    return timedpoints










# Variables with "chktime" are used to start checking the timing
# Variables with "tkntime" hold the amount of time taken

# Overall timing
chktime_overall = time.time()

#### Declare data parameters
chktime_decparam = time.time()
print("Declaring parameters...")

###
# Assumed parameters

# Location of data file
datadir = os.path.join("..", "..", "Data")
# Header for output file
result_info_header = [
                                "dataset", 
                                "event_types",
                                "cell_width", 
                                "eval_date", 
                                "train_len", 
                                "test_len", 
                                "coverage_rate", 
                                "test_events", 
                                "hit_count", 
                                "hit_pct", 
                                "model", 
                                "rand_seed", 
                                "rhs_bandwidth", 
                                "phs_time_unit", 
                                "phs_time_band", 
                                "phs_dist_unit", 
                                "phs_dist_band", 
                                "phs_weight", 
                                ]
date_today = datetime.date.today()
date_today_str = getSixDigitDate(date_today)


###
# Input parameters for the user to provide
# (Unimplemented parameters are commented out)

# Dataset name (to be included in name of output file)
dataset_name= "chicago"
# Crime types
#!!!  May want to change this to be different for train vs test?
crime_type_set = {"BURGLARY"}
#crime_type_set_sweep = [{"BURGLARY"}] # if doing a sweep over different sets
# Size of grid cells
#cell_width_sweep = [100] # if doing a sweep over different cell sizes
cell_width = 100
# Input csv file name
in_csv_file_name = "chi_all_s_BURGLARY_RES_010101_190101_stdXY.csv"
# Format of date in csv file
date_format_csv = "%m/%d/%Y %I:%M:%S %p"
#!!! Should actually force that format when standardizing csv files!!!
#!!! Also change "infeet" issues at standardizing stage too!!!
# Of all planned experiments, earliest start of a test data set
earliest_test_date = "2003-01-01"
# Time between earliest experiment and latest experiment
test_date_range = "1W"
# Length of training data
train_len = "8W"
# Length of testing data
test_len = "1W"
# Time step offset between different experiments
# (If you want non-overlapping, then set test_date_step = test_len)
#test_date_step = "1D"
test_date_step = test_len
# Coverage rates to test
coverage_bounds = [0.01, 0.02, 0.05, 0.10]
# Geojson file
#geojson_file_name = "Chicago_Areas.geojson"
geojson_file_name = "Chicago_South_Side_2790.geojson"
#geojson_file_name = "Durham_27700.geojson"

# There are more parameters to include regarding the models to run
#  but for now let's just try to display the data itself
# Predictive models to run
# models_to_run = ["naivecount","phs"]



###
# Derived parameters

# Printable string of crime types, concatenated by "_" if necessary
crime_types_printable = "_".join(sorted(crime_type_set))
# Full path for input csv file
in_csv_full_path = os.path.join(datadir, in_csv_file_name)
# Nicely-formatted string of test date
earliest_test_date_str = "".join(earliest_test_date.split("-"))[2:]
# Latest start of a test data set, calculated from earliest and length
latest_test_date = generateLaterDate(earliest_test_date, test_date_range)
# List of all experiment dates
start_test_list = generateDateRange(earliest_test_date, \
                                    latest_test_date, \
                                    test_date_step)
# Number of different experiment dates
total_num_exp_dates = len(start_test_list)
# Output csv file name
out_csv_file_name = "results_{}_{}_{}_{}_{}.csv".format(date_today_str, \
                                                        dataset_name, \
                                                        earliest_test_date_str, \
                                                        test_date_range, \
                                                        test_date_step)
# Full path for output csv file
out_csv_full_path = os.path.join(datadir, out_csv_file_name)
# Full path for geojson file
geojson_full_path = os.path.join(datadir, geojson_file_name)

###
## Not sure if we'll need any of these parameters that had been used previously
#chicago_side = "South"
#chicago_load_type = "snapshot"
#if "all" in chicago_file_name:
#    chicago_load_type = "all"
#chicago_file_path = os.path.join(datadir, chicago_file_name)
## Chicago module requires this line to access some data
#chicago.set_data_directory(datadir)


print("...declared parameters.")
tkntime_decparam = time.time() - chktime_decparam


# Open output csv file for writing, write header row
with open(out_csv_full_path, "w") as csvf:
    writer = csv.writer(csvf, delimiter=",", lineterminator="\n")
    writer.writerow(result_info_header)
    
    
    
    
    
    # If we were to do a "crime type sweep", that would go here.
    # But probably simpler to just re-run with a new crime type set instead.
    # Unless we want to actively combine results within this script?
    
    
    
    
    
    ### OBTAIN FULL DATA
    chktime_obtain_data = time.time()
    print("Obtaining full data set and region...")
    
    #!!! Need to change standardization pre-processing so that this input
    #     is always Eastings/Northings, and in meters
    points_crime = loadGenericData(in_csv_full_path, 
                                   crime_type_set=crime_type_set, 
                                   longlat=False, 
                                   infeet=True)
    
    # Obtain polygon from geojson file (which should have been pre-processed)
    region_polygon = gpd.read_file(geojson_full_path).unary_union
    
    # Get subset of input crime that occurred within region
    points_crime_region = intersect_timed_points(points_crime, region_polygon)
    
    # Get grid version of region
    masked_grid_region = mask_grid_by_intersection(region_polygon, 
                                                   DataGrid(xsize=cell_width, 
                                                            ysize=cell_width, 
                                                            xoffset=0, 
                                                            yoffset=0)
                                                   )
    
    # Get "mesh info" of that grid, which is useful for displaying the map
    masked_grid_mesh = masked_grid_region.mesh_info()
    
    
    # Get tuple of all cells in gridded region
    cellcoordlist_region = getRegionCells(masked_grid_region)
    
    print("...obtained full data set and region.")
    tkntime_obtain_data = time.time() - chktime_obtain_data
    print(f'Time taken to obtain data: {tkntime_obtain_data}')
    
    
    
    # How frequently to display a print statement showing the experiment number
    print_exp_freq = 5
    
    # Log of how long each experiment takes to run
    exp_times = []
    
    
    # Each start_test time in the generated list defines the time period for
    #  an experiment (or multiple experiments)
    for exp_date_index, start_test in enumerate(start_test_list):
        
        chktime_exp = time.time()
        
        if exp_date_index % print_exp_freq == 0:
            print(f"Running experiment {exp_date_index+1}/{total_num_exp_dates}...")
        
        # Compute time ranges of training and testing data
        end_train = start_test
        start_train = generateEarlierDate(end_train, train_len)
        end_test = generateLaterDate(start_test, test_len)
        
        
        
        # Obtain training data
        points_crime_region_train = getTimedPointsInTimeRange(points_crime_region, 
                                                              start_train, 
                                                              end_train)
        # Count how many crimes there were in this training data set
        num_crimes_train = len(points_crime_region_train.timestamps)
        
        
        # Obtain testing data
        points_crime_region_test = getTimedPointsInTimeRange(points_crime_region, 
                                                              start_test, 
                                                              end_test)
        # Count how many crimes there were in this test data set
        num_crimes_test = len(points_crime_region_test.timestamps)
        
        
        # Count the number of crimes per cell in test data
        #  This is used for evaluation and also for the "ideal" model
        cells_testcrime_ctr = countPointsPerCell(points_crime_region_test, 
                                                 masked_grid_region)
        
        
        print(f"num_crimes_train: {num_crimes_train}")
        print(f"num_crimes_test: {num_crimes_test}")
        
        
        
        
        
        # Plot training data on plain grid
        plotPointsOnGrid(points_crime_region_train, 
                         masked_grid_region, 
                         region_polygon, 
                         title=f"Train data {exp_date_index+1}", 
                         sizex=8, 
                         sizey=8)
        
        # Plot testing data on plain grid
        plotPointsOnGrid(points_crime_region_test, 
                         masked_grid_region, 
                         region_polygon, 
                         title=f"Test data {exp_date_index+1}", 
                         sizex=8, 
                         sizey=8)
        
        
        
        
        
        
        
        # Count the number of crimes per cell in training data
        #  This is used for the naive model
        cells_traincrime_ctr = countPointsPerCell(points_crime_region_train, 
                                                  masked_grid_region)
        
        
        
        
        naive_data_matrix = np.zeros([masked_grid_region.yextent, masked_grid_region.xextent])
        naive_data_matrix = masked_grid_region.mask_matrix(naive_data_matrix)
        for c in cells_traincrime_ctr:
            naive_data_matrix[c] = cells_traincrime_ctr[c]
        
        
        
        
        plotPointsOnColorGrid(polygon = region_polygon, 
                              points = points_crime_region_train, 
                              mesh_info = masked_grid_mesh, 
                              value_matrix = naive_data_matrix, 
                              cmap_choice = yellow_to_red, 
                              title="Naive map", 
                              sizex=10, 
                              sizey=10)
        
        
        
        
        
        # Still need to replace hard-coded parameters in this call here
        phs_data_matrix = runPhsModel(training_data=points_crime_region_train, 
                                   grid=masked_grid_region, 
                                   cutoff_time=start_test, 
                                   time_unit=np.timedelta64(1, "W"), 
                                   dist_unit=cell_width, 
                                   time_bandwidth=np.timedelta64(4, "W"), 
                                   dist_bandwidth= 4 * cell_width, 
                                   # ^^^ That should be set as a param
                                   weight="classic")
                                   #weight="linear")
        #phs_data_matrix = masked_grid_region.mask_matrix(phs_data_matrix)
        
        
        plotPointsOnColorGrid(polygon = region_polygon, 
                              points = points_crime_region_train, 
                              mesh_info = masked_grid_mesh, 
                              value_matrix = phs_data_matrix, 
                              cmap_choice = yellow_to_red, 
                              title="PHS map", 
                              sizex=10, 
                              sizey=10)
        
        
        
        
        
        
        
        
        sorted_cells_naive = sortCellsByRiskMatrix(cellcoordlist_region, 
                                                   naive_data_matrix)
        
        naive_rank_matrix = rankMatrixFromSortedCells(masked_grid_region, 
                                                      sorted_cells_naive, 
                                                      coverage_bounds)
        
        plotPointsOnColorGrid(polygon = region_polygon, 
                              points = points_crime_region_train, 
                              mesh_info = masked_grid_mesh, 
                              value_matrix = naive_rank_matrix, 
                              cmap_choice = discrete_colors, 
                              title="Naive map ranked by coverage", 
                              sizex=10, 
                              sizey=10)
        
        
        sorted_cells_phs = sortCellsByRiskMatrix(cellcoordlist_region, 
                                                   phs_data_matrix)
        
        phs_rank_matrix = rankMatrixFromSortedCells(masked_grid_region, 
                                                      sorted_cells_phs, 
                                                      coverage_bounds)
        
        plotPointsOnColorGrid(polygon = region_polygon, 
                              points = points_crime_region_train, 
                              mesh_info = masked_grid_mesh, 
                              value_matrix = phs_rank_matrix, 
                              cmap_choice = discrete_colors, 
                              title="PHS map, ranked by coverage", 
                              sizex=10, 
                              sizey=10)
        
        
        
        
        
        
        sorted_cells_ideal = sorted(cellcoordlist_region, 
                                    key=lambda x:cells_testcrime_ctr[x], 
                                    reverse=True)
        
        
        
        ideal_rank_matrix = rankMatrixFromSortedCells(masked_grid_region, 
                                                      sorted_cells_ideal, 
                                                      coverage_bounds)
        
        plotPointsOnColorGrid(polygon = region_polygon, 
                              #points = points_crime_region_train, 
                              points = points_crime_region_test, 
                              mesh_info = masked_grid_mesh, 
                              value_matrix = ideal_rank_matrix, 
                              cmap_choice = discrete_colors, 
                              title="Ideal map, ranked by coverage", 
                              sizex=10, 
                              sizey=10)
        
        
        
        sys.exit(0)
        
        
        
        
        
        
        #hit_rate_list = getHitRateList(sorted_cells, cells_testcrime_ctr)
        
        
        
        
        
        
        
        
        tkn_time_exp = time.time()-chktime_exp
        print(f"time spent on exp: {tkn_time_exp}")
        exp_times.append(tkn_time_exp)
        
    
    
    
    
    
    
    





