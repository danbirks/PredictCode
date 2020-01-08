# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:34:27 2019

@author: lawdfo

riskModelsGeneric.py

Purpose:
    Run spatio-temporal crime risk models and visualise the output. 
    In contrast to previous code, the goal here is to not have any hard-coded 
    parameters that are specific to a given dataset (e.g., the Chicago data 
    we've been using for testing).


"""

import sys
import os
import datetime
import csv
import numpy as np
import time

from collections import Counter, defaultdict
from itertools import product
from copy import deepcopy

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



# Load custom functions that make dealing with datetime and timedelta easier
from crimeRiskTimeTools import generateLaterDate, \
                               generateEarlierDate, \
                               generateDateRange, \
                               getTimedPointsInTimeRange, \
                               getSixDigitDate, \
                               shorthandToTimeDelta
import geodataTools as gdt


###
# Assumed parameters

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
# Columns risk_MODEL and rank_MODEL for each.model will be appended too
risk_info_header = [
                        "dataset", 
                        "event_types",
                        "cell_width", 
                        "eval_date", 
                        "train_len", 
                        "test_len", 
                        "rownum", 
                        "colnum", 
                        "easting", 
                        "northing" ]

date_today = datetime.date.today()
date_today_str = getSixDigitDate(date_today)

recognised_models = ["random", "naive", "phs", "ideal"]








# Custom functions



"""
splitCommaArgs
Given a string as an argument meant to be comma-separated, but possibly
 containing spaces that need to be stripped away, return a list of those
 strings. (This was common enough that it warranted its own function to
 save a bit of space.)
"""
def splitCommaArgs(argstring):
    return [x.strip() for x in argstring.split(",")]


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
getHitCountList

Purpose:
Given a sorted list of cells, and a mapping from cells to number of events
 in those cells, return a list of numbers, of length equal to the given
 list of cells +1, representing the running total of number of events in all
 cells up to that point in the list.
This is a very useful function for measuring the hit rate of a given
 algorithm. An algorithm should output a ranked list of cells, so then we
 can use this function to see how many events the algorithm would have
 detected, given any coverage rate (where the coverage rate corresponds to
 how far along the list we are permitted to search).
Note that the list starts with 0, in the event that the coverage is so low
 that no cells are checked. If the coverage allows you to check n cells, then
 the hit rate will be the value at index n (0-up).
"""
def getHitCountList(sorted_cells, cell_hit_map):
    running_total = 0
    hit_count_list = [0]
    for cell in sorted_cells:
        running_total += cell_hit_map[tuple(cell)]
        hit_count_list.append(running_total)
    return hit_count_list



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
def plotPointsOnGrid(points, 
                     masked_grid, 
                     polygon, 
                     title=None, 
                     sizex=10, 
                     sizey=None, 
                     out_img_file_path=None):
    
    if sizey == None:
        sizey = sizex
    
    fig, ax = plt.subplots(figsize=(sizex,sizey))
    ax.set_aspect(1)
    
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
    
    
    
    
    if out_img_file_path != None:
        fig.savefig(out_img_file_path)
    print(f"Saved image file: {out_img_file_path}")
    
    return



"""
plotPointsOnColorGrid
"""

def plotPointsOnColorGrid(polygon, 
                          points, 
                          mesh_info, 
                          value_matrix, 
                          cmap_choice, 
                          title=None, 
                          sizex=10, 
                          sizey=None, 
                          edge_color = "black", 
                          point_color = "black", 
                          point_shape = "+", 
                          out_img_file_path = None):
    
    
    
    
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
    
    if out_img_file_path != None:
        fig.savefig(out_img_file_path)
    
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
hitRatesFromHitList
"""
def hitRatesFromHitList(hit_count_list, 
                        coverage_rate, 
                        num_cells, 
                        num_crimes_test):
    num_hits = hit_count_list[int(coverage_rate * num_cells)]
    pct_hits = 0
    if num_crimes_test>0:
        pct_hits = num_hits / num_crimes_test
    return num_hits, pct_hits



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
runNaiveModel

Runs the naive model in which we simply count the number of events in the
 training data per cell.
Returns an intensity matrix, where each cell of the (masked) grid is assigned
 a risk score (which is equal to the number of events in the training data).
This is also used for the ideal model, using testing data as training data.

training_data   : 
grid            : 
"""
def runNaiveModel(training_data, grid):
    # Count the number of crimes per cell in training data
    cells_traincrime_ctr = countPointsPerCell(training_data, grid)
    
    naive_data_matrix = np.zeros([grid.yextent, 
                                  grid.xextent])
    naive_data_matrix = grid.mask_matrix(naive_data_matrix)
    for c in cells_traincrime_ctr:
        naive_data_matrix[c] = cells_traincrime_ctr[c]
    
    return naive_data_matrix

def runRandomModel(grid, rand_seed):
    
    np.random.seed(seed=rand_seed)
    random_data_matrix = np.random.rand(grid.yextent, 
                                        grid.xextent)
    random_data_matrix = grid.mask_matrix(random_data_matrix)
    
    return random_data_matrix


"""
runPhsModel

Runs the Prospective Hotspotting model.
Relies on importing open_cp.prohotspot as phs.
Returns an intensity matrix, where each cell of the (masked) grid is assigned
 a risk score.

training_data   : 
grid            : 
cutoff_time     : 
time_unit       : basic unit for time (examples: 3D, 2W, 6M, 1Y)
dist_unit       : basic unit for distance, in meters (ex: 100)
time_bandwidth  : multiple of time_unit used for bandwidth, as an integer
dist_bandwidth  : multiple of dist_unit used for bandwidth, as an integer
weight          : "linear" or "classic"
                  "linear" = use phs.LinearWeightNormalised
                  "classic" = use phs.ClassicWeightNormalised
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
    
    

    
    
    
    # Mask the risk matrix to the relevant region, return it
    return grid.mask_matrix(phs_grid_risk.intensity_matrix)
    


"""
saveModelResultMaps
After running a model, save visualisations of how it mapped risk, both as a
general heat map and in different bins of coverage.
"""
def saveModelResultMaps(model_name, 
                        data_matrix, 
                        rank_matrix, 
                        exp_ident, 
                        file_core, 
                        filedir, 
                        polygon, 
                        points_to_map, 
                        mesh_info, 
                        ):
    
    
    # Define file names
    
    img_file_heat_name = "_".join(["heatmap", 
                                   model_name, 
                                   file_core, 
                                   exp_ident])
    img_file_heat_name += ".png"
    img_file_heat_fullpath = os.path.join(filedir, img_file_heat_name)
    
    
    img_file_cov_name = "_".join(["covmap", 
                                   model_name, 
                                   file_core, 
                                   exp_ident])
    img_file_cov_name += ".png"
    img_file_cov_fullpath = os.path.join(filedir, img_file_cov_name)
    
    model_cap = model_name.capitalize()
    heat_title = f"{model_cap} heat map {exp_ident}"
    cov_title = f"{model_cap} coverage map {exp_ident}"
    
    
    # Save risk heat map
    plotPointsOnColorGrid(polygon = polygon, 
                          points = points_to_map, 
                          mesh_info = mesh_info, 
                          value_matrix = data_matrix, 
                          cmap_choice = yellow_to_red, 
                          title=heat_title, 
                          sizex=10, 
                          sizey=10, 
                          out_img_file_path = img_file_heat_fullpath)
    
    # Save coverage map
    plotPointsOnColorGrid(polygon = polygon, 
                          points = points_to_map, 
                          mesh_info = mesh_info, 
                          value_matrix = rank_matrix, 
                          cmap_choice = discrete_colors, 
                          title=cov_title, 
                          sizex=10, 
                          sizey=10, 
                          out_img_file_path = img_file_cov_fullpath)
    
    return





"""
graphCoverageVsHitRate


"""
def graphCoverageVsHitRate(hit_rates_dict, 
                           model_runs_list, 
                           model_names, 
                           x_limits = None, 
                           title = None, 
                           out_img_file_path = None):
    
    model_hit_rate_pairs = []
    for mn in model_names:
        model_hit_rate_pairs += list(zip(model_runs_list[mn], hit_rates_dict[mn]))
    
    
    # Declare figure
    print("Declaring figure for graphCoverageVsHitRate...")
    fig, ax = plt.subplots(figsize=(12,6))
    
    
    names_for_legend = []
    
    x_axis_size = len(hit_rates_dict[model_names[0]][0])
    x_axis_values = np.linspace(0,1,x_axis_size)
    
    
    for mn in model_names:
        for hr in hit_rates_dict[mn]:
            ax.plot(x_axis_values, hr)
        for mr in model_runs_list[mn]:
            names_for_legend.append(mr)
    
    ax.legend(names_for_legend)
    
    ax.set(xlabel="Coverage", ylabel="Hit rate")
    
    # Set the axes to have a buffer of 500 around the polygon
    if x_limits != None:
        ax.set(xlim=x_limits)
    
    if title != None:
        ax.set_title(title)
    
    if out_img_file_path != None:
        fig.savefig(out_img_file_path)
    
    return
    
    





"""
loadGenericData
"""
def loadGenericData(filepath, 
                    crime_type_set = {"BURGLARY"}, 
                    date_format_csv = "%m/%d/%Y %I:%M:%S %p", 
                    epsg = 4326, # standard lat/long aka WGS84
                    proj=None, 
                    longlat=False, 
                    infeet=True, 
                    has_header=True):
    
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
    
    
    #!!! Need to enforce "%m/%d/%Y %I:%M:%S %p" format for csv files!!!
    #!!! Also change "infeet" issues at standardizing stage too!!!
    
    
    
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




















"""
Former section of setting parameters here;
 it will be removed after confirming it's no longer needed.




###
# Input parameters for the user to provide
# (Unimplemented parameters are commented out)


# Location of data file
datadir = os.path.join("..", "..", "Data")
# Dataset name (to be included in name of output file)
dataset_name= "chicago"
# Crime types
#!!!  May want to change this to be different for train vs test?
#!!!  May want to use multiple sets of types, then combine results?
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
# Of all planned experiments, earliest start of a TEST (not train) data set
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

# How frequently to display a print statement showing the experiment number
print_exp_freq = 5


# Predictive models to run
models_to_run = ["random","naive","phs","ideal"]

num_random = 3

model_param_dict = dict()
model_param_dict["ideal"] = [()]
model_param_dict["naive"] = [()]

# Param list for Random model.
#  For example, if 4, param list is [(0,), (1,), (2,), (3,)]
model_param_dict["random"] = list(product(range(num_random)))



# Param list for PHS model.
phs_time_units = ["1W"]
phs_time_bands = ["4W"]
phs_dist_units = [100]
phs_dist_bands = [400]
#phs_weight = ["linear"]
phs_weight = ["classic"]
model_param_dict["phs"] = list(product(
                            phs_time_units, 
                            phs_time_bands, 
                            phs_dist_units, 
                            phs_dist_bands, 
                            phs_weight))
"""










"""
runModelExperiments

Run a set of models, with various sets of parameters, on multiple time
 window slices of a data set.

Arguments:
    datadir_in : 
        String representing path to data directory.
            The directory should hold the input data for training and
            testing, and the geojson file.
            Ex: "../../Data"
    dataset_name_in : 
        String for part of output file names, to indicate dataset used.
            Ex: "chicago"
    crime_type_set_in : 
        Comma-separated string of crime types as labeled in the input data.
            Exs: "BURGLARY" or "BURGLARY,THEFT"
    cell_width_in : 
        Integer or string for length of cell sides. Will be casted to 
            an integer.
            Ex: 100
    in_csv_file_name_in : 
        Name of input csv file for training and testing.
            Ex: "chi_all_s_BURGLARY_RES_010101_190101_stdXY.csv"
    earliest_test_date_in : 
        String representing date for the cutoff between the training data
            and testing data for the earliest experiment to run.
            Format is "YYYY-MM-DD"
            Ex: "2003-01-31"
    test_date_range_in : 
        String representing the span of time from the date of the earliest
            experiment to be run (earliest_test_date_in) to the date of the
            latest experiment to be run. Format is a timespan shorthand used
            throughout these scripts, as a number followed by D/W/M/Y for
            days/weeks/months/years.
            Exs: "1W" or "6M" or "5Y"
    train_len_in : 
        String in timespan shorthand for size of training data to use.
            Ex: "8W"
    test_len_in : 
        String in timespan shorthand for size of testing data to use.
            Ex: "2W"
    test_date_step_in : 
        String in timespan shorthand for the time step separating one
            experiment's timeframe from the next experiment. Recommended
            default is to set this equal to test_len_in, so that testing
            periods are non-overlapping; can set this to None for this option.
            Exs: "3D" or None
    coverage_bounds_in : 
        Comma-separated string of decimals representing coverage percentages
            at which we would like to evaluate the models.
            Ex: "0.01,0.02,0.05,0.10"
    geojson_file_name_in : 
        File name of relevant geojson file. Should be located in datadir.
            Exs: "Chicago_South_Side_2790.geojson" or "Durham_27700.geojson"
    models_to_run_in : 
        Comma-separated string of names of models to run and evaluate.
            Currently recognised models are: random, naive, phs, ideal.
            Ex: "random,naive,phs,ideal"
    num_random_in : 
        Integer, or string that will be cast to an integer, representing
            the number of different times to run the random model.
            Ex: 3
    phs_time_units_in : 
        Comma-separated string of values for the atomic unit of time to be
            used for the PHS model, in timespan shorthand format.
            Ex: "1W"
    phs_time_bands_in : 
        Comma-separated string of values for the time bandwidths to test
            for the PHS model, in timespan shorthand format.
            Ex: "1W,2W,3W,4W,5W,6W,7W,8W"
    phs_dist_units_in : 
        Comma-separated string of values for the atomic unit of distance to be
            used for the PHS model, in meters.
            Ex: "100"
    phs_dist_bands_in : 
        Comma-separated string of values for the time bandwidths to test
            for the PHS model, in meters.
            Ex: "100,200,300,400,500,600,700,800,900,1000"
    phs_weight_in : 
        Comma-separated string of methods of calculating weight for PHS.
            Current recognised methods: classic, linear.
            classic = 
            linear = 
            Ex: "classic"
    print_exp_freq_in : 
        Integer, or string to be cast to an integer, representing how
            frequently some information about an experiment should be
            sent to stdout to help monitor the script's progress.
"""
def runModelExperiments(
            datadir_in, 
            dataset_name_in, 
            crime_type_set_in, 
            cell_width_in, 
            in_csv_file_name_in, 
            geojson_file_name_in, 
            local_epsg_in, 
            earliest_test_date_in, 
            test_date_range_in, 
            train_len_in, 
            test_len_in, 
            test_date_step_in, 
            coverage_bounds_in, 
            models_to_run_in, 
            coverage_max_in = None, 
            num_random_in = None, 
            phs_time_units_in = None, 
            phs_time_bands_in = None, 
            phs_dist_units_in = None, 
            phs_dist_bands_in = None, 
            phs_weight_in = None, 
            print_exp_freq_in = 1, 
            csv_date_format = "%m/%d/%Y %I:%M:%S %p", 
            csv_longlat = False, 
            csv_epsg = None, 
            csv_infeet = True, 
            csv_has_header = True, 
            ):
    
    
    
    
    # Variables with "chktime" are used to start checking the timing
    # Variables with "tkntime" hold the amount of time taken
    
    # Overall timing
    #chktime_overall = time.time()
    
    #### Declare data parameters
    #chktime_decparam = time.time()
    print("Declaring parameters...")
    

    
    
    
    ###
    # Parameters directly from input, recast as appropriate data types
    
    datadir = os.path.normpath(datadir_in)
    datadir = os.path.expanduser(datadir)
    datadir = os.path.expandvars(datadir)
    print("datadir is defined as:")
    print(datadir)
    if not os.path.isdir(datadir):
        print("Error!")
        print(f"Directory does not exist: {datadir}")
        print("Exiting...")
        sys.exit(1)
    dataset_name = dataset_name_in
    crime_type_set = set(splitCommaArgs(crime_type_set_in))
    cell_width = int(cell_width_in)
    in_csv_file_name = in_csv_file_name_in
    geojson_file_name = geojson_file_name_in
    earliest_test_date = earliest_test_date_in
    test_date_range = test_date_range_in
    train_len = train_len_in
    test_len = test_len_in
    test_date_step = test_date_step_in
    if test_date_step == None:
        test_date_step = test_len
    coverage_bounds = sorted([float(x) for x in \
                       splitCommaArgs(coverage_bounds_in)])
    models_to_run = splitCommaArgs(models_to_run_in)
    if "random" in models_to_run:
        if num_random_in == None:
            num_random = 1
        else:
            num_random = int(num_random_in)
    if "phs" in models_to_run:
        phs_time_units = splitCommaArgs(phs_time_units_in)
        phs_time_bands = splitCommaArgs(phs_time_bands_in)
        phs_dist_units = [int(x) for x in splitCommaArgs(phs_dist_units_in)]
        for du in phs_dist_units:
            if du > cell_width:
                print("PHS distance unit ({du}) cannot exceed cell width ({cell_width})")
                sys.exit(1)
        phs_dist_bands = [int(x) for x in splitCommaArgs(phs_dist_bands_in)]
        phs_weight = splitCommaArgs(phs_weight_in)
    print_exp_freq = int(print_exp_freq_in)
    coverage_max = coverage_bounds[-1]
    if coverage_max_in != None:
        coverage_max = float(coverage_max_in)
    
    
    
    
    
    
    
    
    
    
    
    ###
    # Derived parameters
    
    
    # Mappings of model names to all desired parameter combinations
    # Ideal and Naive models have no additional parameters
    # Random model's only extra parameter is number of times to run
    # PHS has several parameters: time units, time bandwidths, distance units,
    #  distance bandwidths, and weight method
    model_param_dict = dict()
    # Keep separate list of short string representations of parameter
    #  combinations, for ease of printing on graphs, e.g.
    model_params_short = defaultdict(list)
    if "ideal" in models_to_run:
        model_param_dict["ideal"] = [()]
    if "naive" in models_to_run:
        model_param_dict["naive"] = [()]
    # Param list for Random model.
    #  For example, if 4, param list is [(0,), (1,), (2,), (3,)]
    if "random" in models_to_run:
        model_param_dict["random"] = list(product(range(num_random)))
        if num_random > 1:
            model_params_short["random"]=[str(x) for x in range(num_random)]
    # Param list for PHS model.
    if "phs" in models_to_run:
        poss_phs_params = [phs_time_units, 
                           phs_time_bands, 
                           phs_dist_units, 
                           phs_dist_bands, 
                           phs_weight]
        model_param_dict["phs"] = list(product(*poss_phs_params))
        params_to_iter = []
        for p_list in poss_phs_params:
            if len(p_list)>1:
                params_to_iter.append([str(x) for x in p_list])
        model_params_short["phs"] = \
                        ["-".join(p) for p in product(*params_to_iter)]
    model_idents = dict()
    for model_name in models_to_run:
        if len(model_params_short[model_name])<2:
            model_idents[model_name] = [model_name]
        else:
            model_idents[model_name] = \
                [f"{model_name}-{x}" for x in model_params_short[model_name]]
    #model_idents_flat = [iden for mname in model_idents for iden in dd[mname]]
    
    
    
    # Printable string of crime types, concatenated by "_" if necessary
    crime_type_set_dashed = ["-".join(x.split()) for x in sorted(crime_type_set)]
    crime_types_printable = "_".join(sorted(crime_type_set_dashed))
    # Full path for input csv file
    in_csv_full_path = os.path.join(datadir, in_csv_file_name)
    # Nicely-formatted string of test date
    earliest_test_date_str = "".join(earliest_test_date.split("-"))[2:]
    # Latest start of a test data set, calculated from earliest and length
    latest_test_date = generateLaterDate(earliest_test_date, test_date_range)
    # List of all experiment dates
    start_test_list = generateDateRange(start=earliest_test_date, 
                                        end=latest_test_date, 
                                        step=test_date_step)
    # Number of different experiment dates
    total_num_exp_dates = len(start_test_list)
    # If number of experiment dates <= 2, declare this run to be short.
    # This means that for each experiment date, we will create various plots:
    #   - training locations
    #   - test locations
    #   - heat maps for each experiment
    #   - coverage maps for each experiment
    # In this case we might also save off the risk scores computed for each cell
    #  by each experiment, and relative rankings of those cells.
    #  (But that part isn't implemented yet.)
    run_is_short = False
    if total_num_exp_dates <= 2:
        run_is_short = True
    print(f"Number of experiments: {total_num_exp_dates}")
    # String that uniquely identifies this run within a file name
    file_name_core = "_".join([date_today_str, \
                               dataset_name, \
                               earliest_test_date_str, \
                               test_date_range, \
                               test_date_step, \
                               train_len, \
                               test_len])
    # Output csv file name for results summary
    out_csv_file_name_results = f"results_{file_name_core}.csv"
    # Output csv file name for detailed risk info if run is short
    #out_csv_file_name_risks = f"risks_{file_name_core}.csv"
    # Full path for output csv file of results
    out_csv_results_full_path = os.path.join(datadir, out_csv_file_name_results)
    # Full path for output csv file of risk info if run is short
    #out_csv_risks_full_path = os.path.join(datadir, out_csv_file_name_risks)
    # Output geojson training data file
    out_train_geojson_name = f"train_{file_name_core}.geojson"
    # Full path for output geojson training data file
    out_train_geojson_full_path = os.path.join(datadir, out_train_geojson_name)
    # Output geojson testing data file
    out_test_geojson_name = f"test_{file_name_core}.geojson"
    # Full path for output geojson testing data file
    out_test_geojson_full_path = os.path.join(datadir, out_test_geojson_name)
    # Output geojson results file
    out_res_geojson_name = f"results_{file_name_core}.geojson"
    # Full path for output geojson results file
    out_res_geojson_full_path = os.path.join(datadir, out_res_geojson_name)
    # Full path for input geojson file
    in_geojson_full_path = os.path.join(datadir, geojson_file_name)
    # Append risk_MODEL and rank_MODEL to risk_info_header for each model
    if run_is_short:
        for model_name in models_to_run:
            risk_info_header.append(f"risk_{model_name}")
            risk_info_header.append(f"rank_{model_name}")
    
    # Check that data file and geojson file exist
    if not os.path.isfile(in_csv_full_path):
        print("Error!")
        print(f"File does not exist: {in_csv_full_path}")
        print("Exiting...")
        sys.exit(1)
    if not os.path.isfile(in_geojson_full_path):
        print("Error!")
        print(f"File does not exist: {in_geojson_full_path}")
        print("Exiting...")
        sys.exit(1)
    
    
    print("...declared parameters.")
    #tkntime_decparam = time.time() - chktime_decparam
    
    
    
    
    # Open output csv file for writing, write header row
    with open(out_csv_results_full_path, "w") as csvf:
        results_writer = csv.writer(csvf, delimiter=",", lineterminator="\n")
        results_writer.writerow(result_info_header)
        
        
        
        
        
        # If we were to do a "crime type sweep", that would go here.
        # But probably simpler to just re-run a new crime type set instead.
        # Unless we want to actively combine results within this script?
        
        
        
        
        
        ###
        # Obtain input data
        
        chktime_obtain_data = time.time()
        print("Obtaining full data set and region...")
        
        #!!! Need to change standardization pre-processing so that this input
        #     is always Eastings/Northings, and in meters
        points_crime = loadGenericData(in_csv_full_path, 
                                       crime_type_set=crime_type_set, 
                                       date_format_csv = csv_date_format, 
                                       longlat=csv_longlat, 
                                       epsg = csv_epsg, 
                                       infeet=csv_infeet, 
                                       has_header = csv_has_header)
        
        num_crimes_total = len(points_crime.timestamps)
        print(f"Total number of relevant crimes: {num_crimes_total}")
        
        
        # Obtain polygon from geojson file (which should have been pre-processed)
        region_polygon = gpd.read_file(in_geojson_full_path)
        # Convert to relevant CRS for local projection
        
        region_polygon = region_polygon.to_crs({'init': f'epsg:{local_epsg_in}'})
        # Take unary union, which also converts region from being
        #  a GeoDataFrame to a Polygon
        region_polygon = region_polygon.unary_union
        
        # Get subset of input crime that occurred within region
        points_crime_region = intersect_timed_points(points_crime, region_polygon)
        
        
        
        
        num_crimes_region = len(points_crime_region.timestamps)
        print(f"Total number of relevant crimes in area: {num_crimes_region}")
        
        
        
        # Get grid version of region
        masked_grid_region = mask_grid_by_intersection(region_polygon, 
                                                       DataGrid(xsize=cell_width, 
                                                                ysize=cell_width, 
                                                                xoffset=0, 
                                                                yoffset=0)
                                                       )
        
        
        
        # Create GeoDataFrame of relevant cells, where each entry has
        # an associated rectangular geometry corresponding to a cell.
        gdf_cells = gdt.make_cells_frame(
                                        masked_grid_region, 
                                        27700
                                        )
        
        
        
        
        
        # Get "mesh info" of that grid, which is useful for displaying the map.
        #  masked_grid_region is type open_cp.data.MaskedGrid
        #  masked_grid_mesh is type tuple, with 2 elements
        #   1st element is a list of x-coordinates (in Eastings) for grid
        #   2nd element is a list of y-coordinates (in Northings) for grid
        masked_grid_mesh = masked_grid_region.mesh_info()
        
        
        
        
        # Get tuple of all cells in gridded region
        # Each cell is represented as a 2-tuple of 0-up indices for the cells
        #  within the (masked) grid. That is, they are NOT lat/long or E/N
        #  coordinates, but they can be mapped to those coordinates by
        #  using masked_grid_mesh, for example
        cellcoordlist_region = getRegionCells(masked_grid_region)
        
        # Obtain number of cells in the grid that contain relevant geometry
        # (i.e., not the full rectangular grid, only relevant cells)
        num_cells_region = len(cellcoordlist_region)
        
        print("...obtained full data set and region.")
        tkntime_obtain_data = time.time() - chktime_obtain_data
        print(f'Time taken to obtain data: {tkntime_obtain_data}')
        
        
        
        
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
            points_crime_region_train = getTimedPointsInTimeRange(
                                            points_crime_region, 
                                            start_train, 
                                            end_train)
            # Count how many crimes there were in this training data set
            num_crimes_train = len(points_crime_region_train.timestamps)
            
            
            gdf_datapoints_train = \
                        gdt.make_points_frame(points_crime_region_train, 
                                              csv_epsg)
            gdf_datapoints_train.to_file(out_train_geojson_full_path,
                                   driver='GeoJSON')
            
            
            
            # Obtain testing data
            points_crime_region_test = getTimedPointsInTimeRange(
                                            points_crime_region, 
                                            start_test, 
                                            end_test)
            # Count how many crimes there were in this test data set
            num_crimes_test = len(points_crime_region_test.timestamps)
            
            # Count the number of crimes per cell in test data.
            #  This is used for evaluation.
            cells_testcrime_ctr = countPointsPerCell(points_crime_region_test, 
                                                     masked_grid_region)
            
            
            
            gdf_datapoints_test = \
                        gdt.make_points_frame(points_crime_region_test, 
                                              csv_epsg)
            gdf_datapoints_test.to_file(out_test_geojson_full_path,
                                   driver='GeoJSON')
            
            
            # Display number of crimes in training and testing data
            if exp_date_index % print_exp_freq == 0:
                print(f"num_crimes_train: {num_crimes_train}")
                print(f"num_crimes_test: {num_crimes_test}")
            
            
            
            # If we have few experiments (i.e. 1 date), then create and save
            #  map visualisations of the data here. Later we'll also save
            #  visualisations of the models' results.
            if run_is_short:
                
                print("Making image file names for training and testing data")
                
                # Make image file names for training and testing data
                
                img_file_core = "_".join([
                                        date_today_str, 
                                        getSixDigitDate(start_test), 
                                        train_len, 
                                        test_len, 
                                        ])
                
                img_file_train_name = f"trainmap_{img_file_core}.png"
                img_file_test_name = f"testmap_{img_file_core}.png"
                img_file_train_full_path = os.path.join(datadir, 
                                                        img_file_train_name)
                img_file_test_full_path = os.path.join(datadir, 
                                                       img_file_test_name)
                print(f"Training: {img_file_train_name}")
                print(f"Testing: {img_file_test_name}")
                
                # Plot training data on plain grid
                plotPointsOnGrid(points_crime_region_train, 
                                 masked_grid_region, 
                                 region_polygon, 
                                 title=f"Train data {exp_date_index}", 
                                 sizex=10, 
                                 sizey=10, 
                                 out_img_file_path=img_file_train_full_path)
                
                # Plot testing data on plain grid
                plotPointsOnGrid(points_crime_region_test, 
                                 masked_grid_region, 
                                 region_polygon, 
                                 title=f"Test data {exp_date_index}", 
                                 sizex=10, 
                                 sizey=10, 
                                 out_img_file_path=img_file_test_full_path)
                
                
                
                
            
            
            
            # A "data_matrix" contains a risk score for each cell in the 
            #  grid (within the masked region). These scores can then be 
            #  used to rank the cells by risk. Note that different models 
            #  use different scoring systems, so the scores from different 
            #  data_matrices are not directly comparable, unless some 
            #  normalisation process is used.
            
            # A "sorted_cells" list is a ranked list of cells based on the
            #  previously computed data_matrix. Ties are currently broken by
            #  selecting southernmost cells, then westernmost cells.
            
            # A "rank_matrix" is a matrix object where each cell is associated
            #  with its ranking from the sorted_cells list. This is useful for
            #  displaying a map of which cells are covered by various coverage
            #  thresholds.
            
            
            
            # Result objects for various experiments
            # model_name -> exp_num
            data_matrix_dict = defaultdict(list)
            sorted_cells_dict = defaultdict(list)
            rank_matrix_dict = defaultdict(list)
            hit_count_list_dict = defaultdict(list)
            hit_rate_list_dict = defaultdict(list)
            
            
            
            
            
            for model_name in models_to_run:
                if model_name not in recognised_models:
                    print("Error!")
                    print(f"Unrecognised model name: {model_name}")
                    print(f"Recognised models are: {recognised_models}")
                    print("Skipping that model.")
                    continue
                
                model_params = model_param_dict[model_name]
                
                # Generate the data matrix based on the particular model
                if model_name == "random":
                    for rand_seed in [x[0] for x in model_params]:
                        data_matrix_dict[model_name].append(deepcopy(
                                runRandomModel(
                                        grid=masked_grid_region, 
                                        rand_seed = rand_seed
                                        )
                                ))
                elif model_name == "naive":
                    data_matrix_dict[model_name].append(deepcopy(runNaiveModel(
                            training_data=points_crime_region_train, 
                            grid=masked_grid_region
                            )))
                elif model_name == "phs":
                    for params_combo in model_params:
                        # Cast PHS parameters into proper data types
                        phs_time_unit = shorthandToTimeDelta(params_combo[0])
                        phs_time_band = shorthandToTimeDelta(params_combo[1])
                        phs_dist_unit = int(params_combo[2])
                        phs_dist_band = int(params_combo[3])
                        phs_weight = params_combo[4]
                        
                        data_matrix_dict[model_name].append(deepcopy(runPhsModel(
                                training_data=points_crime_region_train, 
                                grid=masked_grid_region, 
                                cutoff_time=start_test, 
                                time_unit=phs_time_unit, 
                                dist_unit=phs_dist_unit, 
                                time_bandwidth=phs_time_band, 
                                dist_bandwidth=phs_dist_band, 
                                weight=phs_weight
                                )))
                elif model_name == "ideal":
                    data_matrix_dict[model_name].append(deepcopy(runNaiveModel(
                            training_data=points_crime_region_test, 
                            grid=masked_grid_region
                            )))
                else:
                    print("Error!")
                    print(f"This model has not been implemented: {model_name}")
                    print("Skipping for now...")
                    continue
                
                
                
                
                
                
                for exp_index, data_matrix in enumerate(
                                                data_matrix_dict[model_name]):
                    
                    if exp_date_index % print_exp_freq == 0:
                        print(f"model_name: {model_name}")
                        print(f"parameter set #: {exp_index}")
                    sorted_cells_dict[model_name].append(deepcopy(
                            sortCellsByRiskMatrix(
                                cellcoordlist_region, 
                                data_matrix)
                            ))
                    
                    
                    hit_count_list_dict[model_name].append(
                        deepcopy(getHitCountList(
                            sorted_cells_dict[model_name][exp_index], 
                            cells_testcrime_ctr))
                        )
                    
                    if num_crimes_test == 0:
                        hit_rate_list_dict[model_name].append(deepcopy(
                            [0 for x in hit_count_list_dict[model_name][-1]]
                        ))
                    else:
                        hit_rate_list_dict[model_name].append(deepcopy(
                            [float(x/num_crimes_test) \
                                 for x in hit_count_list_dict[model_name][-1]]
                        ))
                        
                    
                    
                    
                    for coverage_rate in coverage_bounds:
                        
                        # Get number and % of hits in results
                        #results_hit, results_pct = hitRatesFromHitList(
                        #            hit_count_list_dict[model_name][exp_index], 
                        #            coverage_rate, 
                        #            num_cells_region, 
                        #            num_crimes_test)
                        
                        
                        cov_rate_index = int(coverage_rate * num_cells_region)
                        results_hit = hit_count_list_dict[model_name][exp_index][cov_rate_index]
                        results_pct = hit_rate_list_dict[model_name][exp_index][cov_rate_index]
                        
                        # Standard result info to include for every model
                        result_info = [
                                dataset_name, 
                                crime_types_printable, 
                                cell_width, 
                                start_test, 
                                train_len, 
                                test_len, 
                                coverage_rate, 
                                num_crimes_test, 
                                results_hit, 
                                results_pct, 
                                model_name]
                        # Pre-pad PHS parameters with empty spaces where
                        #  columns for non-PHS parameters are
                        if model_name == "phs":
                            result_info += ["",""]
                        # 
                        result_info += list(model_param_dict[model_name][exp_index])
                        
                        # Pad out rest of row with empty string values
                        while len(result_info) < len(result_info_header):
                            result_info.append("")
                        
                        # Write row to csv file
                        results_writer.writerow(result_info)
                        
                        
                    if run_is_short:
                        print("Plotting 2 maps for each model")
                        # Plot 2 maps for each model
                        # Create rank matrix
                        new_rank_matrix = rankMatrixFromSortedCells(
                                    masked_grid_region, 
                                    sorted_cells_dict[model_name][exp_index], 
                                    coverage_bounds
                                )
                        # Store rank matrix
                        rank_matrix_dict[model_name].append(deepcopy(
                                new_rank_matrix
                                ))
                        
                        
                        
                        
                        # Store data matrix values in GeoDataFrame
                        new_exp_ident = model_idents[model_name][exp_index]
                        new_data_matrix = \
                            data_matrix_dict[model_name][exp_index]
                        gdf_cells[new_exp_ident + "-score"] = \
                            [new_data_matrix[c] for c in cellcoordlist_region]
                        gdf_cells[new_exp_ident + "-rank"] = \
                            [new_rank_matrix[c] for c in cellcoordlist_region]
                        
                        
                        
                        # Save heat map and coverage map as image files
                        saveModelResultMaps(
                            model_name, 
                            data_matrix_dict[model_name][exp_index], 
                            rank_matrix_dict[model_name][exp_index], 
                            exp_ident = model_idents[model_name][exp_index], 
                            file_core = img_file_core, 
                            filedir = datadir, 
                            polygon = region_polygon, 
                            points_to_map = points_crime_region_train, 
                            mesh_info = masked_grid_mesh, 
                            )
                
                
            
            
            # If only a couple experiments, then create a line graph
            #  comparing hit rates of different models (y-axis) based
            #  on coverage rate (x-axis)
            if run_is_short:
                
                # Save GeoDataFrame to file
                print("Saving results to geojson:")
                print(out_res_geojson_full_path)
                gdt.frame_to_json_with_id(gdf_cells, 
                                          out_res_geojson_full_path)
                
                
                print("Plotting graph coverage vs hit rate")
                
                
                
                img_file_core = "_".join([
                        date_today_str, 
                        getSixDigitDate(start_test), 
                        train_len, 
                        test_len, 
                        ])
                
                img_file_graph_name = f"hitrates_{img_file_core}.png"
                img_file_graph_full_path = os.path.join(datadir, 
                                                        img_file_graph_name)
                
                
                
                graphCoverageVsHitRate(
                        hit_count_list_dict, 
                        model_idents, 
                        models_to_run, 
                        x_limits=(0, coverage_max), 
                        title="Hit rate evaluation of models by coverage", 
                        out_img_file_path = img_file_graph_full_path)
                #graphCoverageVsHitRate(hit_rate_list_dict, 
                #                       model_idents, 
                #                       models_to_run)
            
            
            
            
            
            
            tkn_time_exp = time.time()-chktime_exp
            print(f"time spent on exp: {tkn_time_exp}")
            exp_times.append(tkn_time_exp)
            
        
        
        
        
        
        
    
    
    if not run_is_short:
        print("Plotting graph hit rates over time")
        from riskModelsResultsEval import graphHitRatesOverTime
        graphHitRatesOverTime(out_csv_results_full_path)
    
    
    print("Experiment timing info:")
    print("Exp #\tTime")
    for i, t in enumerate(exp_times):
        print(f"{i}\t{t}")
    return

"""
# if run_is_short then write risk_info file at the end here
if run_is_short:
    # Open risk_info csv file for writing, write header row
    with open(out_csv_risks_full_path, "w") as csvrisk:
        risk_writer = csv.writer(csvrisk, delimiter=",", lineterminator="\n")
        risk_writer.writerow(risk_info_header)
        
        # Need to have stored info to write here
    
    

"""




"""
main:

If running this module as a script instead of importing its functions,
 this main function will perform runModelExperiments() with a set of
 default parameters.
"""
def main():
    # Run evaluation function with default arguments
    
    
    
    
    
    # Location of data file
    datadir = "../../Data"
    # Dataset name (to be included in name of output file)
    #dataset_name = "chicago"
    dataset_name = "FantDur"
    # Crime types
    #!!!  May want to change this to be different for train vs test?
    #!!!  May want to use multiple sets of types, then combine results?
    #crime_type_set = "BURGLARY"
    crime_type_set = "Burglary, Vehicle crime"
    #crime_type_set_sweep = [{"BURGLARY"}] # if doing a sweep over different sets
    # Size of grid cells
    #cell_width_sweep = [100] # if doing a sweep over different cell sizes
    cell_width = 500
    # Input csv file name
    #in_csv_file_name = "chi_all_s_BURGLARY_RES_010101_190101_stdXY.csv"
    in_csv_file_name = "Fantasy-Durham-Data_std.csv"
    # Geojson file
    #geojson_file_name = "Chicago_Areas.geojson"
    #geojson_file_name = "Chicago_South_Side_2790.geojson"
    #geojson_file_name = "Durham_27700.geojson"
    geojson_file_name = "Police_Force_Areas_December_2016_Durham_fixed.geojson"
    local_epsg = 27700
    # Of all planned experiments, earliest start of a TEST (not train) data set
    #earliest_test_date = "2016-09-01"
    earliest_test_date = "2019-09-01"
    # Time between earliest experiment and latest experiment
    test_date_range = "1W"
    # Length of training data
    train_len = "4W"
    # Length of testing data
    test_len = "1W"
    # Time step offset between different experiments
    # (If you want non-overlapping, then set test_date_step = test_len)
    #test_date_step = "1D"
    test_date_step = None
    # Coverage rates to test
    coverage_bounds = "0.01,0.02,0.05,0.10"
    
    # Highest coverage we want to plot on line graph, if we make one
    coverage_max = 0.1
    
    
    # Predictive models to run
    #models_to_run = "random,naive,phs,ideal"
    models_to_run = "random,naive,phs,ideal"
    
    num_random = 1
    
    # Param list for PHS model.
    phs_time_units = "1W"
    phs_time_bands = "4W,6W"
    phs_dist_units = "500"
    phs_dist_bands = "500,1000,1500"
    phs_weight = "classic"
    
    
    # How frequently to display a print statement about an experiment
    print_exp_freq = 1
    
    
    
    
    if dataset_name == "chicago":
        """
        crime_type_set = "BURGLARY"
        in_csv_file_name = "chi_all_s_BURGLARY_RES_010101_190101_stdXY.csv"
        geojson_file_name = "Chicago_South_Side_2790.geojson"
        earliest_test_date = "2016-09-01"
        """
        local_epsg = 2790
        csv_date_format = "%m/%d/%Y %I:%M:%S %p"
        csv_longlat = False
        csv_epsg = None
        csv_infeet = True
        csv_has_header = True
    
    if dataset_name == "FantDur":
        """
        crime_type_set = "Burglary, Vehicle crime"
        in_csv_file_name = "Fantasy-Durham-Data_std.csv"
        geojson_file_name = "Durham_27700.geojson"
        earliest_test_date = "2019-09-01"
        """
        
        local_epsg = 27700
        csv_date_format = "%d/%m/%Y"
        csv_longlat = True
        csv_epsg = 27700
        csv_infeet = False
        csv_has_header = True
    
    
    
    
    
    runModelExperiments(
            datadir_in = datadir, 
            dataset_name_in = dataset_name, 
            crime_type_set_in = crime_type_set, 
            cell_width_in = cell_width, 
            in_csv_file_name_in = in_csv_file_name, 
            geojson_file_name_in = geojson_file_name, 
            local_epsg_in = local_epsg, 
            earliest_test_date_in = earliest_test_date, 
            test_date_range_in = test_date_range, 
            train_len_in = train_len, 
            test_len_in = test_len, 
            test_date_step_in = test_date_step, 
            coverage_bounds_in = coverage_bounds, 
            models_to_run_in = models_to_run, 
            coverage_max_in = coverage_max, 
            num_random_in = num_random, 
            phs_time_units_in = phs_time_units, 
            phs_time_bands_in = phs_time_bands, 
            phs_dist_units_in = phs_dist_units, 
            phs_dist_bands_in = phs_dist_bands, 
            phs_weight_in = phs_weight, 
            print_exp_freq_in = print_exp_freq, 
            csv_date_format = csv_date_format, 
            csv_longlat = csv_longlat, 
            csv_epsg = csv_epsg, 
            csv_infeet = csv_infeet, 
            csv_has_header = csv_has_header, 
            )
    
    
    
    




if __name__ == "__main__":
    main()
