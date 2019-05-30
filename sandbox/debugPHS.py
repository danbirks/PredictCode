# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 14:37:28 2019

@author: Dustin
"""



# Some fairly standard modules
import os, csv, lzma
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import descartes
from itertools import product
from collections import Counter
import datetime
import csv
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

# Load custom functions that make dealing with datetime and timedelta easier
from crimeRiskTimeTools import generateDateRange, \
                               generateLaterDate, \
                               generateEarlierDate, \
                               getTimedPointsInTimeRange, \
                               _day

# Useful line for python console
# sys.path.insert(0,os.path.join(os.path.abspath("."),"Documents","GitHub","PredictCode","sandbox"))



#Constants for figures

# Heat-like color mapping
_cdict = {'red':   [(0.0,  1.0, 1.0),
                   (1.0,  1.0, 1.0)],
         'green': [(0.0,  1.0, 1.0),
                   (1.0,  0.0, 0.0)],
         'blue':  [(0.0,  0.2, 0.2),
                   (1.0,  0.2, 0.2)]}

yellow_to_red = matplotlib.colors.LinearSegmentedColormap("yellow_to_red", _cdict)




#Functions


# Functions that help process the data


# Given a TimedPoints object and a Grid (MaskedGrid?) object,
#  return a Counter object that is a mapping from the grid cell
#  coordinates to the number of points within the cell.
#  Note that "grid cell coordinates" refers to which row of cells
#  and which column of cells it's located at, NOT spatial coords.
def countPointsPerCell(points, grid):
    # Get xy coords from TimedPoints
    xcoords, ycoords = points.xcoords, points.ycoords
    # Convert coords to cellcoords
    xgridinds = np.floor((xcoords - grid.xoffset) / grid.xsize).astype(np.int)
    ygridinds = np.floor((ycoords - grid.yoffset) / grid.ysize).astype(np.int)
    # Count the number of crimes per cell
    # NOTE: We do (y,x) instead of (x,y) because cells are (row,col)!!!
    return Counter(zip(ygridinds, xgridinds))


# Given a sorted list of cells, and a mapping from cells to number of events
#  in those cells, return a list of numbers, of length equal to the given
#  list of cells, representing the running total of number of events in all
#  cells up to that point in the list.
def getHitRateList(sorted_cells, cell_hit_map):
    running_total = 0
    hit_rate_list = []
    for cell in sorted_cells:
        running_total += cell_hit_map[cell]
        hit_rate_list.append(running_total)
    return hit_rate_list


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



def getRegionCells(grid):
    # Make sure to do yextent then xextent, because cellcoords
    #  correspond to (row,col) in grid
    all_cells = product(range(grid.yextent), range(grid.xextent))
    return tuple(cc for cc in all_cells 
                 if not grid.mask[cc[0]][cc[1]])





# Functions that run models


# Construct a "random" model by simply randomly sorting the cells
def RandomlySortCells(cells, seed=None):
    if seed != None:
        random.seed(seed)
    cell_list = list(cells)
    random.shuffle(cell_list)
    return cell_list



# Most naive model is to just count the number of events occurring in each
#  cell in the training data, and favor the cells with the most events
def runNaiveCount_model(data_points, grid):
    crime_ctr = countPointsPerCell(data_points, grid)
    grid_cells = getRegionCells(grid)
    print(sorted(grid_cells)[0])
    print(sorted(grid_cells)[-1])
    print(max([x[0] for x in grid_cells]))
    print(max([x[1] for x in grid_cells]))
    return sorted(grid_cells, key=lambda x:crime_ctr[x], reverse=True)



def runRhsModel(training_data, grid, bandwidth=250, rand_seed=None):
    # Set RNG seed if given
    if rand_seed != None:
        np.random.seed(rand_seed)
    
    grid_cells = getRegionCells(grid)
    
    # Obtain model and prediction on grid cells
    rhs_pred = retro.RetroHotSpot()
    rhs_pred.data = training_data
    rhs_pred.weight = retro.Quartic(bandwidth = bandwidth)
    rhs_risk = rhs_pred.predict()
    rhs_grid_risk = open_cp.predictors.grid_prediction(rhs_risk, grid)
    rhs_grid_risk_matrix = rhs_grid_risk.intensity_matrix
    
    # Sort cellcoords by risk in intensity matrix, highest risk first
    return sortCellsByRiskMatrix(grid_cells, rhs_grid_risk_matrix)
















def knox_ratio(statistic, distribution):
    # From "Examples/Chicago Case Study/Knox Statistics" notebook
    # Compute the ratio of the statistic to the 
    #  median of the values in the distribution
    d = np.array(distribution)
    d.sort()
    return statistic / d[len(d)//2]



def all_knox_ratios(result):
    for i, space_bin in enumerate(result.space_bins):
        for j, time_bin in enumerate(result.time_bins):
            yield knox_ratio(result.statistic(i,j), result.distribution(i,j))









def runPhsModel(training_data, grid, time_unit, dist_unit, time_bandwidth, dist_bandwidth, weight="linear"):
    
    
    grid_cells = getRegionCells(grid=grid)
    
    # Obtain model and prediction on grid cells
    phs_predictor = phs.ProspectiveHotSpot(grid=grid)
    phs_predictor.data = training_data
    
    dist_band_in_units = dist_bandwidth/dist_unit
    time_band_in_units = time_bandwidth/time_unit
    
    if weight=="linear":
        phs_predictor.weight = phs.LinearWeightNormalised(space_bandwidth=dist_band_in_units, time_bandwidth=time_band_in_units)
    elif weight=="classic":
        phs_predictor.weight = phs.ClassicWeightNormalised(space_bandwidth=dist_band_in_units, time_bandwidth=time_band_in_units)
    
    phs_predictor.time_unit = time_unit
    phs_predictor.grid = dist_unit
    
    cutoff_time = sorted(training_data.timestamps)[-1] + _day
    
    phs_grid_risk = phs_predictor.predict(cutoff_time, cutoff_time)
    
    
    
    #phs_grid_risk = open_cp.predictors.grid_prediction(phs_risk, grid)
    phs_grid_risk_matrix = phs_grid_risk.intensity_matrix
    print(type(phs_grid_risk))
    print(type(phs_grid_risk_matrix))
    print(phs_grid_risk_matrix)
    print(phs_grid_risk_matrix.shape)
    sys.exit(0)
    
    # Sort cellcoords by risk in intensity matrix, highest risk first
    return sortCellsByRiskMatrix(grid_cells, phs_grid_risk_matrix)
    








# Given a model name and relevant arguments,
#  return a sorted list of cells
def runModelAndSortCells(model_name, model_args):
    
    # We declare our recognised possible models here
    rec_models = ["ideal","random","naivecount","rhs","phs"]
    if model_name not in rec_models:
        print("Unrecognized model name: {}".format(model_name))
        sys.exit(1)
    
    
    
    if model_name=="ideal":
        # We need these variables:
        #   cellcoordlist_region
        #   cells_testcrime_ctr
        
        cellcoordlist_region, cells_testcrime_ctr = model_args
        
        return sorted(cellcoordlist_region, 
                                         key=lambda x:cells_testcrime_ctr[x], 
                                         reverse=True)
        
    if model_name=="random":
        # We need these variables:
        #   cellcoordlist_region
        #   plot_random_seed
        
        cellcoordlist_region, plot_random_seed = model_args
        
        return RandomlySortCells(cellcoordlist_region, seed=plot_random_seed)
        
        
    
    # If the model isn't ideal or random,
    #  then it's naivecount or rhs or phs,
    #  so the first two args should be the data and the region
    
    points_crime_region_train = model_args[0]
    masked_grid_region = model_args[1]
    other_model_args = model_args[2:]
    
    if model_name=="naivecount":
        # We need these variables:
        #   points_crime_region_train
        #   masked_grid_region
        
        return runNaiveCount_model(points_crime_region_train, masked_grid_region)
        
    
    
    if model_name=="rhs":
        # We need these variables:
        #   points_crime_region_train
        #   masked_grid_region
        #   rhs_bandwidth
        #   rhs_random_seed
        
        rhs_random_seed, rhs_bandwidth = other_model_args
        
        return runRhsModel(points_crime_region_train,
                                                   masked_grid_region, 
                                                   bandwidth = rhs_bandwidth, 
                                                   rand_seed=rhs_random_seed)
    
    if model_name=="phs":
        # We need these variables:
        #   points_crime_region_train
        #   masked_grid_region
        #   time_unit (ex: np.timedelta64(1, "W") )
        #   time_bandwidth (ex: np.timedelta64(4, "W") )
        #   dist_unit (ex: 100 )
        #   dist_bandwidth (ex: 500 )
        #   weight (ex: "linear", "classic" )
        
        time_unit, time_bandwidth, dist_unit, dist_bandwidth, weight = other_model_args
        return runPhsModel(
                training_data = points_crime_region_train, 
                grid = masked_grid_region, 
                time_unit = time_unit, 
                dist_unit = dist_unit, 
                time_bandwidth = time_bandwidth, 
                dist_bandwidth = dist_bandwidth, 
                weight = weight)
    





# Start of code



# START TIMER
init_start_time = time.time()



# PARAMETERS




# PARAMETERS TO SWEEP
# Model-invariant parameters:
# > DATA
#   - cell width (assume height is same)
#   - overall data set (currently hard-coded as chicago)
#   - crime type(s), just doing BURGLARY for now
#   - spatial offset? (not doing this now, but maybe worth investigating)
# > DATE
#    - Length of training window
#    - Length of testing window
#    - Date (of, say, testing; other dates calculated from lengths)
# > MODEL
#   - model type(s)
# > EVAL
#   - Coverage (1%, 2%, 5%, 10%)
#   - Hit count or hit % as metric for evaluation
# Model-specific parameters:
# > Random
#   - seed
# > RHS
#   - rhs_bandwidth
#   - seed
# > PHS
#   - time_unit
#   - dist_unit (==cell width, probably?)
#   - time_bandwidth
#   - dist_bandwidth
#   - choice of weight? (linear vs classic)
#   - knox stuff? (some relates to the above)
#   - 


# What do I want to see in output?
# ~~~INPUT DATA
# Data set (e.g. "Chicago South Side")
# Crime type(s) (e.g. "BURGLARY", or "BURGLARY,THEFT")
# Cell width
# Test start date (=train end date)
# Train len
# ~~~EVAL RESULTS
# Test len
# Coverage rate
# Num crimes in test data
# Num crimes in test data captured by model within coverage rate
# % crimes in test data captured by model within coverage rate
# ~~~MODEL
# Model type
# (various model-specific parameters: random seeds, bandwidths, etc, etc)




#models_to_run = ["random", "naivecount", "rhs", "phs", "ideal"]
#models_to_run = ["random", "naivecount", "phs", "ideal"]
models_to_run = ["naivecount", "phs"]
#models_to_run = ["phs"]

model_param_dict = dict()

model_param_dict["ideal"] = [()]

model_param_dict["naivecount"] = [()]

num_random = 1
rand_seeds = range(num_random)
model_param_dict["random"] = list(product(rand_seeds))

num_rhs = 1
rhs_seeds = range(num_rhs)
rhs_bandwidth_sweep = [300]
model_param_dict["rhs"] = list(product(rhs_seeds, rhs_bandwidth_sweep))

phs_time_units = [np.timedelta64(1, "W")]
#phs_time_bands = [np.timedelta64(4, "W")]
#phs_time_bands = [np.timedelta64(x, "W") for x in range(1,7)]
phs_time_bands = [np.timedelta64(x, "W") for x in range(1)]
phs_dist_units = [100]
#phs_dist_bands = [500]
#phs_dist_bands = [x*100 for x in range(1,11)]
phs_dist_bands = [x*100 for x in range(4,5)]
phs_weights = ["linear"]
model_param_dict["phs"] = list(product(
                            phs_time_units, 
                            phs_time_bands, 
                            phs_dist_units, 
                            phs_dist_bands, 
                            phs_weights))






# Parameters for overall data set
dataset_name = "Chicago"
crime_type_set_sweep = [{"BURGLARY"}]
cell_width_sweep = [200]
# If we did spatial offsets, that would belong here too
# Also if there's a convenient way to specify Chicago vs other data set, do that here

# Parameters for time range







# Data parameters
print("Declaring parameters...")
datadir = os.path.join("..", "..", "Data")
chicago_file_name = "chicago_all_old.csv"
chicago_side = "South"
chicago_load_type = "snapshot"
if "all" in chicago_file_name:
    chicago_load_type = "all"
chicago_file_path = os.path.join(datadir, chicago_file_name)
# Chicago module requires this line to access some data
chicago.set_data_directory(datadir)




crime_type_set = {"BURGLARY"}




# Spatial model parameters
#cell_width = 250
#cell_height = cell_width


# Time parameters
# The best single date to define an experiment by is the start of the test
#  data. The training data will be from a given time range fully up to but
#  not including that date, while the test data will be from a given time
#  range starting on that date. If we wish to compare different sizes of
#  training or test data, then the best comparison would be against the
#  other experiments with this same date as the cutoff between the training
#  and testing data, regardless of the sizes of those data sets.

# Of all planned experiments, earliest start of a test data set
earliest_test_date = "2002-01-01"
earliest_test_date_str = "".join(earliest_test_date.split("-"))[2:]
# Time between earliest experiment and latest experiment
test_date_range = "1W"
# Latest start of a test data set, calculated from above 2 variables
latest_test_date = generateLaterDate(earliest_test_date, test_date_range)

# Length of training data
train_len = "8W"
#train_len_sweep = ["4W"] #multi-option not fully implemented
# Length of testing data
test_len = "1D"
#test_len_sweep = ["1D","3D","7D"] #multi-option not fully implemented

# Time step between different experiments
#test_date_step = "1D"
# We have currently decided to step forward the experiment so that test sets
#  do not overlap, the reasoning being roughly: why would we bother evaluating
#  a model on 7 days of data if we're about to retrain the model 1 day later?
#  In the future it is possible this may change if we find a compelling reason
#  otherwise, or may add an option to override this choice.
test_date_step = test_len

# List of all experiment dates
start_test_list = generateDateRange(earliest_test_date, latest_test_date, test_date_step)
# Number of different experiment dates
total_num_exp_dates = len(start_test_list)



coverage_rate_sweep = [0.01, 0.02, 0.05, 0.10]


cell_sampling = 15 #!!! need to find where to use this, for rhs


# Knox statistic parameters
#knox_space_bin_size = 100
#knox_space_bin_count = 5
#knox_space_bins = [(i*knox_space_bin_size,(i+1)*knox_space_bin_size) \
#                   for i in range(knox_space_bin_count)]
#print(knox_space_bins)
#knox_time_bin_size = 3
#knox_time_bin_count = 7
#knox_time_bins = [(i*knox_time_bin_size,(i+1)*knox_time_bin_size) \
#                  for i in range(knox_time_bin_count)]
#print(knox_time_bins)



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





print("...declared parameters.")


date_today = datetime.date.today()
date_today_str = "{:02}{:02}{:02}".format(date_today.year%100, date_today.month, date_today.day)

# Create csv file
#csv_fname = "test190513.csv"
out_csv_fname = "results_{}_{}_{}_{}_{}.csv".format(date_today_str, dataset_name, earliest_test_date_str, test_date_range, test_date_step)
out_csv_full_path = os.path.join(datadir, out_csv_fname)



with open(out_csv_full_path, "w") as csvf:
    writer = csv.writer(csvf, delimiter=",", lineterminator="\n")
    writer.writerow(result_info_header)
    
    
    
    
    # PARAM:crime type
    for crime_type_set in crime_type_set_sweep:
        
        crime_types_printable = "_".join(sorted(crime_type_set))
        
        ### OBTAIN FULL DATA
        print("Obtaining full data set and region...")
        obtain_data_start_time = time.time()
        points_crime = chicago.load(chicago_file_path, crime_type_set, 
                                    type=chicago_load_type)
        
        
        ### OBTAIN GRIDDED REGION
        
        # Obtain polygon shapely object for region of interest
        region_polygon = chicago.get_side(chicago_side)
        
        # Obtain data set within relevant region
        points_crime_region = open_cp.geometry.intersect_timed_points(points_crime, 
                                                                      region_polygon)
        
        
        obtain_data_end_time = time.time()
        print("...obtained full data set and region.")
        print("Time: {}".format(obtain_data_end_time - obtain_data_start_time))
        
        
        
        
        # PARAM:cell width
        for cell_width in cell_width_sweep:
            
            
            print("Obtaining grid for region...")
            obtain_reg_start_time = time.time()
            
            # Obtain grid with cells only overlaid on relevant region
            masked_grid_region = open_cp.geometry.mask_grid_by_intersection(region_polygon, open_cp.data.Grid(xsize=cell_width, ysize=cell_width, xoffset=0, yoffset=0))
            
            # Get a list/tuple of all cellcoords in the region
            cellcoordlist_region = getRegionCells(masked_grid_region)
            
            # Obtain number of cells in the grid that contain relevant geometry
            # (i.e., not the full rectangular grid, only relevant cells)
            num_cells_region = len(cellcoordlist_region)
            coverage_cell_index_map = dict([(c, int(num_cells_region * c)-1) for c in coverage_rate_sweep])
            
            obtain_reg_end_time = time.time()
            print("...obtained grid for region.")
            print("Time: {}".format(obtain_reg_end_time - obtain_reg_start_time))
            
            
            
            # Log of how long an experiment takes to run
            exp_times = []
            
            
            
            # PARAM: Start date
            for exp_date_index, start_test in enumerate(start_test_list):
                
                exp_start_time = time.time()
                
                if exp_date_index % 5 == 0:
                    print("Running experiment {}/{}...".format(exp_date_index, total_num_exp_dates))
                
                # Declare time ranges of training and testing data
                end_train = start_test
                start_train = generateEarlierDate(end_train, train_len)
                end_test = generateLaterDate(start_test, test_len)
                
                #multi-option for train and test data time ranges is not
                # implemented, but would be implemented here
                #start_train_sweep = [generateEarlierDate(end_train, train_len) for train_len in train_len_sweep]
                #end_test_sweep = [generateLaterDate(start_test, test_len) for test_len in test_len_sweep]
                
                
                ### SELECT TRAINING DATA
                
                # Get subset of data for training
                points_crime_region_train = getTimedPointsInTimeRange(points_crime_region, 
                                                                      start_train, 
                                                                      end_train)
                
                
                
                ### TESTING DATA, USED FOR EVALUATION
                # (Also used for Ideal model, which is why we create it here)
                
                # Obtain selection of data for testing
                points_crime_region_test = getTimedPointsInTimeRange(points_crime_region, 
                                                                      start_test, 
                                                                      end_test)
                # Count how many crimes there were in this test data set
                num_crimes_test = len(points_crime_region_test.timestamps)
                
                # Count the number of crimes per cell
                #  This is used for evaluation and also for the "ideal" model
                cells_testcrime_ctr = countPointsPerCell(points_crime_region_test, 
                                                         masked_grid_region)
                
                
                
                
                # PARAM: model
                
                for model_name in models_to_run:
                    
                    args_to_use = []
                    if model_name in ["ideal", "random"]:
                        args_to_use.append(cellcoordlist_region)
                        if model_name == "ideal":
                            args_to_use.append(cells_testcrime_ctr)
                    elif model_name in ["naivecount", "rhs", "phs"]:
                        args_to_use.append(points_crime_region_train)
                        args_to_use.append(masked_grid_region)
                    
                    
                    # PARAM: param sweep for specific model
                    
                    for param_combo_index, param_combo in enumerate(model_param_dict[model_name]):
                        
                        sorted_cells = runModelAndSortCells(model_name, args_to_use + list(param_combo))
                        
                        print(sorted_cells[:10])
                        
                        hit_rate_list = getHitRateList(sorted_cells, cells_testcrime_ctr)
                        
                        
                        # PARAM: coverage
                        
                        for coverage_rate in coverage_rate_sweep:
                            
                            num_hits = hit_rate_list[coverage_cell_index_map[coverage_rate]]
                            pct_hits = 0
                            if num_crimes_test>0:
                                pct_hits = num_hits / num_crimes_test
                            
                            
                            
                            result_info = [
                                    dataset_name, 
                                    crime_types_printable, 
                                    cell_width, 
                                    start_test, 
                                    train_len, 
                                    test_len, 
                                    coverage_rate, 
                                    num_crimes_test, 
                                    num_hits, 
                                    pct_hits, 
                                    model_name]
                            
                            if model_name == "phs":
                                result_info += ["",""]
                            
                            result_info += list(param_combo)
                            
                            while len(result_info)<18:
                                result_info.append("")
                            
                            writer.writerow(result_info)
                        print(f"Wrote results for date {start_test} model {model_name} paramset {param_combo_index}")
                            
                
                
                exp_times.append(time.time() - exp_start_time)
            
            
            exp_times_sorted = sorted(exp_times)
            exp_times_total = sum(exp_times)
            print("Total experiment time: {}".format(exp_times_total))
            print("Average experiment time: {}".format(exp_times_total/len(exp_times)))
            print("Min experiment time: {}".format(exp_times_sorted[0]))
            print("Max experiment time: {}".format(exp_times_sorted[-1]))


print("Total time: {}".format(time.time() - init_start_time))

sys.exit(0)
