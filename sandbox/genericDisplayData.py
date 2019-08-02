# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:34:27 2019

@author: lawdfo
"""

import sys
import os
import datetime
import csv
import numpy as np
import time


import geopandas as gpd

# Elements from PredictCode's custom "open_cp" package
sys.path.insert(0, os.path.abspath(".."))
#import open_cp
#import open_cp.geometry
#import open_cp.sources.chicago as chicago
from open_cp.data import TimedPoints





# Load custom functions that make dealing with datetime and timedelta easier
from crimeRiskTimeTools import generateLaterDate, \
                               generateDateRange, \
                               getSixDigitDate





def onetimeConfirmSameData():
    
    # This code just confirmed that the new method gets same data as the old
    chicago_file_name = "chi_all_s_BURGLARY_RES_010101_190101.csv"
    chicago_file_path = os.path.join(datadir, chicago_file_name)
    import open_cp.sources.chicago as chicago
    points_crime_chiorig = chicago.load(chicago_file_path, crime_type_set, 
                                type="all")
    pc_t = points_crime.timestamps
    pc_co_t = points_crime_chiorig.timestamps
    print(len(pc_t))
    print(len(pc_co_t))
    print(all(pc_t==pc_co_t))
    
    pc_x = points_crime.xcoords
    pc_y = points_crime.ycoords
    pc_co_x = points_crime_chiorig.xcoords
    pc_co_y = points_crime_chiorig.ycoords
    x_diff = []
    y_diff = []
    for i in range(len(pc_t)):
        x_diff.append(abs(pc_x[i]-pc_co_x[i]))
        y_diff.append(abs(pc_y[i]-pc_co_y[i]))
    x_diff = sorted(x_diff)
    y_diff = sorted(y_diff)
    print(x_diff[:10])
    print(x_diff[-10:])
    print(sum(x_diff)/len(pc_t))
    print(y_diff[:10])
    print(y_diff[-10:])
    print(sum(y_diff)/len(pc_t))
    

def onetimeGeojsonManipulation():
    
    ###
    # Here's how a Chicago polygon was acquired:
    
    # Obtain polygon shapely object for South side
    #ss_polygon = chicago.get_side("South")
    
    
    
    geojson_dur_file_name = "durham.geojson"
    geojson_dur_full_path = os.path.join(datadir, geojson_dur_file_name)
    geo_frame_dur = gpd.read_file(geojson_dur_full_path)
    
    
    geo_frame_chi = gpd.read_file(geojson_full_path)
    chicago_side_mapping = {
        "Far North" : [1,2,3,4,9,10,11,12,13,14,76,77],
        "Northwest" : [15,16,17,18,19,20],
        "North" : [5,6,7,21,22],
        "West" : list(range(23, 32)),
        "Central" : [8,32,33],
        "South" : list(range(34,44)) + [60, 69],
        "Southwest" : [56,57,58,59] + list(range(61,69)),
        "Far Southwest" : list(range(70,76)),
        "Far Southeast" : list(range(44,56))
    }
    geo_frame_chi["side"] = geo_frame_chi.area_numbe.map(lambda x : next(key
         for key, item in chicago_side_mapping.items() if int(x) in item) )
    _sides = geo_frame_chi.drop(["area", "area_num_1", "comarea", "comarea_id",
                        "perimeter", "shape_area", "shape_len"], axis=1)
    _sides.crs = {"init": "epsg:4326"} # global coordinate system
    _sides = _sides.to_crs({"init": "epsg:2790"}) # East Illinois-specifc
    
    ss_polygon = _sides[_sides.side == "South"].unary_union
    type(ss_polygon)
    
    
    
    
    # Obtain GeoDataFrame with polygon's geometry
    #  and with CRS epsg:2790
    ss_frame = gpd.GeoDataFrame({"name":["South Side"]})
    ss_frame.geometry = [ss_polygon]
    ss_frame.crs = {"init":"epsg:2790"}
    
    
    geojson_sschi_file_name = "Chicago_South_Side_2790.geojson"
    geojson_sschi_full_path = os.path.join(datadir, geojson_sschi_file_name)
    ss_frame.to_file(geojson_sschi_full_path, driver="GeoJSON")
    
    
    
    
    # Open GeoJSON file
    geojson_dur_file_name = "durham.geojson"
    geojson_dur_full_path = os.path.join(datadir, geojson_dur_file_name)
    geo_frame_dur = gpd.read_file(geojson_dur_full_path)
    
    # Remove unnecessary "Description" column
    geo_frame_dur_nodesc = geo_frame_dur.drop(["Description"], axis=1)
    # Set the name (we might not really need to do this)
    geo_frame_dur_nodesc.iloc[0]["Name"] = "Durham"
    # Transform data to UK-specific EPSG of 27700
    geo_frame_dur_nodesc = geo_frame_dur_nodesc.to_crs({"init": "epsg:27700"})
    
    # Save off new version of geodataframe into GeoJSON file
    geojson_dur_27_file_name = "Durham_27700.geojson"
    geojson_dur_27_full_path = os.path.join(datadir, geojson_dur_27_file_name)
    geo_frame_dur_nodesc.to_file(geojson_dur_27_full_path, driver="GeoJSON")
    
    
    
    """
    Given original file:
        Use ogr2ogr to convert to geojson if necessary
            For example, given durham.kml, use this command:
            > ogr2ogr -f GeoJSON durham.geojson durham.kml
        Use some code from above to read that geojson,
            form the union of certain subsets, if necessary
            convert from EPSG 4326 to a local EPSG (like 2790 or 27700)
            save it to a new geojson
        
    
    
    """












def loadGenericData(filepath, crime_type_set = {"BURGLARY"}, date_format_csv = "%m/%d/%Y %I:%M:%S %p", epsg = None, proj=None, longlat=True, infeet=False):
    
    # Note: Data is expected in a csv file with the following properties:
    # Row 0 = Header
    # Col 0 = Date/time
    # Col 1 = Longitude
    # Col 2 = Latitude
    # Col 3 = Crime type
    # Col 4 = Location type (optional)
    
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






#### Data parameters
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

# Dataset name (to be included in name of output file)
dataset_name= "chicago"
# Crime types
#!!!  May want to change this to be different for train vs test?
crime_type_set = {"BURGLARY"}
crime_types_printable = "_".join(sorted(crime_type_set))
#crime_type_set_sweep = [{"BURGLARY"}] # if doing a sweep over different sets
# Size of grid cells
cell_width_sweep = [100]
# Input csv file name
in_csv_file_name = "chi_all_s_BURGLARY_RES_010101_190101_stdXY.csv"
# Format of date in csv file
date_format_csv = "%m/%d/%Y %I:%M:%S %p"
#!!! Should actually force that format when standardizing csv files!!!
#!!! Also change "infeet" issues at standardizing stage too!!!
# Of all planned experiments, earliest start of a test data set
earliest_test_date = "2003-01-01"
# Time between earliest experiment and latest experiment
test_date_range = "2Y"
# Length of training data
train_len = "8W"
# Length of testing data
test_len = "7D"
# Time step between different experiments
#test_date_step = "1D"
test_date_step = test_len
# Coverage rates to test
coverage_rate_sweep = [0.01, 0.02, 0.05, 0.10]
# Geojson file
#geojson_file_name = "Chicago_Areas.geojson"
geojson_file_name = "Chicago_South_Side_2790.geojson"
geojson_file_name = "Durham_27700.geojson"

# There are more parameters to include regarding the models to run
#  but for now let's just try to display the data itself
# Predictive models to run
# models_to_run = ["naivecount","phs"]




###
# Derived parameters

# Full path for input csv file
in_csv_full_path = os.path.join(datadir, in_csv_file_name)
# Nicely-formatted string of test date
earliest_test_date_str = "".join(earliest_test_date.split("-"))[2:]
# Latest start of a test data set, calculated from earliest and length
latest_test_date = generateLaterDate(earliest_test_date, test_date_range)
# List of all experiment dates
start_test_list = generateDateRange(earliest_test_date, latest_test_date, test_date_step)
# Number of different experiment dates
total_num_exp_dates = len(start_test_list)
# Output csv file name
out_csv_file_name = "results_{}_{}_{}_{}_{}.csv".format(date_today_str, dataset_name, earliest_test_date_str, test_date_range, test_date_step)
# Full path for output csv file
out_csv_full_path = os.path.join(datadir, out_csv_file_name)
# Full path for geojson file
geojson_full_path = os.path.join(datadir, geojson_file_name)

###
## Not sure if we need any of this
#chicago_side = "South"
#chicago_load_type = "snapshot"
#if "all" in chicago_file_name:
#    chicago_load_type = "all"
#chicago_file_path = os.path.join(datadir, chicago_file_name)
## Chicago module requires this line to access some data
#chicago.set_data_directory(datadir)


print("...declared parameters.")


# Open output csv file for writing, write header row
with open(out_csv_full_path, "w") as csvf:
    writer = csv.writer(csvf, delimiter=",", lineterminator="\n")
    writer.writerow(result_info_header)
    
    
    # If we were to do a "crime type sweep", that would go here.
    # But probably simpler to just re-run with a new crime type set instead.
    
    
    
    ### OBTAIN FULL DATA
    print("Obtaining full data set and region...")
    chktime_obtain_data = time.time()
    points_crime = loadGenericData(in_csv_full_path, 
                                   crime_type_set=crime_type_set, 
                                   longlat=False, 
                                   infeet=True)
    
    # Obtain polygon from geojson file (which should have been pre-processed)
    area_polygon = gpd.read_file(geojson_full_path).unary_union
    
    
    
    
    
    
    
    
    
    
    
    
    


