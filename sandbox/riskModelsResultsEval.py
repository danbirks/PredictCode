# -*- coding: utf-8 -*-
"""
Created on Wed May 15 10:22:47 2019

@author: lawdfo
"""


# Some fairly standard modules
import os, csv, lzma
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import descartes
from itertools import product
from collections import Counter, defaultdict
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

"""
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

datadir = os.path.join("..", "..", "Data")
#results_fname = "results_190515_Chicago_160101_1M_1D.csv"
#results_fname = "results_190517_Chicago_020101_1Y_1D.csv"
results_fname = "results_190517_Chicago_020101_1Y_3D.csv"
#results_fname = "results_190522_Chicago_020101_1Y_7D.csv"


results_full_path = os.path.join(datadir, results_fname)

rel_results = []
rel_res_limit = 999999999



header_types = [str, \
                str, \
                int, \
                np.datetime64, \
                str, \
                str, \
                float, \
                int, \
                int, \
                float, \
                str, \
                int, \
                int, \
                str, \
                str, \
                int, \
                int, \
                str]

dates_seen = set()
total_event_count = 0

#datadicts = []
datadicts_by_cov_rate = defaultdict(lambda: defaultdict(list))

with open(results_full_path, newline="") as f:
    reader = csv.reader(f)
    
    # Obtain column names from header in first line
    header = next(reader, None)
    # Now the header is essentially a map from column number to column name.
    # Let's also get the inverse, mapping column name to column number
    #header_map = dict([(h,i) for i,h in enumerate(header)])
    
    # Read each line of data
    for dataline in reader:
        # Create a map from col name to data for this line
        dataline_dict = dict()
        for i,d in enumerate(dataline):
            # Default is empty string
            casted_data = ""
            # Things like int("") don't work, so we catch that here
            if d != "":
                casted_data = header_types[i](d)
            # Transform data into str/int/float/datetime64 before storing it
            dataline_dict[header[i]] = casted_data
        dataline_cov = dataline_dict["coverage_rate"]
        dataline_model = dataline_dict["model"]
        
        if dataline_model == "phs":
            time_band = int(dataline_dict["phs_time_band"].split()[0])
            dist_band = dataline_dict["phs_dist_band"]
            dataline_dict["param_pair"] = (time_band, dist_band)
        
        # Store dict so they're first sorted by coverage then by model type
        datadicts_by_cov_rate[dataline_cov][dataline_model].append(dataline_dict)
        # Keep track of how many eval_date's we've seen
        dataline_date = dataline_dict["eval_date"]
        if dataline_date not in dates_seen:
            total_event_count += dataline_dict["test_events"]
            dates_seen.add(dataline_date)
        

# Determine the number of evaluation dates in the data
# We expect this to equal the number of instances of random/naive/ideal
#  experiments, and also the number of phs experiments when multiplied by
#  the number of phs parameter combinations.
num_dates = len(dates_seen)


for cov, datadicts_by_model in datadicts_by_cov_rate.items():
    print(f"Coverage rate: {cov}")
    
    basic_model_names = ["random","naivecount","ideal"]
    
    for model_name in basic_model_names:
        datalist = datadicts_by_model[model_name]
        if len(datalist) != num_dates:
            print("Error! Unexpected number of results!")
            print(f"Number expected per model: {num_dates}")
            print(f"Number seen for model {model_name}: {len(datalist)}")
            sys.exit(1)
        total_hit_count = sum([d["hit_count"] for d in datalist])
        total_hit_poss = sum([d["test_events"] for d in datalist])
        average_hit_rate = sum(d["hit_pct"] for d in datalist)/num_dates
        print(f"\tModel: {model_name}")
        print(f"\t\tTotal hit rate: {total_hit_count}/{total_hit_poss} = {total_hit_count/total_hit_poss:6.4f}")
        print(f"\t\tAverage hit rate: {average_hit_rate:6.4f}")
    
    phs_list = datadicts_by_model["phs"]
    print("0\t" + "\t\t".join([str(x) + " weeks" for x in range(1,7)]))
    for dist_band in range(100,1100,100):
        toprint_list = [str(dist_band)]
        for time_band in range(1,7):
            total_hit_count = sum(d["hit_count"] for d in phs_list if d["param_pair"]==(time_band, dist_band))
            total_hit_rate = total_hit_count/total_event_count
            hit_rate_sum = sum(d["hit_pct"] for d in phs_list if d["param_pair"]==(time_band, dist_band))
            hit_rate_avg = hit_rate_sum/num_dates
            toprint_list.append(f"{total_hit_rate:6.4f},{hit_rate_avg:6.4f}")
        print("\t".join(toprint_list))
    
    



sys.exit(0)

# What inputs? coverage, time band, dist band
# What outputs? hit count, hit pct, maybe hit poss?

num_results = len(rel_results)
summaries_by_cov = dict()
cov_rates = [0.01,0.02,0.05,0.1]
for cr in cov_rates:
    summaries_by_cov[cr] = defaultdict(list)

for i,r in enumerate(rel_results):
    time_band = int(r["phs_time_band"].split()[0])
    dist_band = int(r["phs_dist_band"])
    cov_rate = float(r["coverage_rate"])
    hit_count = int(r["hit_count"])
    hit_rate = float(r["hit_pct"])
    event_count = int(r["test_events"])
    #if (time_band, dist_band) == (2,100):
    #    print(i)
    #print(f"({time_band}, {dist_band})")
    summaries_by_cov[cov_rate][(time_band,dist_band)].append((hit_count, event_count, hit_rate))
    #print(f'Exp {i} of model {r["model"]} has coverage {r["coverage_rate"]}, \
    #hit count {r["hit_count"]} and success {r["hit_pct"]}, \
    #with bands {r["phs_time_band"]} and {r["phs_dist_band"]}')
    #print("")

for cov_rate in sorted(summaries_by_cov):
    print(cov_rate)
    param_results = summaries_by_cov[cov_rate]
    param_totals = dict()
    for param_pair in param_results:
        totals = [sum([x[vi] for x in param_results[param_pair]]) for vi in range(3)]
        param_totals[param_pair] = totals
    
    print("\t".join([str(x) for x in range(7)]))
    for dist in range(100,1100,100):
        print(str(dist) + "\t" + "\t".join([str(param_totals[(wks, dist)][0]) for wks in range(1,7)]))
    print([len(param_results[pp]) for pp in param_results])
    
    
    



print("Done!")
sys.exit(0)