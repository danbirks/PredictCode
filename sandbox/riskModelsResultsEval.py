# -*- coding: utf-8 -*-
"""
Created on Wed May 15 10:22:47 2019

@author: lawdfo


Purpose:
    Read in the csv results file generated by (e.g.) riskModelsParamSweep.py
    and report back some useful statistics.


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
from copy import deepcopy
import statistics

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
"""
import open_cp.geometry
import open_cp.plot
import open_cp.sources.chicago as chicago
import open_cp.retrohotspot as retro
import open_cp.prohotspot as phs
import open_cp.knox
"""


# Load custom functions that make dealing with datetime and timedelta easier
from crimeRiskTimeTools import generateDateRange, \
                               generateLaterDate, \
                               generateEarlierDate, \
                               getTimedPointsInTimeRange, \
                               getSixDigitDate, \
                               _day




"""
Expected data format of input CSV file, by column:
Header name     Type            Typical contents
dataset         str             Chicago
event_types     str             BURGLARY
cell_width      int             100
eval_date       np.datetime64   2016-03-01
train_len       str             8W
test_len        str             1D
coverage_rate   float           0.01/0.02/0.05/0.1
test_events     int             3/2/5/etc
hit_count       int             1/2/0/etc
hit_pct         float           0.33333 etc
model           str             naivecount/phs/etc
rand_seed       int             
rhs_bandwidth   int             
phs_time_unit   str             1 weeks
phs_time_band   str             4 weeks
phs_dist_unit   int             100
phs_dist_band   int             400
phs_weight      str             linear

"""

csv_data_types = [str, \
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






def splitDataByTimespans(datalist, timespan, dateinfoname="eval_date"):
    print("Performing splitDataByTimespans")
    date_list = sorted(set([d[dateinfoname] for d in datalist]))
    earliest_date = date_list[0]
    latest_date = date_list[-1]
    daterange_list = generateDateRange(earliest_date, latest_date+_day, timespan)
    data_by_daterange = defaultdict(list)
    for d in datalist:
        d_time = d[dateinfoname]
        for t in daterange_list:
            if d_time >= t and d_time < generateLaterDate(t, timespan):
                data_by_daterange[t].append(d)
                break
    print("Ending splitDataByTimespans")
    return data_by_daterange


"""
Each element of output should have this info:
    earliest test date of range
    time band
    dist band
    avg hit rate
"""
def getPhsSpanStats(datalist, timespan):
    print("Performing getPhsSpanStats")
    
    data_by_daterange = splitDataByTimespans(datalist, timespan)
    
    phs_band_rate_summary = []
    
    for daterange in data_by_daterange:
        hit_rates = getPhsHitRates(data_by_daterange[daterange])
        for bp in hit_rates:
            phs_band_rate_summary.append((daterange, bp[0], bp[1], hit_rates[bp]["avg_hit_rate"]))
    
    print("Ending getPhsSpanStats")
    return phs_band_rate_summary


def getModelSpanStats(datalist, timespan, model):
    print("Performing getModelSpanStats")
    recognized_model_list = ["random", "naive", "ideal", "phs"]
    if model not in recognized_model_list:
        print("model required for getModelSpanStats")
        sys.exit(1)
    
    if model=="phs":
        return getPhsSpanStats(datalist, timespan)
    
    data_by_daterange = splitDataByTimespans(datalist, timespan)
    
    model_stats = []
    for daterange, data in data_by_daterange.items():
        model_stats.append((daterange, getAvgHitRates(data)))
    
    print("Ending getModelSpanStats")
    return model_stats




"""
Each element of output should have this info:
    coverage
    earliest test date of range
    time band
    dist band
    avg hit rate
"""
def writeModelSummaryCsv(datalists_by_cov, timespan, model, csvname = "temp.csv"):
    print("Performing writeModelSummaryCsv")
    rate_summaries_by_cov = dict()
    for cov, datalist in datalists_by_cov.items():
        rate_summaries_by_cov[cov] = getModelSpanStats(datalist, timespan, model)
    
    
    
    with open(csvname,"w") as csvf:
        writer = csv.writer(csvf, delimiter=",", lineterminator="\n")
        for cov, rate_summary in rate_summaries_by_cov.items():
            for d in rate_summary:
                writer.writerow([cov] + list(d))
    print("Ending writeModelSummaryCsv")
    sys.exit(0)




def writePhsVariabilityCsv(datalists_by_cov, timespan, csvname = "temp.csv"):
    print("Performing writePhsVariabilityCsv")
    bp_rate_summaries_by_cov = dict()
    for cov, datalist in datalists_by_cov.items():
        bp_rate_summaries_by_cov[cov] = getPhsSpanStats(datalist, timespan)
    
    rates_by_covtimedist = defaultdict(list)
    for cov, rate_summary in bp_rate_summaries_by_cov.items():
        for entry in rate_summary:
            rates_by_covtimedist[(cov, entry[1], entry[2])].append(entry[3])
    covtimedist_trios = sorted(rates_by_covtimedist)
    num_rates_list = [len(rates_by_covtimedist[x]) for x in covtimedist_trios]
    num_rates = num_rates_list[0]
    if not all([x==num_rates for x in num_rates_list]):
        print("Error! Not all (cov, time, dist) trios have same number of results!")
        print(num_rates_list)
        sys.exit(1)
    ratestats_by_covtimedist = dict()
    for covtimedist in covtimedist_trios:
        ratelist = rates_by_covtimedist[covtimedist]
        rate_avg = sum(ratelist)/num_rates
        rate_std = statistics.stdev(ratelist)
        rate_var = statistics.variance(ratelist)
        ratestats_by_covtimedist[covtimedist] = (rate_avg, rate_std, rate_var)
    
    
    with open(csvname,"w") as csvf:
        writer = csv.writer(csvf, delimiter=",", lineterminator="\n")
        for covtimedist, ratestats in ratestats_by_covtimedist.items():
            writer.writerow(list(covtimedist) + list(ratestats))
            print(" ".join([str(x) for x in list(covtimedist) + list(ratestats)]))
    
    print("Ending writePhsVariabilityCsv")
    sys.exit(0)
        





# datalist = list of results for PHS
# timespan = how frequently to check scores. Do we look at the top n models
#             from each day, or averaged over each month, etc
# topnum   = how many of the top models we consider successful. Top 10? Top 1?
def checkPhsConsistency(datalist, timespan, topnum):
    print("Performing checkPhsConsistency")
    
    data_by_daterange = splitDataByTimespans(datalist, timespan)
    
    #best_overallrate_bps = []
    best_avgrate_bps = []
    
    for daterange in data_by_daterange:
        rate_info = getPhsHitRates(data_by_daterange(daterange))
        
        #d_sort_overallrate = sorted(rate_info.items(), key=lambda ri: ri[1]["overall_hit_rate"], reverse=True)
        d_sort_avgrate = sorted(rate_info.items(), key=lambda ri: ri[1]["avg_hit_rate"], reverse=True)
        
        #best_bp_overallrate = d_sort_overallrate[:topnum]
        #for d in d_sort_overallrate[topnum:]:
        #    if d[1]["overall_hit_rate"] < best_bp_overallrate[-1][1]["overall_hit_rate"]:
        #        break
        #    best_bp_overallrate.append(d)
        
        
        best_bp_avgrate = d_sort_avgrate[:topnum]
        for d in d_sort_avgrate[topnum:]:
            if d[1]["avg_hit_rate"] < best_bp_avgrate[-1][1]["avg_hit_rate"]:
                break
            best_bp_avgrate.append(d)
        
        #best_overallrate_bps.append(best_bp_overallrate)
        best_avgrate_bps.append(best_bp_avgrate)
        
    
    #findMinimalPhsBandCovering(best_overallrate_bps, daterange_list)
    findMinimalPhsBandCovering(best_avgrate_bps, data_by_daterange.keys())






def findMinimalPhsBandCovering(best_bps, daterange_list):
    print("Performing findMinimalPhsBandCovering")
    covered_span_list = []
    covered_span_dates = [[]]
    covering_bps = []
    bp_set = set([x[0] for x in best_bps[0]])
    running_span_count = 0
    for i, bp_info in enumerate(best_bps):
        new_bp_set  = bp_set & set(x[0] for x in bp_info)
        # If set is 0, we can no longer cover current time with top choices
        if len(new_bp_set) == 0:
            covered_span_list.append(int(running_span_count+0))
            running_span_count = 1
            covered_span_dates.append([])
            covering_bps.append(deepcopy(bp_set))
            bp_set = set(x[0] for x in bp_info)
        else:
            bp_set = new_bp_set
            running_span_count += 1
        
        covered_span_dates[-1].append(daterange_list[i])
        
    covered_span_list.append(int(running_span_count+0))
    covering_bps.append(deepcopy(bp_set))
    
    print(covered_span_list)
    print(covered_span_dates)
    print(covering_bps)
    sorted_covered_span_list = sorted(covered_span_list)
    print(sorted_covered_span_list[0])
    print(sorted_covered_span_list[-1])
    
    #sys.exit(0)
    
    pass



"""
getPhsHitRates
Input: "datalist" = list where each entry is a dictionary containing the
                     information from a line of the csv results file (casted
                     as the appropriate data type) as well as "param_pair"
                     which is a tuple of the time and dist bandwidths.
       Note: Ideally this datalist is a subset of the full csv data, so that
              hit rates ar calculated over smaller timespans, e.g. monthly
Output: "info_by_band_pair" = dict that maps bandwidth pairs ("bp") to:
            "bands": same as key; can be useful if just grabbing values
            "num_tests": Number of experiments/tests/evaluations performed.
                All bp's within a datalist fed into this function should end
                up with the same number of tests -- I can't think of a reason
                why this wouldn't happen. However, note that this number MAY
                change across multiple runs of this function with different
                data subsets. For example, maybe you calculate over every
                month, but months have different numbers of days.
            "total_events": Total number of events (i.e. crimes) in the data.
                This is calculated by adding the number for the first time
                each date is witnessed. So again, it's important that all bp's
                are tested on all the same days.
            "total_hits": Total number of hits achieved by the bp's model.
            "total_rates": Sum of all daily(?) hit rates. This number is
                essentially useless on its own, but used for calculating avg.
            "avg_hit_rate": Average of all daily hit rates, calculated as
                total_rates/num_tests
            ("overall_hit_rate"): A different average hit rate, being the total
                number of hits divided by the total number of events. This
                was removed from use (commented out) once we decided this
                metric was less useful than avg_hit_rate, since this could be
                swayed by a generally poor model that rarely performs extremely
                well.
"""
def getPhsHitRates(datalist):
    print("Performing getPhsHitRates")
    
    # Obtain set of bandwidths
    band_pair_list = sorted(set([d["param_pair"] for d in datalist]))
    
    # Instantiate info to obtain
    info_by_band_pair = dict()
    for bp in band_pair_list:
        info_by_band_pair[bp] = dict([\
                         ("bands", bp),\
                         ("num_tests", 0),\
                         ("total_events", 0),\
                         ("total_hits", 0),\
                         ("total_rates", float(0))\
                         ])
    
    # Update info via running counts for each bandwidth pair
    for result in datalist:
        bp = result["param_pair"]
        info_by_band_pair[bp]["num_tests"] += 1
        info_by_band_pair[bp]["total_events"] += result["test_events"]
        info_by_band_pair[bp]["total_hits"] += result["hit_count"]
        if result["test_events"] > 0:
            info_by_band_pair[bp]["total_rates"] += result["hit_count"]/result["test_events"]
    
    # Confirm all bandwidth pairs had the same number of tests
    num_tests_per_bp = [info_by_band_pair[bp]["num_tests"] for bp in band_pair_list]
    if len(set(num_tests_per_bp)) != 1:
        print("Error! Some bandwidth pairs have different numbers of tests!")
        print(Counter(num_tests_per_bp))
        sys.exit(1)
    num_tests = num_tests_per_bp[0]
    
    # Compute the average hit rates for each bandwidth pair
    for bp in band_pair_list:
        info_by_band_pair[bp]["avg_hit_rate"] = info_by_band_pair[bp]["total_rates"]/num_tests
        
        # The following deprecated code computes the overall hit rate,
        #  instead of averaging the hit rates
        #if info_by_band_pair[bp]["total_events"] == 0:
        #    info_by_band_pair[bp]["overall_hit_rate"] == 0
        #else:
        #    info_by_band_pair[bp]["overall_hit_rate"] = info_by_band_pair[bp]["total_hits"]/info_by_band_pair[bp]["total_events"]
    
    # Return info
    return info_by_band_pair
    




# Note: 0 hits for 0 events gets counted as a hit rate of 0.
#       Perhaps it should be discarded instead?
#       But then what if the entire span has 0 events?
def getAvgHitRates(datalist):
    print("Performing getAvgHitRates")
    num_tests = len(datalist)
    total_rates = sum([result["hit_count"]/result["test_events"] for result in datalist if result["test_events"]!=0])
    for i, result in enumerate(datalist):
        print(result)
        toprint = [result["hit_count"], result["test_events"]]
        if toprint[1] == 0:
            toprint.append(0)
        else:
            toprint.append(toprint[0]/toprint[1])
        print("\t".join([str(x) for x in toprint]))
    print(total_rates/num_tests)
    return total_rates/num_tests




"""
getDataByCovRate

Given a path to csv results from running risk models,
 return a dictionary where keys are coverage rates and 
 values are the rows of info with that coverage from the csv.
"""
def getDataByCovRate(results_full_path, header_types = csv_data_types):
    
    # Keep track of total number of events (i.e., crimes)
    total_event_count = 0
    
    # Instantiate a mapping from coverage rate to {another mapping of results}.
    # That other mapping will be from model to results.
    # And, those results will be a list of mappings, each entry in the list being
    #  a different row from the csv results
    datadicts_by_cov_rate = defaultdict(lambda: defaultdict(list))
    
    # Open csv output and start reading it
    with open(results_full_path, newline="") as f:
        reader = csv.reader(f)
        
        # Obtain column names from header in first line
        header = next(reader, None)
        
        # Read each line of data
        for dataline in reader:
            # Instantiate a map from col name to data, for this line
            dataline_dict = dict()
            
            # All data is currently in string form.
            # Use header_types to cast the data appropriately.
            for i,d in enumerate(dataline):
                # Default is empty string
                casted_data = ""
                # Things like int("") don't work, so we catch that here
                if d != "":
                    casted_data = header_types[i](d)
                # Transform data into str/int/float/datetime64 before storing it
                dataline_dict[header[i]] = casted_data
            
            
            # Keep track of how many eval_date's we've seen,
            #  and how many events (crimes) there have been in total
            # If date is outside of desired range, continue
            dataline_date = dataline_dict["eval_date"]
            if earliest_eval_date != None and dataline_date < earliest_eval_date:
                continue
            if latest_eval_date != None and latest_eval_date < dataline_date:
                continue
            if dataline_date not in dates_seen:
                total_event_count += dataline_dict["test_events"]
                dates_seen.add(dataline_date)
            
            # Grab coverage and model, since we'll use those a lot
            dataline_cov = dataline_dict["coverage_rate"]
            dataline_model = dataline_dict["model"]
            
            # Grab the bandwidths for PHS results, store them as "param_pair"
            if dataline_model == "phs":
                time_band = int(dataline_dict["phs_time_band"].split()[0])
                dist_band = dataline_dict["phs_dist_band"]
                dataline_dict["param_pair"] = (time_band, dist_band)
            
            # Store dict so they're first sorted by coverage then by model type
            datadicts_by_cov_rate[dataline_cov][dataline_model].append(dataline_dict)
            
    return datadicts_by_cov_rate







datadir = os.path.join("..", "..", "Data")
#results_fname = "results_190515_Chicago_160101_1M_1D.csv"
#results_fname = "results_190517_Chicago_020101_1Y_1D.csv"
#results_fname = "results_190517_Chicago_020101_1Y_3D.csv"
#results_fname = "results_190522_Chicago_020101_1Y_7D.csv"
#results_fname = "results_190621_Chicago_160301_1M_1D.csv"



#results_fname = "results_190621_Chicago_160301_9M_1D.csv"
#results_fname = "temp_results_190621_Chicago_010301_17Y_1D.csv"
#results_fname = "results_190621_Chicago_010301_17Y_1D.csv"
results_fname = "results_190628_Chicago_130101_5Y_1D.csv"


# Only include results of tests later OR EQUAL to this date
earliest_eval_date = np.datetime64("2013-01-01")
# Only include results of tests earlier BUT NOT EQUAL to this date
latest_eval_date = None





results_full_path = os.path.join(datadir, results_fname)


# Keep track of dates seen in the output data
dates_seen = set()

datadicts_by_cov_rate = getDataByCovRate(results_full_path)



# Determine the number of evaluation dates in the data
# We expect this to equal the number of instances of random/naive/ideal
#  experiments, and also the number of phs experiments when multiplied by
#  the number of phs parameter combinations.
num_dates = len(dates_seen)
print(num_dates)
earliest_date_seen =sorted(dates_seen)[0]
latest_date_seen =sorted(dates_seen)[-1]
print(earliest_date_seen)
print(latest_date_seen)

phsdicts_by_cov_rate = dict([(cov, d["phs"]) for cov, d in datadicts_by_cov_rate.items()])

naivedicts_by_cov_rate = dict([(cov, d["naive"]) for cov, d in datadicts_by_cov_rate.items()])


create_naive_csv_summary = True
if create_naive_csv_summary:
    
    timespan = "1M"
    date_today = datetime.date.today()
    date_today_str = getSixDigitDate(date_today)
    earliest_date_str = getSixDigitDate(earliest_date_seen)
    latest_date_str = getSixDigitDate(latest_date_seen)
    sumcsv_base = f"ratesummary_xsr_nai_{date_today_str}_{earliest_date_str}_{latest_date_str}_{timespan}.csv"
    sumcsvname = os.path.join(datadir, sumcsv_base)
    writeModelSummaryCsv(naivedicts_by_cov_rate, timespan, "naive", csvname=sumcsvname)
    
    
    
    sys.exit(0)

create_phs_csv_summary = False
if create_phs_csv_summary:
    
    
    
    timespan = "1M"
    date_today = datetime.date.today()
    date_today_str = getSixDigitDate(date_today)
    earliest_date_str = getSixDigitDate(earliest_date_seen)
    latest_date_str = getSixDigitDate(latest_date_seen)
    phssumcsv_base = f"ratesummary_xsr_phs_{date_today_str}_{earliest_date_str}_{latest_date_str}_{timespan}.csv"
    phssumcsvname = os.path.join(datadir, phssumcsv_base)
    #writePhsSummaryCsv(phs_list, timespan, csvname=phssumcsvname)
    writeModelSummaryCsv(phsdicts_by_cov_rate, timespan, "phs", csvname=phssumcsvname)
    
    
    
    sys.exit(0)


create_phs_csv_var = True
if create_phs_csv_var:
    
    
    
    timespan = "1M"
    date_today = datetime.date.today()
    date_today_str = getSixDigitDate(date_today)
    earliest_date_str = getSixDigitDate(earliest_date_seen)
    latest_date_str = getSixDigitDate(latest_date_seen)
    phssumcsv_base = f"ratevar_{date_today_str}_{earliest_date_str}_{latest_date_str}_{timespan}.csv"
    phssumcsvname = os.path.join(datadir, phssumcsv_base)
    #writePhsSummaryCsv(phs_list, timespan, csvname=phssumcsvname)
    writePhsVariabilityCsv(phsdicts_by_cov_rate, timespan, phssumcsvname)
    
    
    
    sys.exit(0)











all_model_names = ["random", "naivecount", "ideal", "rhs", "phs"]
basic_model_names = all_model_names[:3]

for cov, datadicts_by_model in datadicts_by_cov_rate.items():
    print(f"Coverage rate: {cov}")
    
    # Get overall result summaries for basic models
    for model_name in basic_model_names:
        if model_name in datadicts_by_model:
            
            # Obtain list of results for this (coverage, model) combo
            datalist = datadicts_by_model[model_name]
            
            # Confirm that we have the expected number of results
            if len(datalist) != num_dates:
                print("Error! Unexpected number of results!")
                print(f"Number expected per model: {num_dates}")
                print(f"Number seen for model {model_name}: {len(datalist)}")
                sys.exit(1)
            
            #  ("Hit" = event in testing period within model's top cov% cells)
            # Total number of successful "hits"
            total_hit_count = sum([d["hit_count"] for d in datalist])
            # Total possible number of "hits"
            total_hit_poss = sum([d["test_events"] for d in datalist])
            # Overall hit rate
            total_hit_rate = total_hit_count/total_hit_poss
            # Average of all individual hit rates
            average_hit_rate = sum(d["hit_pct"] for d in datalist)/num_dates
            print(f"\tModel: {model_name}")
            #print(f"\t\tTotal hit rate: {total_hit_count}/{total_hit_poss} = {total_hit_rate:6.4f}")
            print(f"\t\tAverage hit rate: {average_hit_rate:6.4f}")
    
    # Generate a table of results for all PHS bandwidth pairs tested
    if "phs" in datadicts_by_model:
        phs_list = datadicts_by_model["phs"]
        
        # phs_list is what should be fed into checkPhsConsistency etc
        
        #checkPhsConsistency(phs_list, "1M", 10)
        #sys.exit(0)
        #continue
        
        #getPhsStats(phs_list, "1M")
        #sys.exit(0)
        
        
        
        print("0\t" + "\t".join([str(x) + " weeks" for x in range(1,9)]))
        best_sum_dist_time = (-1, 0, 0)
        best_avg_dist_time = (-1, 0, 0)
        for dist_band in range(100,1100,100):
            toprint_list = [str(dist_band)]
            for time_band in range(1,9):
                #total_hit_count = sum(d["hit_count"] for d in phs_list if d["param_pair"]==(time_band, dist_band))
                #total_hit_rate = total_hit_count/total_event_count
                hit_rate_sum = sum(d["hit_pct"] for d in phs_list if d["param_pair"]==(time_band, dist_band))
                hit_rate_avg = hit_rate_sum/num_dates
                #toprint_list.append(f"{total_hit_rate:6.4f},{hit_rate_avg:6.4f}")
                toprint_list.append(f"{hit_rate_avg:6.4f}")
                #if total_hit_rate > best_sum_dist_time[0]:
                #    best_sum_dist_time = (total_hit_rate, dist_band, time_band)
                if hit_rate_avg > best_avg_dist_time[0]:
                    best_avg_dist_time = (hit_rate_avg, dist_band, time_band)
            print("\t".join(toprint_list))
        #print(f"Best total hit rate result: {best_sum_dist_time[0]:6.4f} {best_sum_dist_time[1:]}")
        print(f"Best average hit rate result: {best_avg_dist_time[0]:6.4f} {best_avg_dist_time[1:]}")
    
    


sys.exit(0)



if __name__ == "__main__":
    main()
