# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 13:35:29 2019

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
import random
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
import open_cp.geometry
import open_cp.plot
import open_cp.sources.chicago as chicago
import open_cp.retrohotspot as retro
import open_cp.prohotspot as phs
import open_cp.knox





_cdict = {'red':   [(0.0,  1.0, 1.0),
                   (1.0,  1.0, 1.0)],
         'green': [(0.0,  1.0, 1.0),
                   (1.0,  0.0, 0.0)],
         'blue':  [(0.0,  0.2, 0.2),
                   (1.0,  0.2, 0.2)]}

yellow_to_red = matplotlib.colors.LinearSegmentedColormap("yellow_to_red", _cdict)





class KnoxEntry:
    
    def __init__(self, start_date = None, end_date = None, num_events = -1):
        self.start_date = start_date
        self.end_date = end_date
        self.num_events = num_events
        self.pvals = []
        self.stats = []
        self.dists = []
        self.medians = []
        self.ratios = []
        





def knox_ratio(knox_statistic, distribution):
    """As in the paper, compute the ratio of the statistic to the median
    of the values in the distribution"""
    #d = np.array(distribution)
    #d.sort()
    #return statistic / d[len(d)//2]
    return knox_statistic / statistics.median(distribution)




def significant_cells(pvalue_array, sig_thresh=0.05):
    return pvalue_array < sig_thresh

def contiguous_cells(data_array, origin=(0,0)):
    array_dims = np.shape(data_array)
    need_to_visit_stack = [origin]
    visited_array = np.zeros_like(data_array, dtype=bool)
    contig_array = np.zeros_like(data_array, dtype=bool)
    
    val = data_array[origin]
    
    tempctr = 0
    while len(need_to_visit_stack)>0:
        curr_cell = need_to_visit_stack.pop(0)
        tempctr += 1
        if visited_array[curr_cell]:
            continue
        visited_array[curr_cell] = True
        if data_array[curr_cell] != val:
            continue
        contig_array[curr_cell] = True
        for dim_index, dim_size in enumerate(array_dims):
            cc_index_val = curr_cell[dim_index]
            if cc_index_val>0:
                need_to_visit_stack.append(curr_cell[:dim_index] + (cc_index_val-1,) + curr_cell[dim_index+1:])
            if cc_index_val<dim_size-1:
                need_to_visit_stack.append(curr_cell[:dim_index] + (cc_index_val+1,) + curr_cell[dim_index+1:])
        
    #print(data_array)
    #print(visited_array)
    #print(contig_array)
    return contig_array


def get_bandwidths_from_knox(pvalue_array, selection="contig_to_axis", sig_thresh=0.05):
    if selection not in ["contig_to_axis","contig_anywhere","along_axis"]:
        print("Error, unrecognized selction type: {}".format(selection))
        sys.exit(1)
    signif_array = significant_cells(pvalue_array, sig_thresh=sig_thresh)
    array_dims = np.shape(signif_array)
    #print("pvals")
    #print(pvalue_array)
    #print("signif")
    #print(signif_array)
    if signif_array[0,0] == False:
        print("Warning: Knox statistic of smallest time/space bin is not significant!")
        return (-1,-1)
    if selection == "along_axis":
        row_ind = 0
        while row_ind<array_dims[0]-1 and signif_array[row_ind+1, 0] == True:
            row_ind += 1
        col_ind = 0
        while col_ind<array_dims[1]-1 and signif_array[0, col_ind+1] == True:
            col_ind += 1
        return (row_ind, col_ind)
    contig_signif = contiguous_cells(signif_array, origin=(0,0))
    if selection == "contig_anywhere":
        row_ind = array_dims[0]-1
        while not any(contig_signif[row_ind,:]):
            row_ind -= 1
        col_ind = array_dims[1]-1
        while not any(contig_signif[:,col_ind]):
            col_ind -= 1
        return (row_ind, col_ind)
    if selection == "contig_to_axis":
        row_ind = array_dims[0]-1
        while not contig_signif[row_ind,0]:
            row_ind -= 1
        col_ind = array_dims[1]-1
        while not contig_signif[0,col_ind]:
            col_ind -= 1
        return (row_ind, col_ind)



def plot_signif_knox_ratios(knox_entry: KnoxEntry, sbin_size, tbin_size, p_thresh=0.05):
    
    mappable = plt.cm.ScalarMappable(cmap=yellow_to_red)
    mappable.set_array(np.ravel(knox_entry.ratios))
    mappable.autoscale()
    fig, ax = plt.subplots(figsize=(12,4))
    
    array_dims = np.shape(knox_entry.pvals)
    xmin = 0
    xmax = array_dims[0] * sbin_size
    ax.set(xlim=(xmin, xmax), xlabel="Distance in metres")
    
    ymin = 0
    ymax = array_dims[1] * tbin_size
    ax.set(ylim=(ymin, ymax), ylabel="Time in days")
    
    ax.set_title("Knox, {} events from {} to {}, p={}".format(knox_entry.num_events, knox_entry.start_date, knox_entry.end_date, p_thresh))
    
    for (sbin_index,tbin_index), pval in np.ndenumerate(knox_entry.pvals):
        if pval >= p_thresh:
            continue
        p = matplotlib.patches.Rectangle((sbin_index*sbin_size, tbin_index*tbin_size), sbin_size, tbin_size, fc=mappable.to_rgba(knox_entry.ratios[sbin_index,tbin_index]))
        ax.add_patch(p)
    
    cbar = fig.colorbar(mappable, orientation="vertical")
    cbar.set_label("Knox ratio")
    
    




datadir = os.path.join("..", "..", "Data")
chicago_file_name = "chicago_all_old.csv"
chicago_side = "South"
#crime_type_set = {"THEFT"}
crime_type_set = {"BURGLARY"}
chicago_file_path = os.path.join(datadir, chicago_file_name)
chicago.set_data_directory(datadir)




exp_limit = 100
num_time_bins = 10
num_dist_bins = 20
num_knox_iters = 100
sbin_size = 100
tbin_size = 7

knox_file_path = os.path.join(datadir, 
                              "knox_190330_sschi_burg_cell250_sbin100_tbin7_iter100.txt")
with open(knox_file_path) as kf:
    exp_num = -1
    stype = "lastline"
    sctr = 0
    sdata = []
    knox_data = []
    for lnum, kline in enumerate(kf):
        kline = kline.strip()
        
        #if stype != "distribution":
        #    print("{}:{}:{}".format(exp_num, lnum, kline))
        #if lnum >50:
        #    break
        
        
        if len(kline)==0:
            if stype == "distribution":
                knox_data[-1].dists = np.asarray(sdata).reshape(num_dist_bins, num_time_bins, num_knox_iters)
                knox_data[-1].medians = np.zeros((num_dist_bins, num_time_bins))
                knox_data[-1].ratios = np.zeros((num_dist_bins, num_time_bins))
                for i in range(num_dist_bins):
                    for j in range(num_time_bins):
                        knox_data[-1].medians[i][j] = statistics.median(knox_data[-1].dists[i][j])
                knox_data[-1].ratios = knox_data[-1].stats/knox_data[-1].medians
                
                stype = "lastline"
                continue
            elif stype != "lastline":
                print("Error, unexpected empty line at {}:{}".format(exp_num, lnum))
                sys.exit(1)
            exp_num += 1
            stype = "info"
            sctr = 0
            sdata = []
            if exp_num >= exp_limit:
                break
            continue
        
        if stype == "info":
            if kline == "p values":
                knox_data.append(KnoxEntry(*sdata))
                sdata = []
                sctr = 0
                stype = "p values"
                continue
            if sctr in [0,1]:
                sdata.append(np.datetime64(kline))
            elif sctr == 2:
                sdata.append(int(kline))
            else:
                print("Error, info section incorrect in part {}?".format(exp_num))
                sys.exit(1)
            sctr += 1
            continue
        elif stype == "p values":
            if kline == "statistics":
                knox_data[-1].pvals = sdata
                sdata = []
                sctr = 0
                stype = "statistics"
                continue
            next_row = np.array([float(x) for x in kline.split()])
            if sctr == 0:
                sdata = next_row
            else:
                sdata = np.vstack([sdata, next_row])
            sctr += 1
        elif stype == "statistics":
            if kline == "distribution":
                knox_data[-1].stats = sdata
                sdata = []
                sctr = 0
                stype = "distribution"
                continue
            next_row = np.array([int(float(x)) for x in kline.split()])
            if sctr == 0:
                sdata = next_row
            else:
                sdata = np.vstack([sdata, next_row])
            sctr += 1
        elif stype == "distribution":
            sctr += 1
            
            # [...
            # start of dist for 1st cell in row
            if kline[0]=="[":
                kline = kline[1:]
                
            # ...] [...
            # end of one dist, start of another, in row
            elif "] [" in kline:
                kline = " ".join(kline.split("] ["))
                
            # ...]
            # end of dist for last cell in row
            elif kline[-1] == "]":
                kline = kline[:-1]
                
            # ...
            # typical line of numbers in middle of a dist
            else:
                pass
            
            sdata += [int(float(x)) for x in kline.split()]
        

print("\# events per experiment")
for exp in knox_data:
    print("{}:{}".format(exp.start_date,exp.num_events))

bandwidth_selections = ["along_axis", "contig_to_axis","contig_anywhere"]
bandwidth_pairs_dict = defaultdict(list)
for exp_num, exp in enumerate(knox_data):
    #print("Exp num {}".format(exp_num))
    #print(exp.ratios)
    
    if exp_num<15:
        plot_signif_knox_ratios(exp, sbin_size, tbin_size, 0.05)
    
    
    for band_sel in bandwidth_selections:
        bandwidth_pairs_dict[band_sel].append(get_bandwidths_from_knox(exp.pvals, selection=band_sel))




fig, ax = plt.subplots(figsize=(12,4))
for i, band_sel in enumerate(bandwidth_selections):
    ax.plot(np.linspace(2001, 2017.5, 34), [(x[0]+1+(.1*i))*sbin_size for x in bandwidth_pairs_dict[band_sel]])
ax.legend(bandwidth_selections)
ax.set_title("Spatial bandwidths determined by Knox, in meters")


fig, ax = plt.subplots(figsize=(12,4))
for i, band_sel in enumerate(bandwidth_selections):
    ax.plot(np.linspace(2001, 2017.5, 34), [(x[1]+1+(.1*i))*tbin_size/7 for x in bandwidth_pairs_dict[band_sel]])
ax.legend(bandwidth_selections)
ax.set_title("Temporal bandwidths determined by Knox, in weeks")


