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
from collections import defaultdict
import statistics
import time
from string import ascii_uppercase as asc_up

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
#import open_cp.sources.chicago as chicago
#import open_cp.retrohotspot as retro
#import open_cp.prohotspot as phs
import open_cp.knox


from riskModelsGeneric import splitCommaArgs, \
                              loadGenericData
from crimeRiskTimeTools import generateDateRange, \
                               getTimedPointsInTimeRange, \
                               generateLaterDate, \
                               generateEarlierDate



_cdict = {'red':   [(0.0,  1.0, 1.0),
                   (1.0,  1.0, 1.0)],
         'green': [(0.0,  1.0, 1.0),
                   (1.0,  0.0, 0.0)],
         'blue':  [(0.0,  0.2, 0.2),
                   (1.0,  0.2, 0.2)]}

yellow_to_red = matplotlib.colors.LinearSegmentedColormap("yellow_to_red", _cdict)




"""
KnoxEntry

Object that stores the following values for a Knox run:
    start_date
    end_date
    num_events
    sbins
    tbins
    stats
    medians
    pvals
    ratios
"""
class KnoxEntry:
    
    def __init__(self, 
                 start_date = None, 
                 end_date = None, 
                 window_size = None, 
                 num_events = -1, 
                 sbins = [], 
                 tbins = []):
        self.start_date = start_date
        self.end_date = end_date
        self.window_size = window_size
        self.num_events = num_events
        self.sbins = sbins
        self.tbins = tbins
        self.stats = []
        self.medians = []
        self.pvals = []
        self.ratios = []
        





"""
getKnoxResult
Creates Knox object from open_cp.knox, sets the spatial and temporal bin
 parameters along with the input data, and performs the requisite Knox
 calculations.
"""
def getKnoxResult(data, num_iter, sbins, tbins, tbin_unit="days"):
    
    knox = open_cp.knox.Knox()
    knox.set_time_bins(tbins, unit=tbin_unit)
    knox.space_bins = sbins
    knox.data = data
    result = knox.calculate(iterations=num_iter)
    return result
    

"""
makeBins

Given a size for a "bin" and a desired number of bins, creates that many bins
 of that size, starting with (0, size).
For example, makeBins(10,3) will result in [(0, 10), (10, 20), (20, 30)]
Note that the returned object is a list containing 2-long tuples.
"""
def makeBins(size, num):
    return list( (i*size, (i+1)*size) for i in range(num))







"""
knox_ratio

Calculate the Knox ratios, which are just the Knox statistics divided by the
 median of the distribution.
"""
def knox_ratio(knox_statistic, distribution):
    """As in the paper, compute the ratio of the statistic to the median
    of the values in the distribution"""
    #d = np.array(distribution)
    #d.sort()
    #return statistic / d[len(d)//2]
    return knox_statistic / statistics.median(distribution)



"""
significant_cells

Return array of booleans representing whether each pvalue is less than the
 specified significance threshold
"""
def significant_cells(pvalue_array, sig_thresh=0.05):
    return pvalue_array < sig_thresh


"""
contiguous_cells

Algorithm for determining a contiguous region of cells
"""
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


"""
get_bandwidths_from_knox

Returns a particular time bandwidth and space bandwidth, according to one
 of 3 possible selection methods:
     - contig_to_axis: largest significant value that is along the axis (i.e.
                        is also significant alongside lowest value for other
                        dimension), and can be reached from a contiguous path
                        of significant pairs from the origin
     - contig_anywhere: largest significant value that can be reached farthest
                         from the origin via a contiguous paths of other
                         significant pairs
     - along_axis: largest significant value that can be reached from the
                    origin via a vertical or horizontal contiguous path of
                    significant pairs
"""
def get_bandwidths_from_knox(pvalue_array, selection="contig_to_axis", sig_thresh=0.05):
    if selection not in ["contig_to_axis","contig_anywhere","along_axis"]:
        print("Error, unrecognized selction type: {}".format(selection))
        sys.exit(1)
    signif_array = significant_cells(pvalue_array, sig_thresh=sig_thresh)
    array_dims = np.shape(signif_array)
    
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


"""
get_signif_thresh_index

Given a value (expected 0.0-1.0) and a list of thresholds (also expected in
that range), return the index of the first threshold value that is greater
than the input value.

For example, given a threshold list of [0.02, 0.15, 0.77], an input value
of 0.01 would return 0, 0.08 would return 1, 0.5 would return 2, and 0.99
would return 3.
"""
def get_signif_thresh_index(p_val, thresh_list):
    thresh_list = sorted(thresh_list)
    t_index = 0
    while t_index < len(thresh_list) and p_val >= thresh_list[t_index]:
        t_index += 1
    return t_index

"""
plot_signif_knox_ratios


"""
def plot_signif_knox_ratios(knox_entry: KnoxEntry, 
                            p_thresh=[0.05], 
                            file_path = None):
    
    p_thresh = sorted(p_thresh)
    num_thresh = len(p_thresh)
    cell_texts = [asc_up[x]*(num_thresh-x) for x in range(num_thresh)]
    
    mappable = plt.cm.ScalarMappable(cmap=yellow_to_red)
    mappable.set_array(np.ravel(knox_entry.ratios))
    #mappable.autoscale()
    mappable.set_clim(vmin=0.5, vmax=2.0)
    fig, ax = plt.subplots(figsize=(12,4))
    
    #array_dims = np.shape(knox_entry.pvals)
    
    sbin_size = knox_entry.sbins[0][1] - knox_entry.sbins[0][0]
    tbin_size = knox_entry.tbins[0][1] - knox_entry.tbins[0][0]
    
    xmin = knox_entry.sbins[0][0]
    xmax = knox_entry.sbins[-1][1]
    ax.set(xlim=(xmin, xmax), xlabel="Distance in metres")
    
    ymin = knox_entry.tbins[0][0]
    ymax = knox_entry.tbins[-1][1]
    ax.set(ylim=(ymin, ymax), ylabel="Time in days")
    
    ax.set_title("Knox, {} events from {} to {}, p={}".format(knox_entry.num_events, knox_entry.start_date, knox_entry.end_date, p_thresh))
    
    for (tbin_index,sbin_index), pval in np.ndenumerate(knox_entry.pvals):
        thresh_index = get_signif_thresh_index(pval, p_thresh)
        # Don't plot values that are not significant
        if thresh_index >= num_thresh:
            continue
        # Make a rectangular patch at the position corresponding to the bins,
        # and color it via "fc" corresponding to its Knox ratio
        sbin_val = knox_entry.sbins[sbin_index][0]
        tbin_val = knox_entry.tbins[tbin_index][0]
        p = matplotlib.patches.Rectangle(
             (sbin_val, tbin_val), 
             sbin_size, 
             tbin_size, 
             fc=mappable.to_rgba(knox_entry.ratios[tbin_index,sbin_index]))
        ax.add_patch(p)
        ax.text(sbin_val + (sbin_size * 0.5), 
                tbin_val + (tbin_size * 0.5), 
                cell_texts[thresh_index], 
                horizontalalignment='center',
                verticalalignment='center',
                )
    
    cbar = fig.colorbar(mappable, orientation="vertical")
    cbar.set_label("Knox ratio")
    
    
    if file_path != None:
        fig.savefig(file_path)
    



"""
get_knox_data_from_file

Read in the output from a Knox run with an expected format,
 return usable data.
"""
def get_knox_data_from_file(knox_file_path, exp_limit=0):
    
    
    
    # Expected file format: For each experiment,
    #  - start date
    #  - end date
    #  - time between start and end date, in shorthand
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
    
    
    info_castings = [np.datetime64, 
                     np.datetime64, 
                     str, 
                     int, 
                     None, 
                     eval, 
                     None, 
                     eval]
    
    
    
    with open(knox_file_path) as kf:
        exp_num = -1
        stype = "info"
        sctr = 0
        sdata = []
        knox_data = []
        for lnum, kline in enumerate(kf):
            kline = kline.strip()
            
            
            if stype == "info":
                # If we're done with the info section, store the info we read
                if kline == "Knox Statistics":
                    knox_data.append(KnoxEntry(*sdata))
                    sdata = []
                    sctr = 0
                    stype = "Knox Statistics"
                    continue
                cast_type = info_castings[sctr]
                if cast_type != None:
                    try:
                        sdata.append(cast_type(kline))
                    except:
                        print(f"Error, info section incorrect in part {exp_num}?")
                        sys.exit(1)
                sctr += 1
                continue
            
            
            elif stype == "Knox Statistics":
                if kline == "Monte Carlo Medians":
                    knox_data[-1].stats = sdata
                    sdata = []
                    sctr = 0
                    stype = "Monte Carlo Medians"
                    continue
                next_row = np.array([int(float(x)) for x in kline.split()])
                if sctr == 0:
                    sdata = next_row
                else:
                    sdata = np.vstack([sdata, next_row])
                sctr += 1
            
            
            
            elif stype == "Monte Carlo Medians":
                if kline == "P Values":
                    knox_data[-1].medians = sdata
                    sdata = []
                    sctr = 0
                    stype = "P Values"
                    continue
                next_row = np.array([float(x) for x in kline.split()])
                if sctr == 0:
                    sdata = next_row
                else:
                    sdata = np.vstack([sdata, next_row])
                sctr += 1
            
            
            
            elif stype == "P Values":
                if kline == "":
                    knox_data[-1].pvals = sdata
                    sdata = []
                    sctr = 0
                    stype = "info"
                    exp_num += 1
                    if exp_limit > 0 and exp_num >= exp_limit:
                        print(f"Warning! Reached experiment limit of {exp_limit}")
                        break
                    continue
                next_row = np.array([float(x) for x in kline.split()])
                if sctr == 0:
                    sdata = next_row
                else:
                    sdata = np.vstack([sdata, next_row])
                sctr += 1
                
        
        
        for exp_index, exp in enumerate(knox_data):
            #print(exp.stats)
            #print(exp.medians)
            with np.errstate(divide='ignore', invalid='ignore'):
                knox_data[exp_index].ratios = np.true_divide(exp.stats,
                                                     exp.medians, 
                                                     )
            
            knox_data[exp_index].ratios = \
                np.nan_to_num(knox_data[exp_index].ratios)
            #print(knox_data[exp_index].ratios)
    
    return knox_data
    
    
    





def make_knox_info_file(datadir, 
                        in_csv_file_name, 
                        out_knox_file_name, 
                        geojson_file_name, 
                        local_epsg_in, 
                        crime_types, 
                        num_knox_iterations, 
                        knox_sbin_size, 
                        knox_sbin_num, 
                        knox_tbin_size, 
                        knox_tbin_num, 
                        earliest_exp_time, 
                        num_exp, 
                        time_step, 
                        time_len, 
                        csv_date_format = "%m/%d/%Y %I:%M:%S %p", 
                        csv_longlat = False, 
                        csv_epsg = None, 
                        csv_infeet = True, 
                        csv_col_names = None, 
                        ):
    
    
    
    # Normalised and derived parameters
    
    # Normalised data directory
    datadir = os.path.expanduser(os.path.normpath(datadir))
    
    # Full paths to files
    in_csv_full_path = os.path.join(datadir, in_csv_file_name)
    in_geojson_full_path = os.path.join(datadir, geojson_file_name)
    
    # Set of relevant crime types in the data
    crime_type_set = set(splitCommaArgs(crime_types))
    
    # Spatial and temporal bandwidth bins
    knox_sbins = makeBins(knox_sbin_size, knox_sbin_num)
    knox_tbins = makeBins(knox_tbin_size, knox_tbin_num)
    
    earliest_start_time = generateEarlierDate(earliest_exp_time, time_len)
    print(f"First time window is from \
{earliest_start_time} to {earliest_exp_time}")
    start_times = generateDateRange(start=earliest_start_time, 
                                    step=time_step, 
                                    num=num_exp)
    
    
    
    out_file_path = os.path.join(datadir, out_knox_file_name)
    
    print(f"outfile: {out_file_path}")
    
    
    
    
    # Obtain crime data points, and region polygon
    
    # Obtain all crimes (of relevant types) from input data
    points_crime = loadGenericData(in_csv_full_path, 
                                   crime_type_set=crime_type_set, 
                                   date_format_csv = csv_date_format, 
                                   longlat=csv_longlat, 
                                   epsg = csv_epsg, 
                                   infeet=csv_infeet, 
                                   col_names = csv_col_names
                                   )
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
    points_crime_region = open_cp.geometry.intersect_timed_points(points_crime, region_polygon)
    
    total_num_events = len(points_crime_region.timestamps)
    
    print(f"Successfully obtained data, with {total_num_events} events.")
    
    
    
    
    
    
    # Do Knox runs and store info in file
    
    print(f"Opening file {out_file_path} for writing.")
    with open(out_file_path,"w") as fout:
        
        chkpt_0 = time.time()
        for exp_index, start_time in enumerate(start_times):
            
            chkpt_1 = time.time()
            
            end_time = generateLaterDate(start_time, time_len)
            
            print(f"Time span: {start_time} to {end_time}")
            
            ### SELECT TRAINING DATA
            
            chkpt_2 = time.time()
            print(f"Getting data subset...")
            # Get subset of data for training
            points_crime_region_train = getTimedPointsInTimeRange(points_crime_region, 
                                                              start_time, 
                                                              end_time)
            print(f"...Got data subset. ({time.time()-chkpt_2:.4f})")
            
            
            
            num_events = len(points_crime_region_train.timestamps)
            
            print(f"Number of events in timespan: {num_events}")
            
            chkpt_3 = time.time()
            print("Calculating Knox...")
            knox_result = getKnoxResult(points_crime_region_train, 
                                        num_knox_iterations, 
                                        knox_sbins, 
                                        knox_tbins)
            print(f"...Calculated Knox. ({time.time()-chkpt_3:.4f})")
            
            
            chkpt_4 = time.time()
            print(f"Writing to file {out_file_path} ...")
            fout.write(str(start_time))
            fout.write("\n")
            fout.write(str(end_time))
            fout.write("\n")
            fout.write(str(time_len))
            fout.write("\n")
            fout.write(str(num_events))
            fout.write("\n")
            fout.write("Spatial bins (columns):")
            fout.write("\n")
            fout.write(str(knox_sbins))
            fout.write("\n")
            fout.write("Temporal bins (rows):")
            fout.write("\n")
            fout.write(str(knox_tbins))
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
            print(f"...Wrote to file. ({time.time()-chkpt_4:.4f})")
            print(f"Time for this run: {time.time()-chkpt_1:.4f}")
    
    print(f"Number of runs: {len(start_times)}")
    print(f"Number of bins per run: {len(knox_sbins) * len(knox_tbins)}")
    print(f"Overall time: {time.time()-chkpt_0:.4f}")
    
    
    
    
    
    








def make_graphs_from_knox_file(datadir, 
                               knoxrun_file_name, 
                               signif_cutoff = [0.05], 
                               exp_limit = 0, 
                               jitter_factor = 0.02, 
                               knox_out_custom=None, 
                               graph_best_bands=False):
    
    # Derived parameters
    datadir = os.path.expanduser(os.path.normpath(datadir))
    knoxrun_file_path = os.path.join(datadir, knoxrun_file_name)
    
    # Ensure the significance cutoff is a list object, even if it only has
    #  one element.
    if type(signif_cutoff) != list:
        signif_cutoff = [signif_cutoff]
    
    
    # Retrieve data from saved file
    knox_data = get_knox_data_from_file(knoxrun_file_path, exp_limit=exp_limit)
    
    
    # 3 methods of selecting a bandwidth
    bandwidth_selections = ["along_axis", "contig_to_axis","contig_anywhere"]
    
    # Instantiate a dict that maps from each of the above bandwidth selection
    #  methods to the bandwidths they determine.
    bandwidth_pairs_dict = defaultdict(list)
    
    # Determine the size of each spatial and temporal bandwidth bin.
    # All bins of the same type are assumed to be the same size, so we just
    #  look at the sie of the first bin.
    sbin_size = knox_data[0].sbins[0][1] - knox_data[0].sbins[0][0]
    tbin_size = knox_data[0].tbins[0][1] - knox_data[0].tbins[0][0]
    
    
    for exp_num, exp in enumerate(knox_data):
        
        knox_grid_file_base = "knox_grid_"
        if knox_out_custom != None:
            knox_grid_file_base += knox_out_custom + "_"
        knox_grid_file_base += f"{exp.end_date}_{exp.window_size}.png"
        knox_grid_file_path = os.path.join(datadir, knox_grid_file_base)
        
        # Create grids that illustrate the statistically significant
        #  bandwidth bins, coloured based on their Knox ratios
        plot_signif_knox_ratios(exp, 
                                signif_cutoff, 
                                file_path=knox_grid_file_path)
        
        
        if graph_best_bands:
            # For each bandwidth selection method,
            for band_sel in bandwidth_selections:
                # Determine the largest significant bandwidth for space and time
                band_indices = get_bandwidths_from_knox(exp.pvals, selection=band_sel, sig_thresh=signif_cutoff[-1])
                # Store that result pair in the dictionary
                bandwidth_pairs_dict[band_sel].append(band_indices)
    
    
    
    if graph_best_bands:
        
        
        plot_file_base = "_"
        if knox_out_custom != None:
            plot_file_base += knox_out_custom + "_"
        plot_file_base += f"{knox_data[0].end_date}"
        plot_file_base += f"_{knox_data[0].window_size}"
        plot_file_base += f"_{len(knox_data)}"
        plot_file_base += f".png"
        
        
        xcoords = [exp.end_date for exp in knox_data]
        
        fig, ax = plt.subplots(figsize=(12,4))
        max_y = max([x[0]+1 for x in bandwidth_pairs_dict[band_sel] for band_sel in bandwidth_selections]) * sbin_size
        adjust_y = jitter_factor * max_y
        for i, band_sel in enumerate(bandwidth_selections):
            ycoords = [(x[0]+1)*sbin_size + (adjust_y * i) for x in bandwidth_pairs_dict[band_sel]]
            ax.scatter(xcoords, ycoords)
        ax.legend(bandwidth_selections)
        ax.set_title("Spatial bandwidths determined by Knox")
        ax.set(xlabel="End date of test")
        y_axis_min = min(-1,0-(max_y*jitter_factor))
        y_axis_max = max(1,max_y*(1+len(bandwidth_selections)*jitter_factor))
        ax.set(ylim=(y_axis_min, y_axis_max), ylabel="Meters")
        
        plot_file_path = os.path.join(datadir, 
                                      "knox_timeplot" + plot_file_base)
        fig.savefig(plot_file_path)
        
        
        fig, ax = plt.subplots(figsize=(12,4))
        max_y = max([x[1]+1 for x in bandwidth_pairs_dict[band_sel] for band_sel in bandwidth_selections]) * tbin_size
        adjust_y = jitter_factor * max_y
        for i, band_sel in enumerate(bandwidth_selections):
            ycoords = [(x[1]+1)*tbin_size + (adjust_y * i) for x in bandwidth_pairs_dict[band_sel]]
            ax.scatter(xcoords, ycoords)
        ax.legend(bandwidth_selections)
        ax.set_title("Temporal bandwidths determined by Knox")
        ax.set(xlabel="End date of test")
        y_axis_min = min(-1,0-(max_y*jitter_factor))
        y_axis_max = max(1,max_y*(1+len(bandwidth_selections)*jitter_factor))
        ax.set(ylim=(y_axis_min, y_axis_max), ylabel="Days")
        
        plot_file_path = os.path.join(datadir, 
                                      "knox_spaceplot" + plot_file_base)
        fig.savefig(plot_file_path)
    
    
    
    











"""
main:

If running this module as a script instead of importing its functions,
 this main function will perform a standard analysis with a set of
 default parameters.
"""
def main():
    
    
    dataset = "fantdur"
    
    if dataset == "chicago":
        
        # Location of data file
        datadir = "../../Data"
        
        # Input csv file name
        in_csv_file_name = "chi_all_s_BURGLARY_RES_010101_190101_stdXY.csv"
        
        # Output file for Knox info
        knox_file_name = "knoxtestingA.txt"
        
        # Geojson file
        geojson_file_name = "Chicago_South_Side_2790.geojson"
        
        crime_types = "BURGLARY"
        
        num_knox_iterations = 200
        
        #sbin in meters
        knox_sbin_size = 100
        knox_sbin_num = 10
        #tbin in days
        knox_tbin_size = 7
        knox_tbin_num = 8
        
        # Dates in format YYYY-MM-DD
        first_test_end = "2017-05-01"
        time_window_size = "4M"
        time_step = "1M"
        num_experiments = 4
        
        
        np.random.seed(seed=0)
        
        
        
        csv_date_format = "%m/%d/%Y %I:%M:%S %p"
        csv_longlat = False
        csv_epsg = None
        csv_infeet = True
        csv_has_header = True
        
    
    if dataset == "fantdur":
        
        # Location of data file
        datadir = "../../Data"
        
        # Input csv file name
        in_csv_file_name = "Fantasy-Durham-Data_std.csv"
        
        # Output file for Knox info
        knox_file_name = "knoxtestingFD3.txt"
        
        # Geojson file
        geojson_file_name = "Durham_27700.geojson"
        
        crime_types = "Burglary, Vehicle crime"
        
        num_knox_iterations = 200
        
        #sbin in meters
        knox_sbin_size = 200
        knox_sbin_num = 10
        #tbin in days
        knox_tbin_size = 7
        knox_tbin_num = 4
        
        # Dates in format YYYY-MM-DD
        first_test_end = "2019-09-01"
        time_window_size = "1M"
        time_step = "1W"
        num_experiments = 1
        
        
        np.random.seed(seed=0)
        
        
        csv_date_format = "%d/%m/%Y"
        csv_longlat = True
        csv_epsg = 27700
        csv_infeet = False
        csv_has_header = True
    
    
    
    
    make_knox_info_file(datadir=datadir, 
                        in_csv_file_name=in_csv_file_name, 
                        out_knox_file_name=knox_file_name, 
                        geojson_file_name=geojson_file_name, 
                        crime_types=crime_types, 
                        num_knox_iterations=num_knox_iterations, 
                        knox_sbin_size=knox_sbin_size, 
                        knox_sbin_num=knox_sbin_num, 
                        knox_tbin_size=knox_tbin_size, 
                        knox_tbin_num=knox_tbin_num, 
                        earliest_exp_time=first_test_end, 
                        num_exp=num_experiments, 
                        time_step=time_step, 
                        time_len=time_window_size, 
                        csv_date_format = csv_date_format, 
                        csv_longlat = csv_longlat, 
                        csv_epsg = csv_epsg, 
                        csv_infeet = csv_infeet, 
                        csv_has_header = csv_has_header, 
                        )
    
    print("Finished making Knox info file.")
    print("Next, reading the file and making graphs.")
    
    
    
    # Run evaluation function with default arguments
    
    # Additional input parameters
    
    # String for use in output image files
    knox_out_name = "fantasydurham"
    # Significance bands we're interested in
    signif_cutoff = [0.01, 0.05, 0.1]
    # If you only want to look at the first n results, set that here
    #  0 will look at all results.
    exp_limit = 0
    # Whether you want to generate scatterplots that attempt to pick the best
    #  spatial and temporal bandwidths, in a few different ways
    graph_best_bands = True
    
    make_graphs_from_knox_file(datadir, 
                               knox_file_name, 
                               signif_cutoff=signif_cutoff, 
                               exp_limit=exp_limit, 
                               knox_out_custom=knox_out_name, 
                               graph_best_bands=graph_best_bands)
    
    
    
    
    




if __name__ == "__main__":
    main()
