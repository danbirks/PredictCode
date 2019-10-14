# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 14:22:19 2019

@author: lawdfo

Purpose:
    Some generally useful functions for processing the relevant data files.



"""

import numpy as np
import csv
import sys
import os
import datetime


sys.path.insert(0, os.path.abspath(".."))
# Elements from PredictCode's custom "open_cp" package
import open_cp
import open_cp.geometry
import open_cp.sources.chicago as chicago
from open_cp.data import TimedPoints



"""
isSubsetFile
Purpose:
    Given 2 files, check whether each line of the first (shorter) file is
    present in the longer file, with order preserved but possibly with other
    lines interspersed throughout. This was useful at one point for checking
    that certain data subsets were being generated correctly.
Input: path to shorter file, path to longer file
Output: boolean of whether shorter is subset of longer.
        if false, also return earliest line in shorter that is not in longer.
"""
def isSubsetFile(shortfile, longfile):
    with open(shortfile) as sf:
        with open(longfile) as lf:
            for shortline in sf:
                foundit = False
                longline = lf.readline()
                while longline:
                    if longline == shortline:
                        foundit = True
                        break
                    longline = lf.readline()
                if foundit:
                    break
                return False, shortline
            return True, None

"""
combineNaivePhsCsv

"""
def combineNaivePhsCsv(naivefile, phsfile, newfile="temp.csv"):
    
    naivedata = dict()
    with open(naivefile) as naivef:
        
        naivereader = csv.reader(naivef)
        # Read each line of data
        for dataline in naivereader:
            cov = float(dataline[0])
            date = np.datetime64(dataline[1])
            hitrate = float(dataline[2])
            naivedata[(cov,date)] = hitrate
    
    with open(phsfile) as phsf:
        phsreader = csv.reader(phsf)
        with open(newfile,"w") as newf:
            writer = csv.writer(newf, delimiter=",", lineterminator="\n")
            for dataline in phsreader:
                cov = float(dataline[0])
                date = np.datetime64(dataline[1])
                timeband = int(dataline[2])
                distband = int(dataline[3])
                hitrate = float(dataline[4])
                naiverate = naivedata[(cov,date)]
                writer.writerow([cov, date, timeband, distband, hitrate, naiverate, hitrate-naiverate])
    
    sys.exit(0)



"""
loadGenericData

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
        # Remove header
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









def trialLoadGenericDataOLD(filepath):
    
    
    sys.path.insert(0, os.path.abspath(".."))
    # Elements from PredictCode's custom "open_cp" package
    import open_cp
    import open_cp.geometry
    import open_cp.sources.chicago as chicago
    from open_cp.data import TimedPoints
    
    
    crime_type_set = {"BURGLARY"}

    
    datadir = os.path.join("..", "..", "Data")
    #chicago_file_name = "chicago_all_old.csv"
    #chicago_file_name = "chi_all_s_BURGLARY_010101_190101.csv"
    chicago_file_name = "chi_all_s_BURGLARY_RES_010101_190101.csv"
    chicago_side = "South"
    chicago_load_type = "snapshot"
    if "all" in chicago_file_name:
        chicago_load_type = "all"
    chicago_file_path = os.path.join(datadir, chicago_file_name)
    # Chicago module requires this line to access some data
    chicago.set_data_directory(datadir)
    
    points_crime = chicago.load(chicago_file_path, crime_type_set, type=chicago_load_type)
    
    
    
    
    
    
    
    std_field_names = ["_DESC_FIELD", "_X_FIELD", "_Y_FIELD", "_TIME_FIELD"]
    
    custom_field_names = ['Primary Type', 'X Coordinate', 'Y Coordinate', 'Date']
    
    field_name_map = dict(zip(std_field_names, custom_field_names))
    
    date_format_csv = "%m/%d/%Y %I:%M:%S %p"
    def dt_convert(date_string, date_format=date_format_csv):
        return datetime.datetime.strptime(date_string, date_format)
    
    data = []
    
    with open(filepath) as f:
        reader = csv.reader(f)
        header = next(reader)
        header_num_map = dict(zip(header, range(len(header))))
        field_num_map = dict([(x, header_num_map[field_name_map[x]]) for x in std_field_names])
        
        
        for row in reader:
            desc = row[field_num_map["_DESC_FIELD"]].strip()
            if desc not in crime_type_set:
                continue
            x = row[field_num_map["_X_FIELD"]].strip()
            y = row[field_num_map["_Y_FIELD"]].strip()
            t = row[field_num_map["_TIME_FIELD"]].strip()
            
            data.append((dt_convert(t), float(x), float(y)))
        
    
    
    data.sort(key = lambda triple : triple[0])
    xcoords = np.empty(len(data))
    ycoords = np.empty(len(data))
    for i, triple in enumerate(data):
        xcoords[i], ycoords[i] = triple[1], triple[2]
    times = [triple[0] for triple in data]
    
    to_meters = True
    _FEET_IN_METERS = 3937 / 1200
    
    if to_meters:
        xcoords /= _FEET_IN_METERS
        ycoords /= _FEET_IN_METERS
    
    timedpoints = TimedPoints.from_coords(times, xcoords, ycoords)
    
    
    
    print(len(points_crime.timestamps))
    print(type(points_crime))
    
    print(len(timedpoints.timestamps))
    print(type(timedpoints))
    
    print(points_crime == timedpoints)
    print(points_crime.timestamps == timedpoints.timestamps)
    print(all(points_crime.timestamps == timedpoints.timestamps))
    
    print(points_crime.xcoords == timedpoints.xcoords)
    print(points_crime.ycoords == timedpoints.ycoords)
    print(all(points_crime.xcoords == timedpoints.xcoords))
    print(all(points_crime.ycoords == timedpoints.ycoords))
    
    
    print(points_crime.xcoords == timedpoints.xcoords)
    
    
    print(points_crime.bounding_box == timedpoints.bounding_box)
    print(points_crime.coords == timedpoints.coords)
    print(all((points_crime.coords == timedpoints.coords).flatten()))
    
    sys.exit(0)
    
    
    ### OBTAIN GRIDDED REGION
    
    # Obtain polygon shapely object for region of interest
    region_polygon = chicago.get_side(chicago_side)
    
    # Obtain data set within relevant region
    points_crime_region = open_cp.geometry.intersect_timed_points(points_crime, region_polygon)
    
    
    
    print(len(points_crime_region.timestamps))
    print(type(points_crime_region))
    
    sys.exit(0)
    
    return points_crime_region





def compareCoordConversion():
    
    import sys
    
    proj=None
    #epsg=7405
    epsg=3435
    
    xcoord_orig = 1185226
    ycoord_orig = 1864727
    lat_orig = 41.783960406
    long_orig = -87.596431561
    
    xcoord_orig = 1163097
    ycoord_orig = 1908381
    lat_orig = 41.904243234
    long_orig = -87.676340272
    
    
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
    
    #newx, newy = proj(lat_orig, long_orig)
    #print(newx,newy)
    newx, newy = proj(long_orig, lat_orig)
    print(newx,newy)
    
    _FEET_IN_METERS = 3937 / 1200
    newx = xcoord_orig / _FEET_IN_METERS
    newy = ycoord_orig / _FEET_IN_METERS
    print(newx,newy)
    
    
    print()
    sys.exit(0)
    return


def writeCsvCols(origcsvfile, colnums, newcsvfile=None):
    
    if newcsvfile == None:
        if origcsvfile[-4:]==".csv":
            newcsvfile = origcsvfile[:-4] + "_std.csv"
    
    with open(origcsvfile) as origf:
        csvreader = csv.reader(origf)
        
        with open(newcsvfile, "w") as newf:
            csvwriter = csv.writer(newf, delimiter=",", lineterminator="\n")
            
            # Read each line of data
            for dataline in csvreader:
                csvwriter.writerow([dataline[x] for x in colnums])
    
    return




"""
naivename = "ratesummary_xsr_nai_190701_130101_171231_1M.csv"
naivepath = os.path.join(datadir, naivename)
phsname = "ratesummary_xsr_phs_190701_130101_171231_1M.csv"
phspath = os.path.join(datadir, phsname)
newcsvname = "ratesummary_xsr_naiVphs_190701_130101_171231_1M.csv"
newcsvpath = os.path.join(datadir, newcsvname)
combineNaivePhsCsv(naivepath, phspath, newcsvpath)
"""

"""
chicago_file_name = "chi_all_s_BURGLARY_RES_010101_190101.csv"
datadir = os.path.join("..", "..", "Data")
chicago_file_path = os.path.join(datadir, chicago_file_name)
loadGenericData(chicago_file_path)
"""

"""
0 Date
1 Longitude
2 Latitude
3 Crime type
4 Location type
"""

#compareCoordConversion()


"""
# Converting chicago-all file into "standardised" format of just a few cols
datadir = os.path.join("..", "..", "Data")
origcsvname = "chi_all_s_BURGLARY_RES_010101_190101.csv"
origcsvfile = os.path.join(datadir, origcsvname)
writeCsvCols(origcsvfile, [2,20,19,5,7]) # using long and lat
#writeCsvCols(origcsvfile, [2,15,16,5,7]) # using XYcoords
"""


"""
# Confirming that use of long-lat and x-y have negligible difference
datadir = os.path.join("..", "..", "Data")
chicsvname_ll = os.path.join(datadir, "chi_all_s_BURGLARY_RES_010101_190101_std.csv")
chicsvname_xy = os.path.join(datadir, "chi_all_s_BURGLARY_RES_010101_190101_stdXY.csv")
ll_data = loadGenericData(chicsvname_ll, crime_type_set = {"BURGLARY"}, date_format_csv = "%m/%d/%Y %I:%M:%S %p", epsg = 3435, proj=None, longlat=True)
xy_data = loadGenericData(chicsvname_xy, crime_type_set = {"BURGLARY"}, date_format_csv = "%m/%d/%Y %I:%M:%S %p", epsg = 3435, proj=None, longlat=False)

ll_times = ll_data.timestamps
xy_times = xy_data.timestamps
ll_xs = ll_data.xcoords
ll_ys = ll_data.ycoords
xy_xs = xy_data.xcoords
xy_ys = xy_data.ycoords
print(len(ll_times))
print(len(xy_times))
if len(ll_times) != len(xy_times):
    sys.exit(1)
data_len = len(ll_times)
dist_diffs = []
for i in range(data_len):
    if ll_times[i] != xy_times[i]:
        print(" ".join([str(x) for x in [i, ll_times[i], xy_times[i]]]))
    dist_diffs.append(abs(ll_xs[i] - xy_xs[i]) + abs(ll_ys[i] - xy_ys[i]))
print("ok got it")
print(sum(dist_diffs))
print(sum(dist_diffs)/data_len)
sorted_dist_diffs = sorted(dist_diffs)
print(sorted_dist_diffs[:5])
print(sorted_dist_diffs[-5:])
print("ok done")
# Results:
# Average taxicab difference: .00004849564592251813
# Maximum taxicab difference: .00009684235556051135
"""


"""
# Testing contrived Durham data
datadir = os.path.join("..", "..", "Data")
#durfname = os.path.join(datadir, "chi_all_s_BURGLARY_RES_010101_190101_std.csv")
durfname = os.path.join(datadir, "contrived_durham_test.csv")
fake_durham_data = loadGenericData(durfname, crime_type_set = {"BURGLARY"}, date_format_csv = "%m/%d/%Y %I:%M:%S %p", epsg = 27700, proj=None, longlat=True)
fd_times = fake_durham_data.timestamps
fd_xs = fake_durham_data.xcoords
fd_ys = fake_durham_data.ycoords
fd_len = len(fd_times)
for i, d in enumerate(zip(fd_times, fd_xs, fd_ys)):
    print(f'{i} {d[0]} {d[1]} {d[2]}')
print("ok done")
"""

"""
# Converting durham.kml Force Boundaries file from
#   https://data.police.uk/data/boundaries/
# ogr2ogr -f GeoJSON durham.json durham.kml
import geopandas as gpd
datadir = os.path.join("..", "..", "Data")
durhamfilename = "durham.json"
durhamfilepath = os.path.join(datadir, durhamfilename)
frame = gpd.read_file(durhamfilepath)
frame.crs == {"init": "epsg:4326"}
frame.head()

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10,14))
_ = frame.plot(ax=ax)
_ = frame.plot(ax=ax, column="Name")

import tilemapbase
extent = tilemapbase.extent_from_frame(frame, 1000, 5)
fig, ax = plt.subplots(figsize=(14,15))
extent.plot(ax, tilemapbase.tiles.OSM, alpha=1)
frame.plot(ax=ax, column="region", alpha=0.5)
for _, row in frame.iterrows():
    x, y = row.geometry.centroid.coords[0]
    ax.text(x, y, row.area_numbe, horizontalalignment='center', verticalalignment='center')
ax.set_title("Chicago regions")
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
ax.set_aspect(1)
None
"""


datadir = os.path.join("..", "..", "Data")
#chicsvname_ll = os.path.join(datadir, "chi_all_s_BURGLARY_RES_010101_190101_std.csv")
chicsvname_xy = os.path.join(datadir, "chi_all_s_BURGLARY_RES_010101_190101_stdXY.csv")
chidata = loadGenericData(chicsvname_xy, longlat=False, infeet=True)
ts = chidata.timestamps
print(len(ts))



