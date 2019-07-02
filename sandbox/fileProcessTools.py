# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 14:22:19 2019

@author: lawdfo
"""

import numpy as np
import csv
import sys
import os

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


def loadGenericData(filepath):
    
    
    
    import open_cp.sources.chicago as chicago
    
    
    sys.path.insert(0, os.path.abspath(".."))
    # Elements from PredictCode's custom "open_cp" package
    import open_cp
    import open_cp.geometry
    
    
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
    
    
    ### OBTAIN GRIDDED REGION
    
    # Obtain polygon shapely object for region of interest
    region_polygon = chicago.get_side(chicago_side)
    
    # Obtain data set within relevant region
    points_crime_region = open_cp.geometry.intersect_timed_points(points_crime, region_polygon)
    
    return points_crime_region




"""
datadir = os.path.join("..", "..", "Data")
naivename = "ratesummary_xsr_nai_190701_130101_171231_1M.csv"
naivepath = os.path.join(datadir, naivename)
phsname = "ratesummary_xsr_phs_190701_130101_171231_1M.csv"
phspath = os.path.join(datadir, phsname)
newcsvname = "ratesummary_xsr_naiVphs_190701_130101_171231_1M.csv"
newcsvpath = os.path.join(datadir, newcsvname)
combineNaivePhsCsv(naivepath, phspath, newcsvpath)
"""