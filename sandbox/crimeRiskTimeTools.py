# -*- coding: utf-8 -*-
"""
Created on Tue May 14 10:32:44 2019

@author: lawdfo
"""


import numpy as np
import sys
import csv
import os


import open_cp
import open_cp.geometry
import open_cp.plot
import open_cp.sources.chicago as chicago

_day = np.timedelta64(1,"D")



"""
standardizeTimeStep(step)
Split a string of the form "3D", "2W", "6M", "1Y", etc, into an int and a
letter representing a unit of time. If the unit is "W" (weeks), it is instead
converted to "D" (days).
"""
def standardizeTimeStep(step):
    step_num = int(step[:-1])
    step_unit = step[-1].upper()
    if step_unit not in "YMWD":
        print("This function only recognizes Y/M/W/D as time units")
        sys.exit(1)
    # Weeks seem to cause weirdness, espeically when times get "rounded"
    #  to the nearest week. So here we just convert it to days.
    if step_unit == "W":
        step_num *= 7
        step_unit = "D"
    return step_num, step_unit

"""
shorthandToTimeDelta
Convert a shorthand string like "3D", "2W", "6M", or "1Y" to
 an np.timedelta64 object.
"""
def shorthandToTimeDelta(step):
    return np.timedelta64(*standardizeTimeStep(step))

#
# start: "2018-09-01"
# end: "2019-03-01"
# step: "3D", "2W", "6M", "1Y", etc
# output is always array containing type datetime64[D]
def generateDateRange(start, end, step):
    step_num, step_unit = standardizeTimeStep(step)
    start_datetime = np.datetime64(start)
    end_datetime = np.datetime64(end)
    date_range = np.arange(start_datetime, \
                           end_datetime, \
                           step=np.timedelta64(step_num, step_unit), \
                           dtype="datetime64[{}]".format(step_unit))
    return np.array(date_range, dtype="datetime64[D]")



def generateLaterDate(start, step):
    step_num, step_unit = standardizeTimeStep(step)
    converted_start = np.datetime64(start, step_unit)
    later_date = converted_start + np.timedelta64(step_num, step_unit)
    return np.datetime64(later_date, "D")

# generateEarlierDate is identical to generateLaterDate, just with a
#  negative step. Still nice to have a clearly-named function.
def generateEarlierDate(end, step):
    return generateLaterDate(end, "-"+step)


# Given a TimedPoints object that contains spatial coordinates paired
#  with timestamps, return a similar object that only contains the
#  data whose timestamps are within the given range
def getTimedPointsInTimeRange(points, start_time, end_time):
    return points[(points.timestamps >= start_time)
                  & (points.timestamps <= end_time)]


def getSixDigitDate(in_date):
    str_date = str(in_date)
    return "".join(x[-2:] for x in str_date.split("-"))
    #return "{:02}{:02}{:02}".format(in_date.year%100, in_date.month, in_date.day)


def saveCrimeSubsetCsv(infile, start_time, end_time, crime_list, loc_list=None, outfile=None, nrec=None, chicago_side=None):
    
    
    
    datadir = os.path.join("..", "..", "Data")
    chicago.set_data_directory(datadir)
    
    
    outdir = os.path.dirname(infile)
    if outfile == None:
        crime_str = "-".join(crime_list)
        if loc_list != None:
            crime_str += "_RES"
        chicago_side_short = "a"
        if chicago_side.upper() == "SOUTH":
            chicago_side_short = "s"
        outfilebase = "_".join(["chi", "all", chicago_side_short, crime_str, getSixDigitDate(start_time), getSixDigitDate(end_time)])
        if nrec != None:
            outfilebase += "_" + str(nrec)
        outfilebase += ".csv"
        outfile = os.path.join(outdir, outfilebase)
    print(outfile)
    row_ctr = 0
    with open(outfile, "w") as outf:
        outf_writer = csv.writer(outf, lineterminator='\n')
        
        
        
        # Obtain data as TimedPoints object, and corresponding csv rows
        points_crime, csv_rows, csv_header = chicago.load(infile, crime_list, type="all", withcsvrows=True)
        outf_writer.writerow(csv_header)
        points_csv_dict = dict(zip([tuple(p) for p in points_crime], csv_rows))
        
        # If a region of Chicago has been specified ("South"), do that here
        points_crime_region = points_crime
        if chicago_side != None:
            region_polygon = chicago.get_side(chicago_side)
            points_crime_region = open_cp.geometry.intersect_timed_points(points_crime, region_polygon)
        csv_rows_region = [points_csv_dict[tuple(p)] for p in points_crime_region]
        
        
        for row in csv_rows_region:
            if row[5] not in crime_list:
                print("Wait how'd that happen?")
                print(row)
                print(row[5])
                print(crime_list)
                sys.exit(1)
            
            if loc_list==None or row[7] in loc_list:
                
                
                row_m, row_d, row_y = row[2].split()[0].split("/")
                row_date = np.datetime64("-".join([row_y, row_m, row_d]))
                if start_time <= row_date and row_date < end_time:
                    outf_writer.writerow(row)
                    row_ctr += 1
                    if nrec != None and row_ctr >= nrec:
                        break
    print(f"Wrote header and {row_ctr} rows to {outfile}")
    return

def onetimemake2yrburglaryfile():
    datadir = os.path.join("..", "..", "Data")
    chicago_file_name = "chicago_all_old.csv"
    chicago_file_path = os.path.join(datadir, chicago_file_name)
    saveCrimeSubsetCsv(chicago_file_path, np.datetime64("2016-01-01"), np.datetime64("2018-01-01"), ["BURGLARY"])
    return

def onetimemakeburglaryfilenrecords(nrec):
    datadir = os.path.join("..", "..", "Data")
    chicago_file_name = "chicago_all_old.csv"
    chicago_file_path = os.path.join(datadir, chicago_file_name)
    saveCrimeSubsetCsv(chicago_file_path, np.datetime64("2016-01-01"), np.datetime64("2018-01-01"), ["BURGLARY"], nrec=nrec)
    return

def makeBurgSSXfile(dstart="2016-01-01", dend="2016-02-01", nrec=None):
    datadir = os.path.join("..", "..", "Data")
    chicago_file_name = "chicago_all_old.csv"
    chicago_file_path = os.path.join(datadir, chicago_file_name)
    saveCrimeSubsetCsv(chicago_file_path, np.datetime64(dstart), np.datetime64(dend), ["BURGLARY"], nrec=nrec, chicago_side="South")
    return

def makeBurgResSSXfile(dstart="2016-01-01", dend="2016-02-01", nrec=None):
    datadir = os.path.join("..", "..", "Data")
    chicago_file_name = "chicago_all_old.csv"
    chicago_file_path = os.path.join(datadir, chicago_file_name)
    saveCrimeSubsetCsv(chicago_file_path, np.datetime64(dstart), np.datetime64(dend), ["BURGLARY"], loc_list=["APARTMENT", "RESIDENCE", "RESIDENCE-GARAGE", "CHA APARTMENT"], nrec=nrec, chicago_side="South")
    return




#print("making 1")
#onetimemakeburglaryfilenrecords(1)
#print("made 1")
#print("making 2")
#onetimemakeburglaryfilenrecords(2)
#print("made 2")
#onetimemakeburglaryfilenrecords(3)
#print("made 3")
#onetimemakeburglaryfilenrecords(4)
#print("made 4")
#onetimemakeburglaryfilenrecords(5)
#print("made 5")
#onetimemakeburglaryfilenrecords(10)
#print("made 10")
# chi_s_BURGLARY_160101_160201_10 
#makeBurgSSXfile(dstart="2016-01-01", dend="2016-02-01", nrec=10)
#makeBurgSSXfile(dstart="2016-01-01", dend="2016-02-01", nrec=None)
#makeBurgSSXfile(dstart="2001-01-01", dend="2019-01-01", nrec=None)
#makeBurgSSXfile(dstart="2016-01-01", dend="2017-01-01", nrec=None)
#makeBurgResSSXfile(dstart="2001-01-01", dend="2019-01-01", nrec=None)

