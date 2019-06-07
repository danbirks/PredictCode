# -*- coding: utf-8 -*-
"""
Created on Tue May 14 10:32:44 2019

@author: lawdfo
"""


import numpy as np
import sys
import csv
import os


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


def saveCrimeSubsetCsv(infile, start_time, end_time, crime_list, outfile=None, nrec=None):
    outdir = os.path.dirname(infile)
    if outfile == None:
        crime_str = "-".join(crime_list)
        outfilebase = "_".join(["chi", "all", crime_str, getSixDigitDate(start_time), getSixDigitDate(end_time)])
        if nrec != None:
            outfilebase += "_" + str(nrec)
        outfilebase += ".csv"
        outfile = os.path.join(outdir, outfilebase)
    print(outfile)
    row_ctr = 0
    with open(outfile, "w") as outf:
        outf_writer = csv.writer(outf, lineterminator='\n')
        with open(infile) as inf:
            inf_reader = csv.reader(inf)
            outf_writer.writerow(next(inf_reader))
            for row in inf_reader:
                if row[5] in crime_list:
                    row_m, row_d, row_y = row[2].split()[0].split("/")
                    row_date = np.datetime64("-".join([row_y, row_m, row_d]))
                    if start_time <= row_date and row_date <= end_time:
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