# -*- coding: utf-8 -*-
"""
Created on Tue May 14 10:32:44 2019

@author: lawdfo
"""


import numpy as np
import sys


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

