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

"""
shorthandToTimeDelta
Convert a shorthand string like "3D", "2W", "6M", or "1Y" to
 an np.timedelta64 object.
"""
def shorthandToTimeDelta(step):
    return np.timedelta64(*standardizeTimeStep(step))

"""
multiply_shorthand_timestep
"""
def multipy_shorthand_timestep(step, mult):
    step_num = int(step[:-1])
    step_unit = step[-1].upper()
    return f"{step_num * mult}{step_unit}"

"""
generateDateRange
"""
# start: "2018-09-01"
# end: "2019-03-01"
# step: "3D", "2W", "6M", "1Y", etc
# output is always array containing type datetime64[D]
def generateDateRange(start=None, 
                      end=None, 
                      step=None,
                      num=None):
    
    # 3 values are needed to generate the 4th
    if [start,end,step,num].count(None) != 1:
        print("Error! Unexpected number of arguments for generateDateRange.")
        print("Expected exactly 3 of these arguments: ")
        print(f" start: {start}")
        print(f" end: {end}")
        print(f" step: {step}")
        print(f" num: {num}")
        sys.exit(1)
    
    # Currently, step is required
    if step==None:
        print("Error! No step argument given to generateDateRange")
        sys.exit(1)
    step_num, step_unit = standardizeTimeStep(step)
    
    start_datetime = None
    end_datetime = None
    
    
    # If no given start, generate it via end and step and num
    if start == None:
        end_datetime = np.datetime64(end)
        start_datetime = generateEarlierDate(end_datetime, 
                                        multipy_shorthand_timestep(step, num))
    # If no given end, generate it via start and step and num
    elif end == None:
        start_datetime = np.datetime64(start)
        end_datetime = generateLaterDate(start_datetime, \
                                        multipy_shorthand_timestep(step, num))
    # If given start and end, set them as the start and end datetimes
    else:
        start_datetime = np.datetime64(start)
        end_datetime = np.datetime64(end)
    
    # Generate date range via start and end and step
    date_range = np.arange(start_datetime, \
                           end_datetime, \
                           step=np.timedelta64(step_num, step_unit), \
                           dtype=f"datetime64[{step_unit}]")
    
    print(f"date_range: {date_range}")
    
    
    
    date_range_array = np.array(date_range, dtype="datetime64[D]")
    if num!= None and len(date_range_array) != num:
        print("Warning! Different number of experiments than expected.")
        print(f"  Specified number: {num}")
        print(f"  True number: {len(date_range_array)}")
        print(f"  Earliest: {date_range_array[0]}")
        print(f"  Latest: {date_range_array[-1]}")
    return date_range_array



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
# Note that this returns points whose times are GREATER OR EQUAL to 
#  the start time, but STRICTLY LESS THAN the end time.
def getTimedPointsInTimeRange(points, start_time, end_time):
    return points[(points.timestamps >= start_time)
                  & (points.timestamps < end_time)]


def getSixDigitDate(in_date):
    str_date = str(in_date)
    return "".join(x[-2:] for x in str_date.split("-"))
    #return "{:02}{:02}{:02}".format(in_date.year%100, in_date.month, in_date.day)

