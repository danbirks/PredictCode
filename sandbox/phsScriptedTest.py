# -*- coding: utf-8 -*-
"""
Created on Thu May  2 13:40:04 2019

@author: lawdfo
"""

# Import libraries we need
import csv
import open_cp.data
import geopandas as gpd
import pyproj
import dateutil.parser
import open_cp.scripted as scripted
import datetime


# Load the input data; see `scripted_intro.md`

def row_to_datetime(row):
    datetime_string = row[1]
    return dateutil.parser.parse(datetime_string)

chicago_projection = pyproj.Proj({"init" : "epsg:2790"})

def row_to_coords(row, projection=chicago_projection):
    x = float(row[5])
    y = float(row[4])
    return projection(x, y)

def load_points():
    with open("example.csv") as file:
        reader = csv.reader(file)
        header = next(reader)
        # Assume the header is correct
        times = []
        xcoords = []
        ycoords = []
        for row in reader:
            times.append(row_to_datetime(row))
            x, y = row_to_coords(row)
            xcoords.append(x)
            ycoords.append(y)
      
    # Maybe `times` is not sorted.
    times, xcoords, ycoords = open_cp.data.order_by_time(times, xcoords, ycoords)
      
    return open_cp.data.TimedPoints.from_coords(times, xcoords, ycoords)

def load_geometry():
    frame = gpd.read_file("SouthSide")
    return frame.geometry[0]




def stringToDatetime(s):
    return datetime.datetime(*(s.split("-")))


earliest_time = "2016-01-01"
start_time = "2016-10-01"
end_time = "2017-01-01"



# Perform the predictions; see `scripted_intro.md`

with scripted.Data(load_points, load_geometry,
        start=stringToDatetime(earliest_time)) as state:

    time_range = scripted.TimeRange(stringToDatetime(start_time),
            stringToDatetime(end_time), datetime.timedelta(days=1))

    state.add_prediction(scripted.NaiveProvider, time_range)

    state.score(scripted.HitRateEvaluator)
    state.score(scripted.HitCountEvaluator)

    state.process(scripted.HitRateSave("ratesTest.csv"))
    state.process(scripted.HitCountSave("countsTest.csv"))
