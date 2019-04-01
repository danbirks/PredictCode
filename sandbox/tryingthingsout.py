# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 13:42:31 2019

@author: Dustin
"""

import sys
import matplotlib.pyplot as plt
import os, csv, lzma
import numpy as np
import open_cp.sources.chicago
import geopandas as gpd
import pyproj
import shapely.geometry

print("Establishing basics of os.path")
print(os.path.sep)
print(os.path.join("a","b"))
print("Trying with just __file__")
print(__file__)
print(os.path.split(__file__))
print(os.path.split(__file__)[0])
print(os.path.join(os.path.split(__file__)[0],"uk_police.csv"))
print("Trying with abspath")
print(os.path.abspath(__file__))
print(os.path.split(os.path.abspath(__file__)))
print(os.path.split(os.path.abspath(__file__))[0])
print(os.path.join(os.path.split(os.path.abspath(__file__))[0],"uk_police.csv"))
sys.exit(0)

datadir = os.path.join("..", "..", "Data")
open_cp.sources.chicago.set_data_directory(datadir)

def get_csv_data(shortfilename="chicago.csv"):
    filename = os.path.join(datadir, shortfilename)
    if shortfilename.endswith(".csv.xz"):
        with lzma.open(filename, "rt") as f:
            yield from csv.reader(f)
    elif shortfilename.endswith(".csv"):
        with open(filename, "rt") as f:
            yield from csv.reader(f)
    else:
        yield None

rows = get_csv_data("chicago.csv")
print(next(rows))
print(next(rows))


polygon = open_cp.sources.chicago.get_side("South")
#polygon
#print(polygon)


frame = gpd.GeoDataFrame({"name":["South Side"]})
frame.geometry = [polygon]
frame.crs = {"init":"epsg:2790"}
print(frame)
display(frame)
#display(HTML(frame.to_html()))
