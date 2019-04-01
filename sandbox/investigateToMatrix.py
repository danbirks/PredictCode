# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 14:37:05 2019

@author: Dustin
"""


# Some fairly standard modules
import os, csv, lzma
import numpy as np
import matplotlib.pyplot as plt
import scipy

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


import open_cp
import open_cp.naive as naive

np.random.seed(1)

import open_cp.sources.ukpolice as ukpolice
points = ukpolice.default_burglary_data()

print("ukpolice, len-points")
print(inspect.getfile(ukpolice))
print(len(points.timestamps))

# Use pyproj to make a more properly projected visualization of the data
projected_points = open_cp.data.points_from_lon_lat(points, epsg=7405)
points = projected_points
bbox = points.bounding_box
fig, ax = plt.subplots(figsize=(10, 10 * bbox.aspect_ratio))
ax.scatter(points.xcoords, points.ycoords, s=10, alpha=0.2)
region = open_cp.RectangularRegion(bbox.xmin,bbox.xmax, bbox.ymin,bbox.ymax)

print("bbox, bbox type, region, region type")
print(bbox)
print(type(bbox))
print(region)
print(type(region))
print("")


predictor = naive.ScipyKDE()
predictor.data = points
prediction = predictor.predict()

print("predictor, points, prediction")
print(type(predictor))
print(type(points))
print(type(prediction))
print("prediction.samples, .cell_width, .cell_height, .xoffset, .yoffset")
print(prediction.samples)
print(prediction.cell_width)
print(prediction.cell_height)
print(prediction.xoffset)
print(prediction.yoffset)
print("")



import numpy as _np


def from_continuous_prediction_region(prediction, region, cell_width, cell_height=None):
    """Construct an instance from an instance of
    :class:`ContinuousPrediction` using the region and passed cell sizes.

    :param prediction: An instance of :class:`ContinuousPrediction` to
      sample from
    :param region: The :class:`RectangularRegion` the grid
    :param cell_width: Width of each cell in the resulting grid
    :param cell_height: Optional; height of each cell in the resulting
      grid; defaults to `cell_width`
    """
    if cell_height is None:
        cell_height = cell_width
    width = int(_np.rint((region.xmax - region.xmin) / cell_width))
    height = int(_np.rint((region.ymax - region.ymin) / cell_height))
    print("Width = " + str(width))
    print("Height = " + str(height))
    print(cell_width)
    print(cell_height)
    print(region.xmin)
    print(region.ymin)
    print("prediction samples:")
    print(prediction.samples)
    print("Now to make newpred via rebase")
    newpred = prediction.rebase(cell_width, cell_height, region.xmin, region.ymin)
    print("newpred type:")
    print(type(newpred))
    print("newpred samples:")
    print(newpred.samples)
    print("prediction samples:")
    print(prediction.samples)
    sys.exit(0)
    return from_continuous_prediction(newpred, width, height)


def from_continuous_prediction(prediction, width, height):
    """Construct an instance from an instance of
    :class:`ContinuousPrediction` using the grid size and offset specified
    in that instance.  This is more efficient as we sample each grid cell
    once and then store the result.

    :param prediction: An instance of ContinuousPrediction to sample from
    :param width: Width of the grid, in number of cells
    :param height: Height of the grid, in number of cells
    """
    print("gonna make matrix...")
    matrix = to_matrix(prediction, width, height)
    print("gonna make GridPredictionArray...")
    return GridPredictionArray(prediction.cell_width, prediction.cell_height,
        matrix, prediction.xoffset, prediction.yoffset)


def to_matrix(self, width, height):
    """Sample the risk at each grid point from `(0, 0)` to
    `(width-1, height-1)` inclusive.  Optimised."""
    if self.samples < 0:
        return self._to_matrix_grid(width, height)
    matrix = _np.empty((height, width))
    for gy in range(height):
        if gy%1==0:
            print("gy = "+str(gy))
        print("a")
        y = (gy + _np.random.random(size=self.samples * width)) * self.cell_height + self.yoffset
        # 0,1,...,width-1, 0,1,...,width-1  with the block repeated self.sample times
        print("b")
        gx = _np.broadcast_to(_np.arange(width), (self.samples, width)).ravel()
        print("c")
        x = (gx + _np.random.random(self.samples * width)) * self.cell_width + self.xoffset
        print("d")
        print("let's make a risk array...")
        print(type(x))
        print(len(x))
        print(type(y))
        print(len(y))
        selfriskarray = _risk_array(self,x, y)
        print("e")
        selfsampleswidth = (self.samples, width)
        print("f")
        reshaped = _np.reshape(selfriskarray, selfsampleswidth)
        print("g")
        matrix[gy] = _np.mean(reshaped, axis=0)
        print("h")
    return matrix


def _risk_array(self, x, y):
    # Like `return self.risk(x,y)` but do in blocks of at most 50 to avoid
    # excessive memory usage
    assert len(x.shape) == 1
    out = _np.empty_like(x)
    offset = 0
    length = x.shape[0]
    print("length = " + str(length))
    while offset < length:
        if offset%10000==0:
            print(offset)
        end = min(offset + 50, length)
        xx, yy = x[offset : end], y[offset : end]
        out[offset : end] = self.risk(xx, yy)
        offset = end
    return out


gridpred = from_continuous_prediction_region(prediction, region, 2500)
print("all done! ha ha")

"""
So here's what happens

from_continuous_prediction_region(prediction {type KernelRiskPredictor}, 
                                  region {type RectangularRegion},
                                  2500
                                  )

Within that:
    newpred {type ContinuousPrediction}
    =
    prediction.rebase(cell_width {2500}, 
                      cell_height {2500}, 
                      region.xmin, 
                      region.ymin
                      )
    Here's what that does
        instance = ContinuousPrediction(cell_width, 
                                        cell_height, 
                                        xoffset,
                                        yoffset, 
                                        samples
                                        )
        instance.risk = self.risk
        return instance

And then it returns this:
    from_continuous_prediction(newpred {ContinuousPrediction}, 
                               width, 
                               height
                               )
    Note: "width" is NOT the full range of the variable (x or y);
        It is that range DIVIDED BY the cell_width. (Same with height.)
        So "width" is more like, how many cells there are in each row

So now let's figure out what from_continuous_prediction does.

def from_continuous_prediction(prediction {"newpred", type ContinuousPrediction}, 
                               width, {# cells in a row}
                               height {# cells in a col}
                               ):
    matrix = to_matrix(prediction {"newpred", type ContinuousPrediction}, 
                       width,  {# cells in a row}
                       height  {# cells in a col}
                       )
    return GridPredictionArray(prediction.cell_width, 
                               prediction.cell_height, 
                               matrix, 
                               prediction.xoffset, 
                               prediction.yoffset)

So what does to_matrix do?
def to_matrix(self, {"newpred", type ContinuousPrediction}
              width, {# cells in a row}
              height {# cells in a col}
              ):
    # The comments say:
    # Sample the risk at each grid point from `(0, 0)` to
    #    (width-1, height-1)` inclusive.  Optimised.
    if self.samples < 0:
        return self._to_matrix_grid(width, height)
    matrix = _np.empty((height, width))
    for gy in range(height):
        if gy%1==0:
            print("gy = "+str(gy))
        y = (gy + _np.random.random(size=self.samples * width)) * self.cell_height + self.yoffset
        # 0,1,...,width-1, 0,1,...,width-1  with the block repeated self.sample times
        gx = _np.broadcast_to(_np.arange(width), (self.samples, width)).ravel()
        x = (gx + _np.random.random(self.samples * width)) * self.cell_width + self.xoffset
        # HERE IS THE LINE THAT TAKES FOREVER
        selfriskarray = _risk_array(self,x, y)
        selfsampleswidth = (self.samples, width)
        reshaped = _np.reshape(selfriskarray, selfsampleswidth)
        matrix[gy] = _np.mean(_np.reshape(self._risk_array(x, y), (self.samples, width)), axis=0)
    return matrix

"""




