# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 14:49:11 2019

@author: lawdfo
"""


#import sys
import os
import geopandas as gpd
import matplotlib.pyplot as plt
from itertools import product
from descartes import PolygonPatch
from shapely.geometry import Point, Polygon
import geojson
from collections import OrderedDict
import numpy




def plotPointsOnGrid(points, 
                     masked_grid, 
                     polygon, 
                     title=None, 
                     sizex=10, 
                     sizey=None, 
                     out_img_file_path=None):
    
    if sizey == None:
        sizey = sizex
    
    fig, ax = plt.subplots(figsize=(sizex,sizey))
    
    ax.add_patch(PolygonPatch(polygon, fc="none", ec="Black"))
    ax.add_patch(PolygonPatch(polygon, fc="Blue", ec="none", alpha=0.2))
    #ax.scatter(points.xcoords,
    #           points.ycoords,
    #           marker="+", color="red")
    
    xmin, ymin, xmax, ymax = polygon.bounds
    
    # Set the axes to have a buffer of 500 around the polygon
    ax.set(xlim=[xmin-500,xmax+500], ylim=[ymin-500,ymax+500])
    
    #pc = patches_from_grid(masked_grid)
    
    
    
    #ax.add_collection(PatchCollection(pc, facecolor="None", edgecolor="black"))
    
    if title != None:
        ax.set_title(title)
    
    
    #if out_img_file_path != None:
    #    fig.savefig(out_img_file_path)
    #print(f"Saved image file: {out_img_file_path}")
    
    return



"""

Create shapely.geometry.Polygon object of rectangle with bottom-left
 corner at given point, and given xsize and ysize
"""
def rect_from_coords_and_size(xcoord, 
                              ycoord, 
                              xsize, 
                              ysize
                              ):
    
    rect_corners = [(xcoord, ycoord), 
                    (xcoord, ycoord+ysize), 
                    (xcoord+xsize, ycoord+ysize), 
                    (xcoord+xsize, ycoord), 
                    ]
    return Polygon(rect_corners)
    

"""

"""
def make_points_frame(points, 
                      in_epsg, 
                      out_epsg=4326, 
                      ):
    
    gdf_dict = dict()
    xcoords = points.xcoords
    ycoords = points.ycoords
    gdf_dict["time"] = [str(x) for x in points.timestamps]
    
    
    gdf_dict["x"] = xcoords
    gdf_dict["y"] = ycoords
    gdf_dict["geometry"] = [Point(xy) for xy in zip(xcoords, ycoords)]
    
    
    gdf_datapoints = gpd.GeoDataFrame(gdf_dict)
    
    gdf_datapoints.crs = {'init' :f'epsg:{in_epsg}'}
    gdf_datapoints = gdf_datapoints.to_crs({'init': f'epsg:{out_epsg}'})
    
    
    return gdf_datapoints


def make_cells_frame(masked_grid, 
                     in_epsg, 
                     out_epsg=4326, 
                     ):
    
    
    id_list = []
    x_index_list = []
    y_index_list = []
    x_coord_list = []
    y_coord_list = []
    geometry_list = []
    
    mask = masked_grid.mask
    x_off = masked_grid.xoffset
    y_off = masked_grid.yoffset
    x_ext = masked_grid.xextent
    y_ext = masked_grid.yextent
    x_size = masked_grid.xsize
    y_size = masked_grid.ysize
    x_coords = list(range(x_off, x_off + (x_size * x_ext), x_size))
    y_coords = list(range(y_off, y_off + (y_size * y_ext), y_size))
    
    for yx_indices in product(range(y_ext), range(x_ext)):
        if mask[yx_indices[0]][yx_indices[1]]:
            continue
        
        x_coord = x_coords[yx_indices[1]]
        y_coord = y_coords[yx_indices[0]]
        cell = rect_from_coords_and_size(x_coord, y_coord, x_size, y_size)
        id_list.append(f"{yx_indices}")
        x_index_list.append(yx_indices[1])
        y_index_list.append(yx_indices[0])
        x_coord_list.append(x_coord)
        y_coord_list.append(y_coord)
        geometry_list.append(cell)
    
    
    gdf_dict = dict()
    gdf_dict= {"x_index":x_index_list, 
               "y_index":y_index_list, 
               "x_coord":x_coord_list, 
               "y_coord":y_coord_list, 
               "geometry":geometry_list,
               }
    
    gdf_cells = gpd.GeoDataFrame(gdf_dict)
    
    gdf_cells.crs = {'init' :f'epsg:{in_epsg}'}
    gdf_cells = gdf_cells.to_crs({'init': f'epsg:{out_epsg}'})
    
    return gdf_cells


"""
frame_to_json_with_id

Save a geodataframe as a geojson file, while also adding an "id" field
to each feature.

Note that this function uses the "Feature" constructor in the "geojson" 
module.

This function is necessary due to some unusual interplay between the
ipyleaflet module's Choropleth constructor and the format of geojson
files. The Choropleth constructor requires that its input geojson have
an "id" assigned to each feature. The geopandas function to_file was
previously used to save the GeoJSON file, with option "driver='GeoJSON'",
but I could not find a way to assign an id to each feature -- the closest
I could manage was assigning an id 'property' to each feature, but this
was not compatible with the Choropleth function of ipyleaflet.

If the Choropleth function is not used, possibly due to discovery of a
superior way to generate heatmaps on maps, then I recommend reverting
back to using the to_file function with driver='GeoJSON' for saving the
GeoJSON file, instead of using this function.
"""
def frame_to_json_with_id(frame, json_name):
    feat_coll_list = []
    for i in range(frame.shape[0]):
        this_geometry = frame["geometry"][i]
        this_properties = OrderedDict()
        for col_name in frame.columns:
            if col_name == "geometry":
                continue
            value = frame[col_name][i]
            if type(value)==numpy.int64:
                value = int(value)
            this_properties[col_name] = value
        this_id = str((this_properties["y_index"],this_properties["x_index"]))
        feat_coll_list.append(geojson.Feature(
                                geometry = this_geometry,
                                properties = this_properties,
                                id = this_id))
    feat_coll = geojson.FeatureCollection(feat_coll_list)
    
    with open(json_name,"w") as f:
        geojson.dump(feat_coll, f)
    return



def list_risk_model_properties(geojson_file):
    cell_results = gpd.read_file(geojson_file)
    ignore_properties = ["id", \
                         "x_index", \
                         "y_index", \
                         "x_coord", \
                         "y_coord", \
                         "geometry"]
    return [p for p in cell_results.columns if p not in ignore_properties]







"""
main:

If running this module as a script instead of importing its functions,
 this main function will perform () with a set of
 default parameters.
"""
def main():
    
    # Data directory
    datadir = "../../Data"
    
    # Geojson file
    #geojson_file_name = "Chicago_Areas.geojson"
    #geojson_file_name = "Chicago_South_Side_2790.geojson"
    #geojson_file_name = "Durham_27700.geojson"
    #geojson_file_name = "Police_Force_Areas_December_2016_Durham.geojson"
    geojson_file_name = "Police_Force_Areas_December_2016_Durham_fixed.geojson"
    
    geojson_full_path = os.path.join(datadir, geojson_file_name)
    
    region_polygon = gpd.read_file(geojson_full_path)
    
    union_polygon = region_polygon.unary_union
    
    print(type(region_polygon))
    print(type(union_polygon))
    
    print(region_polygon.crs)
    print(region_polygon.has_z)
    print(type(region_polygon.has_z))
    
    
    
    altered_region = region_polygon.to_crs({'init': 'epsg:27700'})
    
    plotPointsOnGrid([],[],altered_region.unary_union)
    
    
    
    



if __name__ == "__main__":
    main()
