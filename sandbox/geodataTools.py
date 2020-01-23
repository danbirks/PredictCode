# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 14:49:11 2019

@author: lawdfo
"""


import sys
import os
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from itertools import product
from descartes import PolygonPatch
from shapely.geometry import Point, Polygon
import geojson
import json
from collections import OrderedDict, defaultdict
import numpy
import ipyleaflet
from copy import deepcopy


"""
plotPointsOnGrid
"""
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
make_points_frame
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


"""
make_cells_frame
"""
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
        this_id = str((this_properties["y_index"],
                       this_properties["x_index"]))
        feat_coll_list.append(geojson.Feature(
                                geometry = this_geometry,
                                properties = this_properties,
                                id = this_id))
    feat_coll = geojson.FeatureCollection(feat_coll_list)
    
    with open(json_name,"w") as f:
        geojson.dump(feat_coll, f)
    return



def list_risk_model_properties(geojson_file_name=None,
                             geojson_file_contents=None):
    
    ignore_properties = ["id", \
                         "x_index", \
                         "y_index", \
                         "x_coord", \
                         "y_coord", \
                         "geometry"]
    
    cell_results = geojson_file_contents
    if cell_results == None:
        if geojson_file_name != None:
            with open(geojson_file_name) as eg:
                cell_results = json.load(eg)
            return [p for p in cell_results.columns \
                            if p not in ignore_properties]
        else:
            print("Error! Must supply either geojson_file_name or "+\
                  "geojson_file_contents to list_risk_model_properties")
            sys.exit(1)
    
    return [p for p in cell_results["features"][0]["properties"] \
                if p not in ignore_properties]


"""
sort_geojson_features
Input can be the name of the geojson file as a string,
 or a GeoDataFrame object itself from a geojson file.
"""
def sort_geojson_features(geojson, sort_property):
    cell_results = deepcopy(geojson)
    if type(geojson)==str:
        cell_results = gpd.read_file(geojson)
    return cell_results.sort_values(by=sort_property, ascending=False)
    

def top_geojson_features(geojson, sort_property, top_portion=0.01):
    sorted_frame = sort_geojson_features(geojson, sort_property)
    num_features = len(sorted_frame.index)
    top_num_features = int(top_portion * num_features)
    return sorted_frame[:top_num_features]





def json_dict_to_geoframe(jsondict):
    
    # Make geodataframe with columns for id, geometry, and all properties
    
    inverted_dict = defaultdict(list)
    for feat in jsondict['features']:
        try:
            inverted_dict['id'].append(feat['id'])
        except KeyError:
            print(feat)
            sys.exit(1)
        for prop, val in feat['properties'].items():
            inverted_dict[prop].append(val)
        geom_info = feat['geometry']
        shape_type_str = geom_info['type']
        coords = geom_info['coordinates']
        if shape_type_str=="Polygon" and len(coords)==1:
            coords = coords[0]
        to_eval = shape_type_str + "(" + str(coords) + ")"
        geom_obj = eval(to_eval)
        inverted_dict['geometry'].append(geom_obj)
    dataframe = pd.DataFrame(inverted_dict)
    geoframe = gpd.GeoDataFrame(dataframe)
    return geoframe









def marker_cluster_from_data(geojson_file):
    with open(geojson_file) as eg:
        datapoints = json.load(eg)
    marker_list = []
    for feat in datapoints['features']:
        loc = feat['geometry']['coordinates'][::-1]
        marker_list.append(ipyleaflet.CircleMarker(location=loc))
    marker_cluster = ipyleaflet.MarkerCluster(markers=marker_list)
    return marker_cluster


"""
normalise_geojson_features

geojson_data : If this is a string, then assume it is the name of a
                geojson file, and read in the data.
               If this is not a string, assume it is the data read from
                a geojson file (which would be a dict object).
norm_property : The feature property to normalise
norm_type : If "max" then scale the values via linear transformation so that
             the largest is 1 and smallest is 0.
             This is achieved by:
                 1st) Subtract smallest value from all values.
                 2nd) Divide all values by (now) largest value
            If "sum" then scale the values via linear transformation so that
             the sum of all values totals 1.
             This is achieved by dividing all values by the sum of all
              the original values.
"""
def normalise_geojson_features(geojson_data,
                               norm_property, 
                               norm_type="max"
                               ):
    
    norm_type = norm_type.lower()
    recognised_norm_types = ["max","sum"]
    if norm_type not in recognised_norm_types:
        print("Error! Unexpected norm_type for normalise_geojson_features!")
        print(f"Given norm_type: {norm_type}")
        print(f"Recognised options: {recognised_norm_types}")
        sys.exit(1)
    
    # If the data is a string, then assume it's a file name to read.
    #  Otherwise, assume it's the data (dict) itself
    if type(geojson_data) == str:
        with open(geojson_data) as gf:
            geojson_data = json.load(gf)
    
    
    value_list = [x['properties'][norm_property] for x in geojson_data['features']]
    
    if norm_type == "max":
        min_val = min(value_list)
        max_val = max(value_list)
        norm_data = deepcopy(geojson_data)
        for feat in norm_data['features']:
            feat['properties'][norm_property] -= min_val
            feat['properties'][norm_property] /= max_val
    elif norm_type == "sum":
        val_sum = sum(value_list)
        norm_data = deepcopy(geojson_data)
        for feat in norm_data['features']:
            feat['properties'][norm_property] /= val_sum
    
    return norm_data
    


"""
combine_geojson_features
"""
def combine_geojson_features(geojson_data_list, 
                             combine_property_list, 
                             multiplier_list = 1, 
                             new_property_name=None, 
                             new_file_name=None, 
                             normalisation=None, 
                             ):
    
    combine_property_list = [p.strip() for p in \
                                 combine_property_list.split(",")]
    
    # Determine how long all the input lists should be.
    # If multiple actual lists exist, check they're all the same length
    lists_to_check = (geojson_data_list, 
                      combine_property_list, 
                      multiplier_list)
    real_lists = [x for x in lists_to_check if type(x)==list and len(x)>1]
    lists_len = 1  # Default length for all 3 input lists
    if len(real_lists)>0:
        lens_seen = set([len(x) for x in real_lists])
        if len(lens_seen) != 1:
            print("Error! Not all lists are the same length!")
            list_names = ["geojson_data_list", 
                          "combine_property_list", 
                          "multiplier_list"]
            for ln, ltc in zip(list_names, lists_to_check):
                data_len = "not a list"
                if type(ltc)==list:
                    data_len = len(ltc)
                print(f"{ln} : {data_len}")
            sys.exit(1)
        lists_len = lens_seen.pop()
    
    
    # If an input list is not a list, then we create a list that
    #  just repeats it the right number of times
    if type(geojson_data_list) != list:
        geojson_data_list = [geojson_data_list]*lists_len
    if type(combine_property_list) != list:
        combine_property_list = [combine_property_list]*lists_len
    if type(multiplier_list) != list:
        multiplier_list = [multiplier_list]*lists_len
    
    
    
    # If no specified new property name, then set it to be
    #  all the original names, deduplicated, and concatenated by "_"
    if new_property_name == None:
        prop_set = set([p.lower() for p in combine_property_list])
        new_property_name = "_".join(sorted(prop_set))
    
    
    
    # If the data is a string, then assume it's a file name to read.
    #  Otherwise, assume it's the data (dict) itself
    for i, gd in enumerate(geojson_data_list):
        if type(gd) == str:
            gd = os.path.normpath(gd)
            gd = os.path.expanduser(gd)
            gd = os.path.expandvars(gd)
            if new_file_name == None:
                new_dir = os.path.dirname(gd)
                new_base = f"results_{new_property_name}.geojson"
                new_file_name = os.path.join(new_dir, new_base)
                new_file_name = new_file_name.replace("\\","/")
            with open(gd) as gf:
                geojson_data_list[i] = json.load(gf)
    
    
    
    # Normalise if specified
    if type(normalisation)==str:
        normalisation = normalisation.lower()
    recognised_norm_types = ["max","sum","none",None,False]
    if normalisation not in recognised_norm_types:
        print("Error! Unexpected norm_type for normalise_geojson_features!")
        print(f"Given norm_type: {normalisation}")
        print(f"Recognised options: {recognised_norm_types}")
        sys.exit(1)
    if normalisation in ["max","sum"]:
        for i, gd in enumerate(geojson_data_list):
            geojson_data_list[i] = normalise_geojson_features(gd,
                                                    combine_property_list[i], 
                                                    norm_type=normalisation
                                                    )
    
    
    
    save_properties = ['x_index', 'y_index', 'x_coord', 'y_coord']
    combo_data = deepcopy(geojson_data_list[0])
    feat_id_map = dict()
    for i, feat in enumerate(combo_data['features']):
        feat_id_map[feat['id']] = i
        # Clear all properties besides basic ID ones
        props_to_del = [p for p in feat['properties'] \
                                    if p not in save_properties]
        for prop in props_to_del:
            del combo_data['features'][i]['properties'][prop]
        feat['properties'][new_property_name] = 0
    
    for i, gd in enumerate(geojson_data_list):
        for feat in gd['features']:
            add_value = feat['properties'][combine_property_list[i]] * \
                                                multiplier_list[i]
            feat_id = feat['id']
            combo_feat_index = feat_id_map[feat_id]
            combo_feat = combo_data['features'][combo_feat_index]
            combo_feat['properties'][new_property_name] += add_value
    
    
    with open(new_file_name, "w") as f:
        f.write(str(combo_data))
    
    return combo_data, new_property_name, new_file_name



def get_superclass(c):
    return [b.__name__ for b in c.__class__.__bases__]


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
