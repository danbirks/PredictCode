# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 14:49:11 2019

@author: lawdfo
"""



import sys
import os
import geopandas as gpd
import matplotlib.pyplot as plt
from descartes import PolygonPatch
#from open_cp.plot import patches_from_grid
#from matplotlib.collections import PatchCollection








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
