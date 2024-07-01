
"""
Shading Horizon Calculation

This piece of script defines a function that can modify the shading_horzion
attribute in the surface objects. 

shading horizon is calculated using the relative distance in polar coordinates 
between each two vertices. the shading coefficient then can be calculated
by comparing the solar position in a certain moment and the height angle that 
the shading horzion covers. There are two main functions in this script:
    
create_shading_horizon:
    input: list of all surfaces (surface Class) that are considered in shading calculation
    function: modification of the attribute _shading_horizon in the surface objects
    based on the centroid of the surface
    
    this function is used in under city.py to add the horizon list to each surface

calculate_shading_coefficient:

---
Created by: Mohamad
Betalab - DII, University of Padua
---
"""


import numpy as np

from scipy.spatial import cKDTree
import pandas as pd
from eureca_building.config import CONFIG
import warnings
warnings.filterwarnings("ignore")


"""
---auxillary functions---
"""
def interpolate_points_between_arrays(arr1, arr2, num_points):
    interpolated_points = []
    for i in range(num_points):
        t = i / (num_points - 1)  # Interpolation parameter
        interpolated_point = (1 - t) * arr1 + t * arr2
        interpolated_points.append(interpolated_point)
    return np.array(interpolated_points)
def interpolate_points_along_path(points_list, num_points):
    interpolated_points = []
    for i in range(len(points_list) - 1):
        interpolated_points.extend(
            interpolate_points_between_arrays(points_list[i], points_list[i+1], num_points)
        )
    return np.array(interpolated_points)

"""
---ccreate shading horizon---
"""
def create_shading_horizon(surfaces,surface_to_calculate,indices):
    #setup the angles at which the shading horizon will be calculated
   angle_sections=CONFIG.azimuth_subdivisions
   angle_values=np.linspace(0, 360, num=angle_sections+1)
   # # find the centroid of each surface
   # # list_of_centroids = [obj._centroid for obj in surfaces]
   # # # find the closest surfaces to each centroid
   # # plg_kdtree=cKDTree(list_of_centroids,kdtree)
   # max_considered_shading_surface=50
   # max_number_of_neighborhoods = len(surfaces) if len(surfaces) < max_considered_shading_surface else max_considered_shading_surface
   # # _, filtered_indices = kdtree.query(centroids, k=max_number_of_neighborhoods)    
   # indices_list = indices.tolist()
   # i=0
   # #do the calculation for each centroid point:
   # for surface_index in range(len(surfaces)):
   #set up the current surface that is going to be analyzed
   vertices_to_calculate=surface_to_calculate._Surface__vertices
   vertices_to_calculate=[np.array(tuple_) for tuple_ in vertices_to_calculate] #change to array for convention
   #Create a surface similar to the current surface but a little bit
   #behind it so it can show that sunlight does not get to behind the surface
   normal_vector=surface_to_calculate._normal
   self_shade_surface_vertices=[array - normal_vector*0.1 for array in vertices_to_calculate] 
   self_shader_surface=surface_to_calculate
   setattr(self_shader_surface,'_Surface_vertices',self_shade_surface_vertices)
   #create shader surfaces
   shader_surfaces=surfaces.copy()
   # filter the closest ones
   shader_surfaces=[shader_surfaces[tree_index] for tree_index in indices]
   #remove the surface in study from the batch 
   shader_surfaces.remove(surface_to_calculate)
   #add self shader
   # shader_surfaces.append(self_shader_surface)
        
   #initialize calculation of horizon
   point_to_calculate=surface_to_calculate._centroid
   shading_horizon_start=[np.array([0,0]),np.array([360,0])]
   clear_horizon=pd.DataFrame({'azimuth':angle_values,'tilt':np.zeros_like(angle_values)})
        
        
   shading_horizon=shading_horizon_start
   for shading_surface_index in range(len(shader_surfaces)):
      #create shading surface vertices as numpy arrays
      shading_surface_vertices=shader_surfaces[shading_surface_index]._Surface__vertices
      shading_surface_vertices=[np.array(tuple_) for tuple_ in shading_surface_vertices]
      #create points alongside the edges to increase the precision
      num_interpolated_points=2 
      shading_surface_vertices_interpolated=interpolate_points_along_path(shading_surface_vertices, num_interpolated_points)
      #calculate the vertices in polar coordinate (height and azimuth) relative to the 
      #point being studied
      shading_surface_vertices_relative=[array-point_to_calculate for array in shading_surface_vertices_interpolated if array[2]>=0]
      shading_surface_vertices_polar=[np.array([((180+np.degrees(np.arctan2(arr[1],arr[0])))//angle_sections*angle_sections),
                                    (max(0,np.degrees(np.arctan2(arr[2],np.sqrt(arr[0]**2+arr[1]**2)))))]) for arr in shading_surface_vertices_relative]
      for array in shading_surface_vertices_polar:
          array[0] =float(array[0])//angle_sections*angle_sections
          shading_vertices_polar=np.vstack(shading_surface_vertices_polar)
          shading_horizon=np.vstack((shading_horizon,shading_vertices_polar))
    #process to find maximum height in each angle
   Horizon_dataframe=pd.DataFrame(shading_horizon,columns=['azimuth','tilt'])
   Horizon_dataframe_processed=Horizon_dataframe.groupby('azimuth')['tilt'].max().reset_index()
   Horizon_dataframe_merged=pd.merge(clear_horizon,Horizon_dataframe_processed,
                                        on='azimuth',suffixes=('_A', '_B'), how='left')
        
   Horizon_dataframe_merged['tilt_A'] = Horizon_dataframe_merged['tilt_B'].fillna(Horizon_dataframe_merged['tilt_A'])    
   Horizon_profile_list=Horizon_dataframe_merged['tilt_A'].tolist()
   surface_to_calculate._shading_horizon=Horizon_profile_list
   pass

"""
---calculate shading coefficient---
"""
def calculate_shading_coefficient(hourly_data,surface):
    #initialize the calculation setup
    angle_sections=CONFIG.azimuth_subdivisions
    angle_values=np.linspace(0, 360, num=angle_sections+1)
    horizon=pd.Series(surface._shading_horizon, index=angle_values)
    #interpolate the horizon data to match current solar azimuth
    horizon_elevation_data = np.interp(
    hourly_data['solar_position_azimuth'], 
    horizon.index, 
    horizon)
    #check if the centroid sees the light
    is_lit = (hourly_data['solar_position_apparent_elevation'] > horizon_elevation_data).astype(int)
    return is_lit



    