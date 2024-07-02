# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 13:15:33 2024

@author: khajmoh18975
"""
import pickle 
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import pandas as pd
import pvlib
import warnings
warnings.filterwarnings("ignore")


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



with open("std1.pkl", "rb") as file:
    surfaces = pickle.load(file)
    
latitude, longitude = 39.76, -105.22
tz = 'MST'

# Set times in the morning of the December solstice.
times = pd.date_range(
    '2020-12-20 0:00', '2020-12-21 0:00', freq='1T', tz=tz
)
# Create location object, and get solar position and clearsky irradiance data.
location = pvlib.location.Location(latitude, longitude, tz)
solar_position = location.get_solarposition(times)
clearsky = location.get_clearsky(times)
angle_sections=30
x_values = np.linspace(0, 360, num=360//angle_sections+1)
# for obj in surfaces:
#     # Convert each tuple within the attribute B to a list of NumPy arrays
#     obj._Surface__vertices = [np.array(tuple_) for tuple_ in obj._Surface__vertices]
list_of_centroids = [obj._centroid for obj in surfaces]
plg_kdtree=cKDTree(list_of_centroids)
maxis=50
max_number_of_neighborhoods = len(surfaces) if len(surfaces) < maxis else maxis
_, filtered_indices = plg_kdtree.query(list_of_centroids, k=max_number_of_neighborhoods)
indices_list = filtered_indices.tolist()
THE=0
for i in range(50):
    if (i%100==0):
        print(f"{i} from {len(surfaces)}")
    current_surface=surfaces[i]
    current_vertices=current_surface._Surface__vertices
    current_vertices= [np.array(tuple_) for tuple_ in current_vertices]
    current_surface_shader=current_surface
    normal_vec=current_surface._normal
    new_vertices = [arr-normal_vec*1 for arr in current_vertices]
    # current_surface_shader._Surface__vertices=new_vertices
    setattr(current_surface_shader,'_Surface__vertices',new_vertices)
    other_surfaces=surfaces.copy()
    other_surfaces=[other_surfaces[j] for j in indices_list[i]]
    other_surfaces.remove(current_surface)
    other_surfaces.append(current_surface_shader)
    Centroid=current_surface._centroid
    ABs=[Centroid]
    normal_vec=current_surface._normal
    a=normal_vec[0]
    b=normal_vec[1]
    c=normal_vec[2]
    if (a==0 and b==0):
        shader_starter=[np.array([0,0]),np.array([360,0])]
    else:
        shader_starter=[np.array([x,max(0,np.degrees(np.arctan(-b/(a*np.cos(np.deg2rad(x))+c*np.sin(np.deg2rad(x))))))]) for x in x_values]

    shader_starter=[np.array([0,0]),np.array([360,0])]
    
    for vertix in ABs:
    # for vertix in current_surface._vertices:
    # vertix= current_surface._Surface__vertices[2]
        # shaders=np.array([0,0])
        shaders=shader_starter



        


        horizon = pd.DataFrame({'x': x_values, 'y': np.zeros_like(x_values)})
        for j in range(len(other_surfaces)+1):

            if j==0:
                obj=other_surfaces[1]
                Surface_vertices=new_vertices
            else:
                obj=other_surfaces[j-1]
                Surface_vertices=obj._Surface__vertices
                
            Surface_vertices= [np.array(tuple_) for tuple_ in Surface_vertices]
            x_range=max(arr[0] for arr in Surface_vertices)-min(arr[0] for arr in Surface_vertices)
            y_range=max(arr[1] for arr in Surface_vertices)-min(arr[1] for arr in Surface_vertices)
            ranges=max(x_range,y_range)
            THE=max(THE,ranges)
            num_interpolated_points=int(ranges/40)+2
            Surface_vertices = interpolate_points_along_path(Surface_vertices, num_interpolated_points)
            
            q=Surface_vertices
            Surface_vertices=[arr-vertix for arr in Surface_vertices]
            
            p=Surface_vertices

            Polar_vertices=[np.array([((180+np.degrees(np.arctan2(arr[1],arr[0])))//angle_sections*angle_sections),
                                           (max(0,np.degrees(np.arctan2(arr[2],np.sqrt(arr[0]**2+arr[1]**2)))))]) for arr in Surface_vertices]
            
            

            if x_range>0:
                # num_interpolated_points=int(x_range)//angle_sections+3
                # Polar_vertices = interpolate_points_along_path(Polar_vertices, num_interpolated_points)
                for arr in Polar_vertices:
                    arr[0] =float(arr[0])//angle_sections*angle_sections
                Polar_vertices=np.vstack( Polar_vertices)
                e=Polar_vertices
                shaders=np.vstack((shaders,e))
        if (len(shaders)<3):
            shaders=[np.array([0,0]),np.array([0,0])]
        # print((shaders))
        # print("\n")
        df=pd.DataFrame(shaders,columns=['x', 'y'])
        df = df.groupby('x')['y'].max().reset_index()
        horizon = pd.merge(horizon, df, on='x', suffixes=('_A', '_B'), how='left')
        # Drop the 'y_B' column
        horizon['y_A'] = horizon['y_B'].fillna(horizon['y_A'])
        horizon.drop(columns=['y_B'], inplace=True)

        # Rename the 'y_A' column to 'y'
        horizon.rename(columns={'y_A': 'y'}, inplace=True)
        
        horizon_profile=horizon.set_index('x')['y']
        
        # if current_surface._surface_type=="Roof":
        #     plt.figure()
        #     plt.plot(horizon['x'], horizon['y'], color='red', linestyle='-')
            # plt.ylim([0,90])
        plt.figure()
        plt.plot(horizon['x'], horizon['y'], color='red', linestyle='-')
        plt.ylim([0,90])
        
        


        # Assign variable names for easier reading.
        surface_tilt = current_surface._height_round
        surface_azimuth =  current_surface._azimuth
        solar_azimuth = solar_position.azimuth
        solar_zenith = solar_position.apparent_zenith
        solar_elevation = solar_position.apparent_elevation
        dni = clearsky.dni
        ghi = clearsky.ghi
        dhi = clearsky.dhi
        
        horizon_elevation_data = np.interp(
            solar_azimuth, horizon_profile.index, horizon_profile
            )

        # Convert to Pandas Series for easier usage.
        horizon_elevation_data = pd.Series(horizon_elevation_data, times)

        # Adjust DNI based on data - note this is returned as numpy array
        dni_adjusted = np.where(solar_elevation > horizon_elevation_data, dni, 0)
        dni_adjusted=pd.Series(dni_adjusted,times)
        # Adjust GHI and set it to DHI for time-periods where 'dni_adjusted' is 0.
        # Note this is returned as numpy array
        ghi_adjusted = np.where(dni_adjusted == 0, dhi, ghi)
        ghi_adjusted=pd.Series(ghi_adjusted,times)
        # plt.figure(figsize=(10, 6))  # Adjust the figure size if needed
        # plt.plot(ghi, label='Series 1')
        # plt.plot(ghi_adjusted, label='Series 2')
        irrad_pre_adj = pvlib.irradiance.get_total_irradiance(
            surface_tilt, surface_azimuth, solar_zenith, solar_azimuth, dni, ghi, dhi
            )

        irrad_post_adj = pvlib.irradiance.get_total_irradiance(
            surface_tilt, surface_azimuth, solar_zenith, solar_azimuth, dni_adjusted,
            ghi_adjusted, dhi
            )

        plt.figure(figsize=(10, 6))  # Adjust the figure size if needed
        plt.plot(irrad_pre_adj.poa_global, label='Series 1')
        plt.plot(irrad_post_adj.poa_global, label='Series 2')

    