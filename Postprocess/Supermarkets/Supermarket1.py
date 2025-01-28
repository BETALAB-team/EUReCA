import pandas as pd
import json


# Path to the JSON file
json_file_path = '4.0_20_MPEurospinIperrossettoViaSarpi.json'

# Read the JSON file
with open(json_file_path, "r") as file:
    json_data = json.load(file)

# Extract the components
columns = json_data["columns"]  # List of column names
data = json_data["data"]        # List of rows (as lists)
index = json_data["index"]      # Row indices

# Reconstruct the DataFrame
df_raw = pd.DataFrame(data, columns=columns, index=index)
columns_to_change = ["Solar Thermal Production [Wh]","PV production [Wh]",
                     "Given to Grid [Wh]",'Solar Thermal Surplus Without Tank [Wh]',
                     "Solar Thermal Surplus With Tank [Wh]"]

#%% create zero solar
import pandas as pd
import numpy as np


df_zero=df_raw.copy()
mask = ~df_zero['End Use'].isin(['supermarket', 'parking'])
for solar_column in columns_to_change:
    df_zero.loc[mask, solar_column] = df_zero.loc[mask, solar_column].apply(lambda arr: np.zeros_like(arr))

df_zero.loc[mask, "DHW Intake With Tank [Wh]"] = df_zero.loc[mask, "DHW Demand [Wh]"]
df_zero.loc[mask, "DHW Intake Without Tank [Wh]"] = df_zero.loc[mask, "DHW Demand [Wh]"]
#%% create 50% solar
from sklearn.utils import resample
import pandas as pd

random_seed = 17
df_half=df_raw.copy()
non_supermarket_parking = df_half[~df_half['End Use'].isin(['supermarket', 'parking'])]

sampled_indices = resample(non_supermarket_parking.index, 
                           replace=False, 
                           n_samples=len(non_supermarket_parking) // 2, 
                           random_state=random_seed)

for solar_column in columns_to_change:
    df_half.loc[sampled_indices, solar_column] = df_half.loc[sampled_indices, solar_column].apply(lambda arr: np.zeros_like(arr))
df_half.loc[mask, "DHW Intake With Tank [Wh]"] = df_half.loc[mask, "DHW Demand [Wh]"]
df_half.loc[mask, "DHW Intake Without Tank [Wh]"] = df_half.loc[mask, "DHW Demand [Wh]"]
#%%

import numpy as np
import pandas as pd

def filter_by_distance_to_supermarket(df):


    def haversine(lat1, lon1, lat2, lon2):
        R = 6371.0  # Earth's radius in kilometers
        # Convert degrees to radians
        lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
        # Differences in coordinates
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        # Haversine formula
        a = (dlat**2+dlon**2)**0.5
        return a

    # 1. Find the midpoint of rows with End Use as "supermarket"
    supermarket_rows = df[df['End Use'] == 'supermarket']
    if supermarket_rows.empty:
        raise ValueError("No rows with 'End Use' as 'supermarket' found in the dataframe.")
    
    midpoint_lat = supermarket_rows['latitude'].mean()
    midpoint_lon = supermarket_rows['longitude'].mean()
    # 2. Calculate distance to the midpoint for all rows
    df['Distance to Supermarket'] = df.apply(
        lambda row: haversine(row['latitude'], row['longitude'], midpoint_lat, midpoint_lon), axis=1
    )

    # 3. Filter rows where distance is less than the median
    median_distance = df['Distance to Supermarket'].median()
    print(median_distance)
    filtered_df = df[df['Distance to Supermarket'] < median_distance]
    filtered_df = filtered_df.drop(columns=['Distance to Supermarket'])
    return filtered_df


df_zero_small=filter_by_distance_to_supermarket(df_zero)
df_half_small=filter_by_distance_to_supermarket(df_half)
df_raw_small=filter_by_distance_to_supermarket(df_raw)
#%%
import os
name=json_file_path
output_dir = "Technologies"
os.makedirs(output_dir, exist_ok=True)
df_zero_small.to_json(os.path.join(output_dir,"NoSolar_Small_"+name), orient="split")
df_zero.to_json(os.path.join(output_dir,"NoSolar_Large_"+name), orient="split")
df_half_small.to_json(os.path.join(output_dir,"HalfSolar_Small_"+name), orient="split")
df_half.to_json(os.path.join(output_dir,"HalfSolar_Large_"+name), orient="split")
df_raw.to_json(os.path.join(output_dir,"AllSolar_Small_"+name), orient="split")
df_raw_small.to_json(os.path.join(output_dir,"AllSolar_Large_"+name), orient="split")