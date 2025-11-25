import os 
import geopandas as gpd
import pandas as pd
import numpy as np


Results_folder=os.path.join("..","..","eureca_ubem","Input","geojson_corr")

file_list = os.listdir(Results_folder)

# Initialize a list to store GeoDataFrames with the 'index' attribute
geo_dfs_with_index = []

for file_name in file_list:
    if file_name.endswith(".geojson"):
        gdf = gpd.read_file(os.path.join(Results_folder,file_name))
            
        if 'index' in gdf.columns:
            geo_dfs_with_index.append(gdf)
                

if geo_dfs_with_index:
    gdf = pd.concat(geo_dfs_with_index, ignore_index=True)
else:
    print("No Summary geojson with index attribute exists")
    
gdf["latitude"] = gdf.centroid.y
gdf["longitude"] = gdf.centroid.x
gdf = gdf[['index', 'End Use','latitude', 'longitude']]

def sum_same_sign_arrays(df, column_list, new_column_name, sign):

    def process_row(row, columns):
        # Initialize an array of zeros with the same shape as the arrays in the columns
        summed_array = np.zeros_like(row[columns[0]])
        
        for col in columns:
            # Replace any negative values in the array with zero
            if sign=="positive":
                cleaned_array = np.maximum(row[col], 0)
            if sign=="negative":
                cleaned_array = -np.minimum(row[col], 0)
            # Add the cleaned array to the summed_array
            summed_array += cleaned_array
        
        return summed_array
    
    # Apply the `process_row` function to each row in the dataframe and store the result in the new column
    df[new_column_name] = df.apply(lambda row: process_row(row, column_list), axis=1)

    return df
ts_per_hour=2
#%%

list_of_variables=['Solar Thermal Production [Wh]','Non-Renewable DHW [Wh]','Solar Surplus [Wh]','TZ DHW demand [W]',
                   'TZ DHW demand [W]','Solar Thermal Production [Wh]','Solar Thermal Production [Wh]',
                   'TZ sensible load [W]','TZ latent load [W]','TZ AHU pre heater load [W]','TZ AHU post heater load [W]',
                   'PV production [Wh]','Given to Grid [Wh]','Taken from the Gird [Wh]',
                   'Appliances electric consumption [Wh]','Refrigerator Heat Absorbed [Wh]','Refrigerator Heat Rejected [Wh]']
gdf[list_of_variables] = None
gdf["Emitter Temperature [C]"]=50
for i in range(len(gdf)):
    idx = gdf.iloc[i]['index']
    csv_file_path = os.path.join(Results_folder, "Results Bd "+idx+".csv")
    
    try:
        # Load the CSV file as a pandas DataFrame
        csv_df = pd.read_csv(csv_file_path,sep=";",low_memory=False)
        for value in list_of_variables:
            
            if value in csv_df.columns:
                add_value = csv_df[value][1:].astype(float)
                add_value=add_value.to_numpy()
                gdf.at[i, value] = add_value

            
    except FileNotFoundError:
        print(f"CSV file for index {idx} not found.")


gdf.loc[gdf['End Use'] == 'residential', 'Emitter Temperature [C]'] = 70
gdf.loc[gdf['End Use'] == 'supermarket', 'Emitter Temperature [C]'] = 40
sum_same_sign_arrays(gdf,['TZ sensible load [W]','TZ latent load [W]','TZ AHU pre heater load [W]','TZ AHU post heater load [W]'],
                     'Space Heating Load [W]',"positive")
gdf["Space Heating Load [Wh]"]=gdf["Space Heating Load [W]"]/ts_per_hour


sum_same_sign_arrays(gdf,['TZ sensible load [W]','TZ latent load [W]','TZ AHU pre heater load [W]','TZ AHU post heater load [W]'],
                     'Space Cooling Load [W]',"negative")
gdf["Space Cooling Load [Wh]"]=gdf["Space Cooling Load [W]"]/ts_per_hour

gdf.drop(['TZ sensible load [W]','TZ latent load [W]','TZ AHU pre heater load [W]','TZ AHU post heater load [W]','Space Cooling Load [W]',"Space Heating Load [W]"],axis=1,inplace=True)
gdf["Solar Thermal Surplus With Tank [Wh]"]=gdf["Solar Surplus [Wh]"]
gdf["DHW Demand [Wh]"]=gdf["TZ DHW demand [W]"]/ts_per_hour
gdf["Solar Thermal Surplus Without Tank [Wh]"]=gdf['Solar Thermal Production [Wh]']-gdf["TZ DHW demand [W]"]/ts_per_hour
sum_same_sign_arrays(gdf,["Solar Thermal Surplus Without Tank [Wh]"],"Solar Thermal Surplus Without Tank [Wh]","positive")
gdf["DHW Intake With Tank [Wh]"]=gdf['Non-Renewable DHW [Wh]']
gdf["DHW Intake Without Tank [Wh]"]=-gdf['Solar Thermal Production [Wh]']+gdf["TZ DHW demand [W]"]/ts_per_hour
sum_same_sign_arrays(gdf,["DHW Intake Without Tank [Wh]"],"DHW Intake Without Tank [Wh]","positive")
gdf["Refrigerator Electric Load [Wh]"]=-gdf['Refrigerator Heat Absorbed [Wh]']+gdf['Refrigerator Heat Rejected [Wh]']
gdf["Plug in Load [Wh]"]=gdf["Refrigerator Electric Load [Wh]"]+gdf["Appliances electric consumption [Wh]"]
gdf["Taken from the Grid [Wh]"]=gdf["Taken from the Gird [Wh]"]
gdf.drop(["Taken from the Gird [Wh]","TZ DHW demand [W]","Solar Surplus [Wh]",'Non-Renewable DHW [Wh]','Refrigerator Electric Load [Wh]', 'Appliances electric consumption [Wh]'],axis=1,inplace=True)



#%% 
df=gdf.copy()
df.to_json("output.json", orient="split")

#%%
for col in df.columns:
    print(f"{col}:")