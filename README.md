# EUReCA 

![Insert caption here](https://research.dii.unipd.it/betalab/wp-content/uploads/sites/33/2021/03/EUReCA_logo_300x300.jpg)

The Energy Urban Resistance Capacitance Approach provides an efficient and reliable Urban Building Energy Modeling platform, entirely developed in Python, aiming at simulating and predicting cities and urban areas energy consumption. The tool exploits a bottom-up modeling methodology, creating simple and useful dynamic building energy models.

This research project has been developed within the [BETALAB](https://research.dii.unipd.it/betalab/) research group of the University of Padua, Italy.

## Python environment set up
The tool is distributed via the GitHub repository. As first step, you must create a new conda or venv environment. You can name it eureca.

`conda create -n eureca python=3.9`

and activate it:

``conda activate eureca``

Then install the following package in the created environment, the user can type:

``
pip install git+https://github.com/BETALAB-team/EUReCA
``

This will install the latest version. 

To install a specific version the user can type something like this:

``
pip install git+https://github.com/BETALAB-team/EUReCA@v1.0.0-beta
``

## Preparing and run a simulation
### Input files

The [eureca_ubem/Input](https://github.com/BETALAB-team/EUReCA/tree/main/eureca_ubem/Input) folder has some examples files to run the simulation. 
To simulate cities energy consumption in EUReCA, some input files must be prepared:
 - A `weather_data.epw` weather file. These files are available at the [EnergyPlus](https://www.energyplus.net/weather) website.
 - A `EnvelopeTypes.xlsx` spreadsheet. It includes the thermo-physic properties of building envelopes. An example is available in the `materials_and_construction_test.xlsx`
 - A `Schedules.xlsx` spreadsheet. It includes the operational schedules of occupancy, appliances, temperature, humidity setpoints, HVAC usage for different end-uses. Example in `Schedules.xlsx`.
 - A `Systems.xlsx` spreadsheet. It includes data about the systems types that the user can apply. Example in `Systems.xlsx`.
 - The `config.json` file, which defines the simulation parameters. Example in [config.json](https://github.com/BETALAB-team/EUReCA/blob/main/eureca_ubem/Input/config.json).
 - The `city.geojson` model. See the next section for further info on the alternatives.

### The JSON city model
Currently, EUReCA can handle geojson shapefiles. 

The required attributes are:
- GeoJSON: 
  ```
  "id": integer, 
  "Name": "name", 
  "End Use": "schedule_archetype_name", 
  "Envelope": "envelope_archetype_name", 
  "Height": float, "Nfloors": integer, 
  "Floors": float, 
  "Heating System": "heating_system_name", 
  "Cooling System": "cooling_system_name",
  "Lower End Use": "schedule_archetype_name for ground floor",   # Optional (altrenative to end use)
  "Upper End Use": "schedule_archetype_name for upper floors",   # Optional (altrenative to end use)
  "Solar technologies": ["PV", "ST", "PV; ST"],                  # Optional
  ```
  
The strings in the city model's attribute table (`End Use` and `Envelope`) must match the labels of the End Uses and Envelope types listed in the `Schedules.xlsx` and `EnvelopeTypes.xlsx`.

`Heating System` and `Cooling System` must match one of the names given to the systems in the `Systems.xlsx` file, or chosen from the following items:
List of available heating systems:
- IdealLoad
- CondensingBoiler
- TraditionalBoiler
- A-W HP Staffel, Centralized, Low Temp Radiator
- A-W HP Staffel, Centralized, High Temp Radiator
- A-W HP Staffel, Centralized, Fan coil
- A-W HP Staffel, Centralized, Radiant surface
- G-W HP Staffel, Centralized, Low Temp Radiator
- G-W HP Staffel, Centralized, High Temp Radiator
- G-W HP Staffel, Centralized, Fan coil
- G-W HP Staffel, Centralized, Radiant surface
- Traditional Gas Boiler, Centralized, Low Temp Radiator
- Traditional Gas Boiler, Single, Low Temp Radiator
- Traditional Gas Boiler, Centralized, High Temp Radiator
- Traditional Gas Boiler, Single, High Temp Radiator
- Traditional Gas Boiler, Centralized, Fan coil
- Traditional Gas Boiler, Single, Fan coil
- Traditional Gas Boiler, Centralized, Radiant surface
- Traditional Gas Boiler, Single, Radiant surface
- Condensing Gas Boiler, Centralized, Low Temp Radiator
- Condensing Gas Boiler, Single, Low Temp Radiator
- Condensing Gas Boiler, Centralized, High Temp Radiator
- Condensing Gas Boiler, Single, High Temp Radiator
- Condensing Gas Boiler, Centralized, Fan coil
- Condensing Gas Boiler, Single, Fan coil
- Condensing Gas Boiler, Centralized, Radiant surface
- Condensing Gas Boiler, Single, Radiant surface
- Oil Boiler, Centralized, High Temp Radiator
- Oil Boiler, Single, High Temp Radiator
- Coal Heater, Centralized, High Temp Radiator
- Coal Heater, Single, High Temp Radiator
- District Heating, Centralized, Low Temp Radiator
- District Heating, Centralized, High Temp Radiator
- District Heating, Centralized, Fan coil
- District Heating, Centralized, Radiant surface
- Stove
- A-W Heat Pump, Centralized, Low Temp Radiator
- A-W Heat Pump, Single, Low Temp Radiator
- A-W Heat Pump, Centralized, Fan coil
- A-W Heat Pump, Single, Fan coil
- A-W Heat Pump, Centralized, Radiant surface
- A-W Heat Pump, Single, Radiant surface
- Electric Heater

List of available cooling systems:
- IdealLoad
- SplitAirCooler
- ChillerAirtoWater
- SplitAirConditioner
- A-A split
- A-W chiller, Centralized, Fan coil
- A-W chiller, Centralized, Radiant surface
- A-W chiller, Single, Fan coil
- A-W chiller, Single, Radiant surface

Input folder provides some example for the city of Padua.

### Simulation

After the set up of all input files, you can run the simulation throughout a python file, as following:

```
import os
import time as tm


# CONFIG FILE LOADING
from eureca_building.config import load_config
load_config("path_to_your_config\\config.json")

from eureca_ubem.city import City

# SET INPUT FILES
weather_file = os.path.join(".","path_to_you_input","weather_file.epw")
schedules_file = os.path.join(".","path_to_you_input","Schedules.xlsx")
materials_file = os.path.join(".","path_to_you_input","materials_and_construction_test.xlsx")
city_model_file = os.path.join(".","path_to_you_input","citymodel.geojson")
systems_file = os.path.join(".","systems.xlsx")

# Creation of the City object and simulation
city_geojson = City(
    city_model=city_model_file,
    epw_weather_file=weather_file,
    end_uses_types_file=schedules_file,
    envelope_types_file=materials_file,
    systems_templates_file=systems_file,
)
city_geojson.loads_calculation()
city_geojson.simulate(print_single_building_results=True, output_type="csv")
```

### Output report
If `output_folder=os.path.join(".","your_output_folder")` is set, outputs are printed in the output folder.
Each file is a csv or a parquet with the main output variables of each building.

### How to cite EUReCA
In case you want to use EUReCA for your own research project, please cite the following paper: 

@article{\
PRATAVIERA2021544,\
title = {EUReCA: An open-source urban building energy modelling tool for the efficient evaluation of cities energy demand},\
journal = {Renewable Energy},\
volume = {173},\
pages = {544-560},\
year = {2021},\
issn = {0960-1481},\
doi = {https://doi.org/10.1016/j.renene.2021.03.144}, \
url = {https://www.sciencedirect.com/science/article/pii/S0960148121005085}, \
author = {Enrico Prataviera and Pierdonato Romano and Laura Carnieletto and Francesco Pirotti and Jacopo Vivian and Angelo Zarrella},\
keywords = {Urban building energy modelling, Lumped-capacitance thermal networks, Semantic georeferenced data, EUReCA, District simulation}\
}
