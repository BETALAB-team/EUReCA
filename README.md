# EUReCA

![Insert caption here](https://research.dii.unipd.it/betalab/wp-content/uploads/sites/33/2021/03/EUReCA_logo_300x300.jpg)

The **E**nergy **U**rban **Re**sistance **C**apacitance **A**pproach provides an efficient and reliable Urban Building Energy Modeling platform, entirely developed in Python, aiming at simulating and predicting cities and urban areas energy consumption. The tool exploits a bottom-up modeling methodology, creating simple and useful dynamic building energy models.

This research project has been developed within the [BETALAB](https://research.dii.unipd.it/betalab/) research group of the University of Padua.

## Python environment set up
The tool is distributed via the GitHub repository. As first step, you must create a new conda or venv environment. You can name it eureca.

> conda create -n eureca python 3.9

and activate it:

> conda activate eureca

Then clone the following package in a separate folder:

> git clone https://github.com/BETALAB-team/EUReCA.git

and install it in the same environement:

> pip install -e *your_path_to_the_folder/eureca-ubem*

## Preparing and run a simulation
### Input files

The eureca_ubem/Input folder contain some examples file to run the simulation. 
To simulate cities energy consumption in EUReCA, some input files must be prepared:
- A `weather_data.epw` weather file. These files are available at the [EnergyPlus](https://www.energyplus.net/weather) website.
- A `EnvelopeTypes.xlsx` spreadsheet. It includes the thermo-physic properties of building envelopes. An example is available in the `materials_and_construction_test.xlsx`
- A `Schedules.xlsx` spreadsheet. It includes the operational schedules of occupancy, appliances, temperature, humidity setpoints, HVAC usage for different end-uses. Example in `Schedules.xlsx`.
- The `config.json` file, which defines the simulation parameters. Example in `config.json`.
- The `city.json` model. See the next section for further info on the alternatives.

### The JSON city model
Currently, EUReCA can handle two typologies of JSON city models. The recommended methodology consists of importing buildings' geometries via semantic [CityJSON](https://www.cityjson.org/) files, but also 2D shapefiles, encoded in GeoJSON format, can be utilized to build up the city.

The required attributes are:
- CityJSON: 
  ```
  "End Use": "schedule_type_name", 
  "Envelope": "envelope_type_name", 
  "Heating System": "heating_system_name", 
  "Cooling System": "cooling_system_name"
  ```
- GeoJSON: 
  ```
  "id": integer, 
  "Name": "name", 
  "End Use": "schedule_archetype_name", 
  "Envelope": "envelope_archetype_name", 
  "Height": float, "Nfloors": integer, 
  "Floors": float, 
  "Heating System": "heating_system_name", 
  "Cooling System": "cooling_system_name"
  ```

Input folder provides some example for the city of Padua.

List of available heating systems:
- CondensingBoiler
- TraditionalBoiler

List of available cooling systems:
- SplitAirCooler
- SplitAirConditioner
- ChillerAirtoWater

### Simulation

After the set up of all input files, you can run the Main.py file:

```
import os
import time as tm


# CONFIG FILE LOADING
from eureca_building.config import load_config
load_config("config.json")

from eureca_ubem.city import City

# SET INPUT FILES
weather_file = os.path.join(".","ITA_Venezia-Tessera.161050_IGDG.epw")
schedules_file = os.path.join(".","Schedules.xlsx")
materials_file = os.path.join(".","materials_and_construction_test.xlsx")
city_model_file = os.path.join(".","PiovegoRestricted_with_holes.geojson")

# Creation of the City object and simulation
city_geojson = City(
    city_model=city_model_file,
    epw_weather_file=weather_file,
    end_uses_types_file=schedules_file,
    envelope_types_file=materials_file,
    shading_calculation=True,
    output_folder=os.path.join(".","geojson")
)
city_geojson.loads_calculation()
city_geojson.simulate()
```

### Output report
If `output_folder=os.path.join(".","name_output")` is set an output file for each building is created in the output folder.
Each file is a csv with the main output variables of each building.

### How to cite EUReCA
In case you want to use EUReCA for your own research project, please cite this paper: 

@article{PRATAVIERA2021544,
title = {EUReCA: An open-source urban building energy modelling tool for the efficient evaluation of cities energy demand},
journal = {Renewable Energy},
volume = {173},
pages = {544-560},
year = {2021},
issn = {0960-1481},
doi = {https://doi.org/10.1016/j.renene.2021.03.144},
url = {https://www.sciencedirect.com/science/article/pii/S0960148121005085},
author = {Enrico Prataviera and Pierdonato Romano and Laura Carnieletto and Francesco Pirotti and Jacopo Vivian and Angelo Zarrella},
keywords = {Urban building energy modelling, Lumped-capacitance thermal networks, Semantic georeferenced data, EUReCA, District simulation}
}
