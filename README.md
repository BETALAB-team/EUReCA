# EUReCA

The **E**nergy **U**rban **Re**sistance **C**apacitance **A**pproach provides an efficient and reliable Urban Building Energy Modeling platform, entirely developed in Python, aiming at simulating and predicting cities and urban areas energy consumption. The tool exploits a bottom-up modeling methodology, creating simple and useful dynamic building energy models.

This research project has been developed within the [BETALAB](https://research.dii.unipd.it/betalab/) research group of the University of Padua

## Python environment set up
The tool is distributed via the GitHub repository. It can be freely cloned with git, `git clone https://github.com/BETALAB-team/EUReCA.git`, or clicking on the **code** button and downloading the zip file.

An eureca.yml file is also included in the repository. It provides the python packages needed to run the EUReCA simulation. 

The virtual environment can be easily set up using the [Anaconda](https://www.anaconda.com/products/individual) package manager. With the command line: 
```
conda -env create -f EUReCA_PATH\eureca.yml
```
Or using the Anaconda navigator application, under `Environments -> Import -> Import environment from path` and searching for  `EUReCA_PATH\eureca.yml`. The environment will be set up for you by the conda package manager. 

To activate the RC environment with the command :
```
conda activate eureca
```

## Preparing and run a simulation
### Input files

To simulate cities energy consumption in EUReCA, some input files must be prepared:
- A `weather_data.epw` weather file. These files are available at the [EnergyPlus](https://www.energyplus.net/weather) website.
- A `Envelopes.xlsx` spreadsheet. It includes the thermophysical properties of building envelopes. An example is available in the `Input.Envelopes.xlsx`
- A `Schedules.xlsx` spreadsheet. It includes the operational schedules of occupancy, appliances, temperature, humidity setpoints, HVAC usage for different end-uses. There are two possible ways to set up the end-uses usage: using the daily mode (example in `Input\ScheduleSemp.xlsx`) and the Yearly mode (example in `Input\ScheduleComp.xlsx`).
- A `PlantList.xlsx` spreadsheet. The file includes the input data of many plants model, for different sizes.
- The `SimInput` file, which defines the simulation parameters. `SimInput.txt` or `SimInput.xlsx` can be used as a reference.
- The `city.json` model. See the next section for further info on the alternatives.

### The JSON city model
Currently, EUReCA can handle two typologies of JSON city models. The recommended methodology consists of importing buildings' geometries via semantic [CityJSON](https://www.cityjson.org/) files, but also 2D shapefiles, encoded in GeoJSON format, can be utilized to build up the city.

The required attributes are:
- CityJSON: 
  ```
  "Use": "schedule_archetype_name", 
  "Age": "envelope_archetype_name", 
  "H_Plant": "heating_plant_name", 
  "C_Plant": "cooling_plant_name"
  ```
- GeoJSON: 
  ```
  "id": integer, 
  "Name": "name", 
  "Use": "schedule_archetype_name", 
  "Age": "envelope_archetype_name", 
  "Height": float, "Nfloors": integer, 
  "ExtWallCoeff": float, 
  "VolCoeff": float, 
  "C_Plant": "cooling_plant_name", 
  "H_Plant": "heating_plant_name"
  ```

Input folder provides some example for the city of Padua.

### Simulation

After the set up of all input files, you can run the Main.py file:

```
''' IMPORTING MODULES '''

import os
import numpy as np
from RC_classes.Simulation import Sim


# Creation of the Sim object
city = Sim()

# Loading the input data: just uncomment one of the following lines depending on the simulation input file you filled
city.set_input_from_text_file(os.path.join('.','Input','SimInput'))
# city.set_input_from_excel_file(os.path.join('.','Input','SimInput.xlsx'))

# Loading weather data, envelopes and schedules
city.preprocessing()

# Creation of the district (geometrical processing)
city.city_creation()

# Evaluating Urban shadings between buildings
city.urban_shading()

# Calculation buildings parameters
city.buildings_params_and_loads()

# Design power of buildings and plants creation
city.plants_design_and_creation()

# Annual simulation
city.simulation()

# Output processing
city.output()
```

### Output report
In case you ran the command `city.output()` the folder OutputReport will be created as final step.

The report consists of several `output_variable.csv` files, including many outputs, as buildings' temperature, humidity, and consumption.

A `warning.txt` file is printed, as well. 
