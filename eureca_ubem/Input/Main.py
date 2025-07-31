def main():
    """
    Entry point for the Eureca Building Energy Simulation.

    This function initializes the configuration, loads necessary input files,
    creates the City object, and runs the dynamic simulation.

    Workflow:
    ---------
    1. Load simulation configuration from a JSON file.
    2. Set up file paths for weather, schedules, materials, systems, and city model.
    3. Create a City object using the provided inputs.
    4. Run the simulation (dynamic or quasi-steady-state).

    Notes:
    ------
    - Ensure that all paths are correctly set relative to the script location.
    - The config module must be loaded before importing other modules that rely on it.
    - Outputs are saved in the specified output folder.
    """
    import os
    from eureca_building import config
    # === 1. Load Configuration ===
    config_path = os.path.join(".", "Example_District_Config.json")
    config.load_config(config_path)  # Initializes global CONFIG object

    # === 2. Define Input File Paths ===
    weather_file    = os.path.join(".", "ITA_Venezia-Tessera.161050_IGDG.epw")  # EPW weather file
    schedules_file  = os.path.join(".", "Schedules_total.xlsx")                 # Occupancy/end-use schedules
    materials_file  = os.path.join(".", "Materials.xlsx")                       # Construction materials
    city_model_file = os.path.join(".", "Example_District.geojson")             # Building footprints (GeoJSON)
    systems_file    = os.path.join(".", "systems.xlsx")                         # HVAC system templates
    output_folder   = os.path.join(".", "Output_folder")                        # Where simulation results will go

    # === 3. Create City Object ===
    # Import here to ensure config has been initialized first
    from eureca_ubem.city import City

    city_geojson = City(
        city_model=city_model_file,
        epw_weather_file=weather_file,
        end_uses_types_file=schedules_file,
        envelope_types_file=materials_file,
        systems_templates_file=systems_file,
        shading_calculation=True,    # Enable solar shading (slower)
        building_model="2C",         # Use 2-node thermal model (VDI 6007)
        output_folder=output_folder
    )

    # === 4. Run Simulation ===
    city_geojson.simulate(output_type="csv")  # Dynamic simulation results in CSV format
    # city_geojson.simulate_quasi_steady_state()  # Uncomment for quasi-steady-state alternative


if __name__ == "__main__":
    main()