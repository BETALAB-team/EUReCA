import os
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))  
project_root = os.path.abspath(os.path.join(script_dir, "..", ".."))

output_dir = os.path.join(project_root, "eureca_building", "output")
input_dir = script_dir  

print("\U0001F4C2 Script Directory:", script_dir)
print("\U0001F4C2 Project Root:", project_root)
print("\U0001F4C2 Input Directory:", input_dir)
print("\U0001F4C2 Output Directory:", output_dir)

if project_root not in sys.path:
    sys.path.append(project_root)
    print(f"✅ Added {project_root} a sys.path")

eureca_building_path = os.path.join(project_root, "eureca_building")
if eureca_building_path not in sys.path:
    sys.path.append(eureca_building_path)
    print(f"✅ Added {eureca_building_path} a sys.path")

eureca_ubem_path = os.path.join(project_root, "eureca_ubem")
if eureca_ubem_path not in sys.path:
    sys.path.append(eureca_ubem_path)
    print(f"✅ Added {eureca_ubem_path} a sys.path")

print("🔍 Python sys.path:")
for path in sys.path:
    print(f"   📂 {path}")

try:
    from eureca_building.config import load_config
    config_path = os.path.join(input_dir, "appsettings.json")
    print(f"🟢 Loading config file: {config_path}")
    load_config(config_path)
except ModuleNotFoundError as e:
    print(f"❌ ERRORE: {e}")
    raise

try:
    from eureca_ubem.city import City
    print("✅ Module City imported correctly!")
except ModuleNotFoundError as e:
    print(f"❌ ERRORE: {e}")
    raise

def run_simulation(env_file, end_use_file, climate_file, georeference_file):
    print("\n🟢 Files for simulation:")
    print(f"   📂 Env File: {env_file}")
    print(f"   📂 End Use File: {end_use_file}")
    print(f"   📂 Climate File: {climate_file}")
    print(f"   📂 Georeference File: {georeference_file}")

    
    input_files = {
        "Env File": env_file,
        "End Use File": end_use_file,
        "Climate File": climate_file,
        "Georeference File": georeference_file
    }

    for name, file in input_files.items():
        if not os.path.exists(file):
            print(f"⚠ ATTENTION: The file {name} doesn't exist ({file})!")
            return
        else:
            print(f"Found file: {file}")

    city_geojson = City(
        city_model=georeference_file,
        epw_weather_file=climate_file,
        end_uses_types_file=end_use_file,
        envelope_types_file=env_file,
        shading_calculation=True,
        output_folder=output_dir 
    )

    print(f"📂 Output Folder: {output_dir}")  

    try:
        print("🟢 Start simulation...")
        city_geojson.loads_calculation()
        city_geojson.simulate(print_single_building_results=True)
        print("✅ Simulation completed!")
    except Exception as e:
        print(f"❌ ERROR during simuation: {e}")
        raise


# import os
# import sys
# import time as tm

# script_dir = os.path.dirname(os.path.abspath(__file__))  

# project_root = os.path.abspath(os.path.join(script_dir, "..", ".."))

# output_dir = os.path.join(project_root, "eureca_building\\output")

# input_dir = script_dir  

# print("📂 Script Directory:", script_dir)
# print("📂 Project Root:", project_root)
# print("📂 Input Directory:", input_dir)
# print("📂 Output Directory:", output_dir)

# if project_root not in sys.path:
#     sys.path.append(project_root)
#     print(f"✅ Added {project_root} a sys.path")

# eureca_building_path = os.path.join(project_root, "eureca_building")
# if eureca_building_path not in sys.path:
#     sys.path.append(eureca_building_path)
#     print(f"✅ Added {eureca_building_path} a sys.path")

# eureca_ubem_path = os.path.join(project_root, "eureca_ubem")
# if eureca_ubem_path not in sys.path:
#     sys.path.append(eureca_ubem_path)
#     print(f"✅ Added {eureca_ubem_path} a sys.path")


# print("🔍 Python sys.path:")
# for path in sys.path:
#     print(f"   📂 {path}")

# try:
#     from eureca_building.config import load_config
#     config_path = os.path.join(input_dir, "config.json")
#     print(f"🟢 Loading config file: {config_path}")
#     load_config(config_path)
# except ModuleNotFoundError as e:
#     print(f"❌ ERRORE: {e}")
#     raise

# try:
#     from eureca_ubem.city import City
#     print("✅ Modulo City importato correttamente!")
# except ModuleNotFoundError as e:
#     print(f"❌ ERRORE: {e}")
#     raise

# weather_file = os.path.join(input_dir, "ITA_Venezia-Tessera.161050_IGDG.epw")
# schedules_file = os.path.join(input_dir, "Schedules_total.xlsx")
# materials_file = os.path.join(input_dir, "materials_and_construction_test.xlsx")
# city_model_file = os.path.join(input_dir, "PiovegoRestricted_with_holes_corr_coef.geojson")

# input_files = [config_path, weather_file, schedules_file, materials_file, city_model_file]
# for file in input_files:
#     if not os.path.exists(file):
#         print(f"⚠ ATTENZIONE: Il file {file} non esiste!")
#     else:
#         print(f"✅ File trovato: {file}")

# city_geojson = City(
#     city_model=city_model_file,
#     epw_weather_file=weather_file,
#     end_uses_types_file=schedules_file,
#     envelope_types_file=materials_file,
#     shading_calculation=True,
#     output_folder=output_dir 
# )

# print(f"📂 Output Folder: {output_dir}")  

# try:
#     print("🟢 Inizio simulazione...")
#     city_geojson.loads_calculation()
#     city_geojson.simulate(print_single_building_results=True)
#     print("✅ Simulazione completata con successo!")
# except Exception as e:
#     print(f"❌ ERRORE durante la simulazione: {e}")
#     raise
