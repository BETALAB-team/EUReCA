import geopandas as gpd
import tkinter as tk
from tkinter import ttk



ENVELOPE_LEVELS = [
    "none",
    "shallow",
    "medium",
    "deep"
]

HEATING_OPTIONS = [
    "boiler",
    "dhn",
    "hp_le",
    "hp_me",
    "hp_he"
]

FUEL_OPTIONS = [
    "gas",
    "bio"
]

PV_TYPES = [
    "none",
    "A",
    "B"
]

def building_editor(geojson_path):

    gdf = load_buildings(geojson_path)
    baseline = gdf.copy()
    root, frame = initialize_gui()

    rows = create_building_rows(frame, gdf)

    attach_callbacks(rows)

    result = run_interface(root, rows, gdf)
    interventions = build_intervention_dict(baseline, gdf)


    return result, interventions


def build_intervention_dict(baseline, gdf):

    interventions = {}

    for idx in gdf.index:

        base = baseline.loc[idx]
        new  = gdf.loc[idx]

        changes = {}

        # envelope
        if base["EEdepth"] != new["EEdepth"]:
            changes["envelope"] = (base["EEdepth"], new["EEdepth"])

        # heating
        if base["SHSource"] != new["SHSource"]:
            changes["sh_source"] = (base["SHSource"], new["SHSource"])

        # dhw
        if base["DHWsource"] != new["DHWsource"]:
            changes["dhw_source"] = (base["DHWsource"], new["DHWsource"])

        # PV type
        if base["PVType"] != new["PVType"]:
            changes["PVtype"] = (base["PVType"], new["PVType"])

        # PV %
        if float(base["PVpercentage"]) != float(new["PVpercentage"]):
            changes["PVpercentage"] = (
                float(base["PVpercentage"]),
                float(new["PVpercentage"])
            )

        interventions[idx] = changes

    return interventions

def load_buildings(path):
    """
    Load the building GeoJSON and perform minimal validation.
    
    Parameters
    ----------
    path : str
        Path to the GeoJSON file.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the buildings.
    """

    gdf = gpd.read_file(path)

    required_columns = [
        "Name",
        "EEdepth",
        "SHSource",
        "DHWsource",
        "PVType",
        "PVpercentage"
    ]

    missing = [c for c in required_columns if c not in gdf.columns]

    if missing:
        raise ValueError(f"Missing required columns in GeoJSON: {missing}")

    return gdf



def initialize_gui():
    root = tk.Tk()
    root.title("Building configuration")
    root.geometry("950x600")

    # Main container
    main = ttk.Frame(root)
    main.pack(fill="both", expand=True)

    # Canvas for scrolling
    canvas = tk.Canvas(main)
    canvas.pack(side="left", fill="both", expand=True)

    scrollbar = ttk.Scrollbar(main, orient="vertical", command=canvas.yview)
    scrollbar.pack(side="right", fill="y")

    canvas.configure(yscrollcommand=scrollbar.set)

    # Frame inside canvas
    inner_frame = ttk.Frame(canvas)
    canvas.create_window((0, 0), window=inner_frame, anchor="nw")

    def on_configure(event):
        canvas.configure(scrollregion=canvas.bbox("all"))

    inner_frame.bind("<Configure>", on_configure)

    # Headers
    headers = [
        "Building",
        "Envelope",
        "Heating",
        "DHW",
        "Fuel",
        "PV Type",
        "PV %",
        ""
    ]

    for c, h in enumerate(headers):
        ttk.Label(
            inner_frame,
            text=h,
            font=("Arial", 10, "bold")
        ).grid(row=0, column=c, padx=6, pady=6)

    return root, inner_frame




def create_building_rows(frame, gdf):

    rows = {}

    for r, (idx, b) in enumerate(gdf.iterrows(), start=1):

        ttk.Label(frame, text=b["Name"]).grid(row=r, column=0, padx=4, pady=2)

        # ---------- Envelope ----------
        
        baseline_env = b["EEdepth"]
        
        envelope_order = ["none", "shallow", "medium", "deep"]
        
        idx_env = envelope_order.index(baseline_env)
        
        allowed_env = envelope_order[idx_env:]
        
        
        env_var = tk.StringVar()
        
        env_box = ttk.Combobox(
            frame,
            textvariable=env_var,
            values=allowed_env,
            state="readonly",
            width=10
        )
        
        env_box.grid(row=r, column=1, padx=4)
        
        env_box.set(baseline_env)

        # -------- Heating --------
        baseline_heat = translate_source(b["SHSource"])
        
        hp_order = ["hp_le", "hp_me", "hp_he"]
        
        if baseline_heat in hp_order:
        
            idx_hp = hp_order.index(baseline_heat)
        
            allowed_hp = hp_order[idx_hp:]
        
            allowed_heat = ["dhn"] + allowed_hp
        
        else:
        
            allowed_heat = HEATING_OPTIONS
        
        
        heat_var = tk.StringVar()
        
        heat_box = ttk.Combobox(
            frame,
            textvariable=heat_var,
            values=allowed_heat,
            state="readonly",
            width=10
        )
        
        heat_box.grid(row=r, column=2, padx=4)
        
        heat_box.set(baseline_heat)


        # -------- DHW --------
        baseline_dhw = translate_source(b["DHWsource"])
        
        hp_order = ["hp_le", "hp_me", "hp_he"]
        
        if baseline_dhw in hp_order:
        
            idx_hp = hp_order.index(baseline_dhw)
        
            allowed_hp = hp_order[idx_hp:]
        
            allowed_dhw = ["dhn"] + allowed_hp
        
        else:
        
            allowed_dhw = HEATING_OPTIONS
        
        
        dhw_var = tk.StringVar()
        
        dhw_box = ttk.Combobox(
            frame,
            textvariable=dhw_var,
            values=allowed_dhw,
            state="readonly",
            width=10
        )
        
        dhw_box.grid(row=r, column=3, padx=4)
        
        dhw_box.set(baseline_dhw)


        # -------- Fuel --------
        fuel_default = extract_fuel(b["SHSource"]) or extract_fuel(b["DHWsource"])
        if fuel_default not in FUEL_OPTIONS:
            fuel_default = ""

        fuel_var = tk.StringVar()
        fuel_box = ttk.Combobox(frame, textvariable=fuel_var, values=FUEL_OPTIONS, width=8)
        fuel_box.grid(row=r, column=4, padx=4)

        if fuel_default:
            fuel_box.set(fuel_default)
            fuel_box.configure(state="readonly")
        else:
            fuel_box.set("")
            fuel_box.configure(state="readonly")


        # -------- PV type --------
        pv_default = b["PVType"]
        if pv_default not in PV_TYPES:
            pv_default = PV_TYPES[0]

        pv_var = tk.StringVar()
        pv_box = ttk.Combobox(frame, textvariable=pv_var, values=PV_TYPES, state="readonly", width=6)
        pv_box.grid(row=r, column=5, padx=4)
        pv_box.set(pv_default)


        # -------- PV slider --------
        
        baseline = float(b["PVpercentage"])
        
        pv_slider_var = tk.DoubleVar(value=baseline)
        
        pv_slider = ttk.Scale(
            frame,
            from_=baseline,      # left limit = baseline
            to=100,
            orient="horizontal",
            length=120,
            variable=pv_slider_var
        )
        
        pv_slider.grid(row=r, column=6, padx=4)
        
        pv_label = ttk.Label(frame, text=f"{baseline:.0f}%")
        pv_label.grid(row=r, column=7, padx=4)
        
        
        def update_label(value, label=pv_label):
            label.config(text=f"{float(value):.0f}%")
        
        pv_slider.configure(command=update_label)
        
        
        if pv_default == "none":
            pv_slider.configure(state="disabled")
            pv_slider_var.set(0)
            pv_label.config(text="0%")
            
        rows[r] = dict(
            idx=idx,
        
            env_var=env_var,
            env_box=env_box,
        
            heat_var=heat_var,
            heat_box=heat_box,
        
            dhw_var=dhw_var,
            dhw_box=dhw_box,
        
            fuel_var=fuel_var,
            fuel_box=fuel_box,
        
            pv_var=pv_var,
            pv_box=pv_box,
        
            pv_slider=pv_slider,
            pv_slider_var=pv_slider_var
        )
    return rows

def attach_callbacks(rows):
    pass
    # def update_fuel_state(r):
    
    #     heat = rows[r]["heat_var"].get()
    #     dhw = rows[r]["dhw_var"].get()
    
    #     fuel_box = rows[r]["fuel_box"]
    
    #     if heat == "boiler" or dhw == "boiler":
    #         fuel_box.configure(state="readonly")
    #     else:
    #         fuel_box.configure(state="disabled")


    # def heating_changed(event, r):
    #     update_fuel_state(r)


    # def dhw_changed(event, r):
    #     update_fuel_state(r)


    # for r in rows:

    #     rows[r]["heat_box"].bind(
    #         "<<ComboboxSelected>>",
    #         lambda e, r=r: heating_changed(e, r)
    #     )

    #     rows[r]["dhw_box"].bind(
    #         "<<ComboboxSelected>>",
    #         lambda e, r=r: dhw_changed(e, r)
    #     )

    # # IMPORTANT: initialize after bindings and after defaults exist
    # for r in rows:
    #     update_fuel_state(r)
        
        


def run_interface(root, rows, gdf):

    import pandas as pd
    from tkinter import ttk

    result = {"gdf": None}

    def save():

        for r, row in rows.items():
    
            idx = row["idx"]
    
            env  = row["env_box"].get()
            heat = row["heat_box"].get()
            dhw  = row["dhw_box"].get()
            fuel = row["fuel_box"].get()
    
            pv_type = row["pv_box"].get()
            pv_perc = float(row["pv_slider"].get())
    
            # Envelope
            gdf.at[idx, "EEdepth"] = env
    
            # Heating
            if heat == "boiler":
                gdf.at[idx, "SHSource"] = f"boiler_{fuel}"
            else:
                gdf.at[idx, "SHSource"] = heat
    
            # DHW
            if dhw == "boiler":
                gdf.at[idx, "DHWsource"] = f"boiler_{fuel}"
            else:
                gdf.at[idx, "DHWsource"] = dhw
    
            # PV
            gdf.at[idx, "PVType"] = pv_type
            gdf.at[idx, "PVpercentage"] = max(gdf.at[idx, "PVpercentage"], pv_perc)
    
    
        result["gdf"] = gdf
    
        root.quit()
        root.destroy()

    button_frame = ttk.Frame(root)
    button_frame.pack(pady=10)

    ttk.Button(
        button_frame,
        text="Save configuration",
        command=save
    ).pack()

    root.mainloop()

    return result

def translate_source(source):
    """
    Translate the GeoJSON heating/DHW source to the UI option.

    Examples
    --------
    boiler_gas  -> boiler
    boiler_bio  -> boiler
    hp_he       -> hp_he
    hp_me       -> hp_me
    hp_le       -> hp_le
    dhn         -> dhn
    """

    if not isinstance(source, str):
        return None

    if source.startswith("boiler"):
        return "boiler"

    if source in ["dhn", "hp_le", "hp_me", "hp_he"]:
        return source

    return None

def extract_fuel(source):
    """
    Extract boiler fuel from GeoJSON source.

    Examples
    --------
    boiler_gas -> gas
    boiler_bio -> bio
    """

    if isinstance(source, str) and source.startswith("boiler_"):
        return source.split("_", 1)[1]

    return ""