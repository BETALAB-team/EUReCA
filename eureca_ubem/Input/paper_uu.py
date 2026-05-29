"""
run_feasibility_analysis.py

Main runner for the scenario-swarm feasibility analysis.

This script uses the same path structure as the previous PUBEM orchestrator,
but it calls feasibility_orchestrator.py instead of the market/game workflow.

Run from the project environment where eureca_building, eureca_ubem, and
eureca_pubem are importable.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from time import time

import pandas as pd
import geopandas as gpd

from eureca_building.config import load_config
from eureca_ubem.city import City
from eureca_pubem import scenario_process as sc

from feasibility_orchestrator import (
    CostConfig,
    DiffusionPoint,
    InterventionPolicy,
    TechnicalLimits,
    build_diffusion_grid,
    run_convergence_analysis,
    run_feasibility_swarm,
)

warnings.filterwarnings("ignore")


PATHS = {
    "config": r".\Example_District_Config.json",
    "weather": r".\SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",
    "schedules": r"TransitionMatrixSchedule.csv",
    "materials": r".\tabula_sverige.xlsx",
    "city_model": r".\soderman_limited_reproject.geojson",
    "systems": r".\systems.xlsx",
    "output_folder": r".\feasibility_outputs",

    "distribution": r"C:\Works\EUReCA\EUReCA\eureca_ubem\Input\distribution.geojson",
    "baseline_geojson": r"C:\Works\EUReCA\EUReCA\eureca_ubem\Input\soderman_limited_reproject.geojson",
    "roads": r"C:\Works\EUReCA\EUReCA\eureca_ubem\Input\roads.geojson",
    "weather_abs": r"C:\Works\EUReCA\EUReCA\eureca_ubem\Input\SWE_UP_Uppsala.Univ.024620_TMYx.2009-2023.epw",

    "dhn_pipes": r"C:\Works\EUReCA\EUReCA\eureca_pubem\dhn\_pipe_diameters.json",
    "cables": r"C:\Works\EUReCA\EUReCA\eureca_pubem\grid\swedish_cable.json",
}

def safe_path_name(text: str) -> str:
    """
    Make a string safe for Windows folder/file names.
    """
    bad_chars = '<>:"/\\|?*'
    out = str(text)

    for ch in bad_chars:
        out = out.replace(ch, "_")

    out = out.replace("=", "-")
    out = out.replace(" ", "_")

    return out

def load_buildings(input_data):
    if isinstance(input_data, gpd.GeoDataFrame):
        gdf = input_data.copy()
    elif isinstance(input_data, str):
        gdf = gpd.read_file(input_data)
    else:
        raise ValueError("invalid_input")

    required_columns = [
        "id",
        "Name",
        "EEdepth",
        "SHSource",
        "DHWsource",
        "PVType",
        "PVpercentage",
    ]

    missing = [c for c in required_columns if c not in gdf.columns]

    if missing:
        raise ValueError(f"missing_columns: {missing}")

    return gdf


def simulate_city_envelopes():
    """
    Runs EUReCA-UBEM once for each envelope level and stores the hourly
    building demand states used later by the scenario swarm.
    """

    output_folder = PATHS["output_folder"]
    os.makedirs(output_folder, exist_ok=True)

    load_config(PATHS["config"])

    cities_envelopes = {}
    last_gdf = None

    for depth, envelope in enumerate(["none", "shallow", "medium", "deep"]):
        print(f"\nSimulating envelope state: {envelope}")

        gdf = sc.load_geojson_with_envelope_prefix(PATHS["city_model"], depth)

        scenario_city = City(
            city_model=gdf,
            epw_weather_file=PATHS["weather"],
            end_uses_types_file=PATHS["schedules"],
            envelope_types_file=PATHS["materials"],
            systems_templates_file=PATHS["systems"],
            shading_calculation=True,
            building_model="2C",
            output_folder=output_folder,
        )

        scenario_city.simulate()
        last_gdf, value_store = sc.create_gdf_dictionary(scenario_city)
        cities_envelopes[envelope] = value_store

    return last_gdf, cities_envelopes


def create_baseline_for_feasibility(mycity_gdf, city_demand_states):
    """
    Creates the baseline scenario object needed by the feasibility orchestrator.

    This mirrors the previous create_initial_baseline() logic, but stops before
    market optimization and building NPV response.
    """

    baseline = sc.create_baseline(
        input_gdf=mycity_gdf,
        input_city=city_demand_states,
        distributions_geojson=PATHS["distribution"],
        baseline_gdf=PATHS["baseline_geojson"],
        weather_path=PATHS["weather_abs"],
        street_path=PATHS["roads"],
    )

    baseline_gdf = load_buildings(PATHS["baseline_geojson"])

    return baseline, baseline_gdf


def make_policy() -> InterventionPolicy:
    return InterventionPolicy(
        target_envelope="deep",
        target_hp_source="hp_he",
        target_dhn_source="dhn",
        target_pv_type="A",
        target_pv_percentage=100.0,
        apply_hp_to="both",
        apply_dhn_to="both",
        conflict_mode="disjoint_hp_dhn",
        allow_envelope_downgrade=False,
        allow_pv_decrease=False,
    )


def make_limits() -> TechnicalLimits:
    return TechnicalLimits(
        min_voltage_pu=0.95,
        max_voltage_pu=1.05,
        max_supply_apparent_power=1_000_000.0,
        max_supply_heat=500_000.0,
        max_pipe_pressure=2_000_000.0,
        max_dp_per_m=200.0,
        max_pump_power=1_000.0,
        required_supply_temperature=60.0,
    )


def make_cost_config() -> CostConfig:
    return CostConfig(
        dhn_area_type="urban",
        grid_area_type="town",
        dhn_assumptions={
            "replacement_factor_ground": 0.5,
            "removal_factor": 0.3,
        },
        grid_assumptions={
            "replacement_factor_ground": 1.0,
            "removal_factor": 0.2,
            "ground_share": 0.55,
            "rest_share": 0.45,
        },
        grid_cable_key="name",
    )


def make_smoke_test_diffusion_points():
    """
    Small run to verify that the full pipeline works.
    Use this first.
    """

    return [
        DiffusionPoint(pv=0.0, hp=0.0, dhn=0.0, envelope=0.0),
        DiffusionPoint(pv=0.5, hp=0.5, dhn=0.0, envelope=0.5),
        DiffusionPoint(pv=1.0, hp=1.0, dhn=0.0, envelope=1.0),
        DiffusionPoint(pv=0.5, hp=0.0, dhn=0.5, envelope=0.5),
    ]


def make_full_diffusion_grid():
    """
    Full grid. This can become expensive fast.
    Start with the smoke test before using this.
    """

    return build_diffusion_grid(
        pv_levels=[0.0, 0.25, 0.50, 0.75, 1.0],
        hp_levels=[0.0, 0.25, 0.50, 0.75, 1.0],
        dhn_levels=[0.0, 0.25, 0.50, 0.75, 1.0],
        envelope_levels=[0.0, 0.50, 1.0],
    )


def make_convergence_points():
    """
    Representative non-boundary stress points for convergence analysis.

    For each pair:
    - low-low: 0.2, 0.2
    - low-high: 0.2, 0.8
    - high-low: 0.8, 0.2
    - high-high: 0.8, 0.8

    PV-HP plane:
        vary pv and hp, fix dhn=0 and envelope=0

    DHN-envelope plane:
        vary dhn and envelope, fix pv=0 and hp=0
    """

    return [
        # PV-HP convergence stress points
        DiffusionPoint(pv=0.2, hp=0.2, dhn=0.0, envelope=0.0),
        DiffusionPoint(pv=0.2, hp=0.8, dhn=0.0, envelope=0.0),
        DiffusionPoint(pv=0.8, hp=0.2, dhn=0.0, envelope=0.0),
        DiffusionPoint(pv=0.8, hp=0.8, dhn=0.0, envelope=0.0),

        # DHN-envelope convergence stress points
        DiffusionPoint(pv=0.0, hp=0.0, dhn=0.2, envelope=0.2),
        DiffusionPoint(pv=0.0, hp=0.0, dhn=0.2, envelope=0.8),
        DiffusionPoint(pv=0.0, hp=0.0, dhn=0.8, envelope=0.2),
        DiffusionPoint(pv=0.0, hp=0.0, dhn=0.8, envelope=0.8),
    ]


def select_r_star(convergence_df: pd.DataFrame, r_values: list[int]) -> int:
    """
    Select final R* from all convergence stress points.

    Preferred logic:
    - if convergence_df contains selected_r, take max selected_r
    - else if it contains converged/r, take max r among converged rows
    - otherwise fall back to max tested R
    """

    if convergence_df.empty:
        return max(r_values)

    if "selected_r" in convergence_df.columns:
        valid = pd.to_numeric(convergence_df["selected_r"], errors="coerce").dropna()

        if not valid.empty:
            return int(valid.max())

    if "selected_R" in convergence_df.columns:
        valid = pd.to_numeric(convergence_df["selected_R"], errors="coerce").dropna()

        if not valid.empty:
            return int(valid.max())

    if "converged" in convergence_df.columns and "r" in convergence_df.columns:
        converged = convergence_df[convergence_df["converged"] == True]
        valid = pd.to_numeric(converged["r"], errors="coerce").dropna()

        if not valid.empty:
            return int(valid.max())

    if "converged" in convergence_df.columns and "R" in convergence_df.columns:
        converged = convergence_df[convergence_df["converged"] == True]
        valid = pd.to_numeric(converged["R"], errors="coerce").dropna()

        if not valid.empty:
            return int(valid.max())

    print("\nWARNING: No explicit convergence-selected R found.")
    print(f"Falling back to max tested R = {max(r_values)}")

    return max(r_values)


def run_all_convergence_points(
    *,
    convergence_points,
    r_values,
    baseline_geojson,
    city_demand_states,
    baseline_scenario,
    weatherfile_path,
    policy,
    limits,
    cost_config,
    dhn_pipe_cost_json,
    grid_cable_cost_json,
    output_folder,
    base_seed,
    verbose,
):
    """
    Runs convergence analysis for all selected stress points and concatenates results.
    """

    convergence_tables = []

    for i, dp in enumerate(convergence_points):
        print("\n" + "=" * 100)
        print(f"Running convergence point {i + 1}/{len(convergence_points)}")
        print(f"Diffusion point: {dp}")
        print("=" * 100)

        df = run_convergence_analysis(
            diffusion=dp,
            r_values=r_values,
            baseline_geojson=baseline_geojson,
            city_demand_states=city_demand_states,
            baseline_scenario=baseline_scenario,
            weatherfile_path=weatherfile_path,
            policy=policy,
            limits=limits,
            cost_config=cost_config,
            dhn_pipe_cost_json=dhn_pipe_cost_json,
            grid_cable_cost_json=grid_cable_cost_json,
            indicator_name="combined_stress_mean",
            output_folder=str(Path(output_folder) / "convergence" / safe_path_name(dp.as_key())),
            base_seed=base_seed + i * 10_000,
            verbose=verbose,
        )

        df = df.copy()
        df["convergence_point_index"] = i
        df["diffusion_key"] = dp.as_key()
        df["pv"] = dp.pv
        df["hp"] = dp.hp
        df["dhn"] = dp.dhn
        df["envelope"] = dp.envelope

        convergence_tables.append(df)

    if not convergence_tables:
        return pd.DataFrame()

    return pd.concat(convergence_tables, ignore_index=True)


def main():
    start = time()

    run_mode = "smoke"  # "smoke" first, then change to "full"

    # Smoke keeps the swarm small even if convergence suggests a larger R.
    smoke_n_realisations = 3

    # These are the tested R values for convergence.
    # For final runs, consider list(range(20, 151, 10)) or list(range(20, 201, 10)).
    r_values = [5, 10, 20, 40, 80] if run_mode == "smoke" else list(range(20, 151, 10))

    print("\nStarting feasibility-swarm analysis.")
    print(f"Run mode: {run_mode}")
    print(f"Convergence R values: {r_values}")

    mycity_gdf, city_demand_states = simulate_city_envelopes()
    print("\nEnvelope simulations completed.")

    baseline_scenario, baseline_gdf = create_baseline_for_feasibility(
        mycity_gdf=mycity_gdf,
        city_demand_states=city_demand_states,
    )
    print("Baseline scenario created.")

    if run_mode == "smoke":
        diffusion_points = make_smoke_test_diffusion_points()
    else:
        diffusion_points = make_full_diffusion_grid()

    policy = make_policy()
    limits = make_limits()
    cost_config = make_cost_config()

    output_folder = Path(PATHS["output_folder"])
    output_folder.mkdir(parents=True, exist_ok=True)

    print("\nRunning convergence analysis on representative non-boundary stress points.")

    convergence_points = make_convergence_points()

    convergence_df = run_all_convergence_points(
        convergence_points=convergence_points,
        r_values=r_values,
        baseline_geojson=PATHS["baseline_geojson"],
        city_demand_states=city_demand_states,
        baseline_scenario=baseline_scenario,
        weatherfile_path=PATHS["weather_abs"],
        policy=policy,
        limits=limits,
        cost_config=cost_config,
        dhn_pipe_cost_json=PATHS["dhn_pipes"],
        grid_cable_cost_json=PATHS["cables"],
        output_folder=output_folder,
        base_seed=5000,
        verbose=True,
    )

    print("\nConvergence results:")
    print(convergence_df)

    convergence_df.to_csv(output_folder / "convergence_all_points.csv", index=False)

    r_star = select_r_star(convergence_df, r_values)

    print("\n" + "=" * 100)
    print(f"Selected R_star from convergence = {r_star}")
    print("=" * 100)

    if run_mode == "smoke":
        n_realisations = smoke_n_realisations
        print(f"\nSmoke mode active. Using n_realisations={n_realisations}, not R_star.")
    else:
        n_realisations = r_star
        print(f"\nFull mode active. Using n_realisations=R_star={n_realisations}.")

    print("\nRunning scenario-swarm feasibility analysis.")

    realisations_df, aggregates_df, _ = run_feasibility_swarm(
        diffusion_points=diffusion_points,
        n_realisations=n_realisations,
        baseline_geojson=PATHS["baseline_geojson"],
        city_demand_states=city_demand_states,
        baseline_scenario=baseline_scenario,
        weatherfile_path=PATHS["weather_abs"],
        policy=policy,
        limits=limits,
        cost_config=cost_config,
        dhn_pipe_cost_json=PATHS["dhn_pipes"],
        grid_cable_cost_json=PATHS["cables"],
        output_folder=str(output_folder),
        base_seed=1000,
        verbose=True,
    )

    print("\nDone.")
    print(f"Realisation rows: {len(realisations_df)}")
    print(f"Aggregate rows: {len(aggregates_df)}")
    print(f"Outputs written to: {output_folder.resolve()}")

    end = time()
    print(f"Finished in {(end - start) / 60:.2f} minutes")

    return {
        "city_demand_states": city_demand_states,
        "baseline_scenario": baseline_scenario,
        "baseline_gdf": baseline_gdf,
        "convergence_df": convergence_df,
        "r_star": r_star,
        "realisations_df": realisations_df,
        "aggregates_df": aggregates_df,
    }


if __name__ == "__main__":
    results = main()