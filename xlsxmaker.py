from __future__ import annotations

# -*- coding: utf-8 -*-\
"""
Created on Fri Dec 12 09:27:53 2025

@author: khajmoh18975
"""

#!/usr/bin/env python3
"""
GeoJSON -> Excel (SWEREF 99 TM EPSG:3006) in the format:
Building ID | Axis (X/Y) | Point 1 | Point 2 | ... | Point N

- Supports Polygon, MultiPolygon, LineString, MultiLineString, Point, MultiPoint.
- For Polygon/MultiPolygon: uses the exterior ring of the largest polygon part.
- Drops the duplicated closing coordinate of polygons.
- If input CRS is missing in the file, assumes EPSG:4326 unless you pass --input-crs.
"""


import argparse
from typing import List, Tuple, Optional

import geopandas as gpd
import pandas as pd
from shapely.geometry import (
    Polygon, MultiPolygon,
    LineString, MultiLineString,
    Point, MultiPoint,
)
from shapely.geometry.base import BaseGeometry


SWEREF_EPSG = 3006


def _largest_polygon_part(geom: MultiPolygon) -> Polygon:
    parts = list(geom.geoms)
    return max(parts, key=lambda p: p.area)


def extract_ordered_coords(geom: BaseGeometry) -> List[Tuple[float, float]]:
    """Return ordered (x,y) coords for a geometry."""
    if geom is None or geom.is_empty:
        return []

    if isinstance(geom, Polygon):
        coords = list(geom.exterior.coords)
        # drop duplicate closing point if present
        if len(coords) >= 2 and coords[0] == coords[-1]:
            coords = coords[:-1]
        return [(float(x), float(y)) for x, y in coords]

    if isinstance(geom, MultiPolygon):
        poly = _largest_polygon_part(geom)
        return extract_ordered_coords(poly)

    if isinstance(geom, LineString):
        return [(float(x), float(y)) for x, y in list(geom.coords)]

    if isinstance(geom, MultiLineString):
        # concatenate parts in order (keeps all points)
        out: List[Tuple[float, float]] = []
        for ls in geom.geoms:
            out.extend(extract_ordered_coords(ls))
        return out

    if isinstance(geom, Point):
        return [(float(geom.x), float(geom.y))]

    if isinstance(geom, MultiPoint):
        return [(float(p.x), float(p.y)) for p in geom.geoms]

    # fallback: try generic coords access if possible
    try:
        return [(float(x), float(y)) for x, y in list(geom.coords)]
    except Exception:
        return []


def geojson_to_excel(
    geojson_path: str,
    excel_path: str,
    id_field: Optional[str] = None,
    input_crs: Optional[str] = None,
):
    gdf = gpd.read_file(geojson_path)

    # Ensure CRS is known
    if gdf.crs is None:
        if input_crs is None:
            input_crs = "EPSG:4326"  # common default for GeoJSON
        gdf = gdf.set_crs(input_crs)

    # Reproject to SWEREF 99 TM
    gdf = gdf.to_crs(epsg=SWEREF_EPSG)

    # Decide building IDs
    if id_field and id_field in gdf.columns:
        building_ids = gdf[id_field].astype(str).tolist()
    else:
        # 1..N like your screenshot
        building_ids = [str(i + 1) for i in range(len(gdf))]

    # Extract coords per feature
    coords_per_feat: List[List[Tuple[float, float]]] = [
        extract_ordered_coords(geom) for geom in gdf.geometry
    ]

    max_points = max((len(c) for c in coords_per_feat), default=0)
    point_cols = [f"Point {i}" for i in range(1, max_points + 1)]

    rows = []
    for bid, coords in zip(building_ids, coords_per_feat):
        xs = [x for x, _ in coords]
        ys = [y for _, y in coords]

        row_x = {"Building ID": bid, "Axis": "X"}
        row_y = {"Building ID": bid, "Axis": "Y"}

        for i in range(max_points):
            row_x[f"Point {i+1}"] = xs[i] if i < len(xs) else None
            row_y[f"Point {i+1}"] = ys[i] if i < len(ys) else None

        rows.append(row_x)
        rows.append(row_y)

    df = pd.DataFrame(rows, columns=["Building ID", "Axis"] + point_cols)

    # Write Excel
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="Points_SWEREF3006")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("geojson", help="Input GeoJSON path")
    ap.add_argument("xlsx", help="Output Excel path (.xlsx)")
    ap.add_argument("--id-field", default=None, help="Column name to use as Building ID")
    ap.add_argument(
        "--input-crs",
        default=None,
        help="Only used if the GeoJSON has no CRS. Example: EPSG:4326",
    )
    args = ap.parse_args()

    geojson_to_excel(
        geojson_path=args.geojson,
        excel_path=args.xlsx,
        id_field=args.id_field,
        input_crs=args.input_crs,
    )



if __name__ == "__main__":
    geojson_to_excel(
        geojson_path=r"C:/works/EUReCA/EUReCA/eureca_ubem/Input/grasloken_rip.geojson",
        excel_path=r"C:/works/EUReCA/EUReCA/eureca_ubem/Input/grasloken_rip.xlsx",
        id_field="building_id",
        input_crs=None,  # or "EPSG:4326" if the geojson has no CRS
    )
