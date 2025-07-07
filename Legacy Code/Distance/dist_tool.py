import json
import os
import sys
import time
from typing import Tuple, List

import geopandas as gpd
from shapely.geometry import MultiPolygon, Polygon, Point
from shapely.ops import unary_union

# Import necessary functions and data from dist_updated.py
# Ensure that dist_updated.py is in the same directory as this script
from dist_updated import (
    VALID_COUNTRIES,
    SYNONYM_MAP,
    SAMPLE_SIZE_MAP,
    DEFAULT_SAMPLE_SIZE,
    MICRONATION_COORDS,
    get_sample_size,
    normalize_name,
    distance_polygon_to_polygon,
    distance_multiple_points_to_multiple_points,
    distance_multiple_points_to_polygon,
    geodesic_distance,
    filter_polygons
)

# Constants
SHAPEFILE_NAME = "ne_110m_admin_0_countries.shp"
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_country_geometries(shapefile_path: str) -> dict:
    """
    Load the shapefile and process country geometries.
    Returns a dictionary mapping country names to their geometries.
    """
    print("Loading shapefile...")
    gdf = gpd.read_file(shapefile_path)

    # Ensure geometry is in EPSG:4326
    if gdf.crs and gdf.crs.to_string() != "EPSG:4326":
        print("Reprojecting to EPSG:4326...")
        gdf = gdf.to_crs(epsg=4326)

    # Potential columns that might contain the country name
    possible_name_cols = [
        "ADMIN", "NAME", "NAME_LONG", "SOVEREIGNT", "BRK_NAME",
        "FORMAL_EN", "GEOUNIT", "GU_A3", "ISO_A3", "ISO_A2"
    ]
    found_cols = [col for col in possible_name_cols if col in gdf.columns]
    if not found_cols:
        raise ValueError("No known name columns found in the shapefile! "
                         "Please update 'possible_name_cols' or your data.")

    def extract_valid_country(row_series):
        """
        Try each column in found_cols; return first valid country name or "".
        """
        for col in found_cols:
            raw_val = str(row_series[col]).strip()
            normed = normalize_name(raw_val)
            if normed in VALID_COUNTRIES:
                return normed
        return ""

    print(f"Found {len(gdf)} rows in the shapefile. Processing...")
    country_geoms = {}

    for idx, row in gdf.iterrows():
        country_name = extract_valid_country(row)
        if not country_name:
            continue
        raw_geom = row.geometry
        filtered = filter_polygons(raw_geom, country_name)
        if filtered and not filtered.is_empty:
            if country_name not in country_geoms:
                country_geoms[country_name] = filtered
            else:
                combined = unary_union([country_geoms[country_name], filtered])
                country_geoms[country_name] = combined

    # Ensure all 196 countries are present
    final_polygons = {}
    for ctry in VALID_COUNTRIES:
        if ctry in country_geoms:
            final_polygons[ctry] = country_geoms[ctry]
        else:
            final_polygons[ctry] = None

    return final_polygons


def update_sample_sizes():
    """
    Allow the user to update sample sizes for specific countries.
    """
    print("\nDo you want to override sample sizes? (yes/no)")
    choice = input("Enter choice: ").strip().lower()
    if choice not in ['yes', 'y']:
        return

    while True:
        country_input = input("\nEnter country name to update sample size (or 'done' to finish): ").strip().lower()
        if country_input == 'done':
            break
        normalized = normalize_name(country_input)
        if normalized not in VALID_COUNTRIES:
            print("Invalid country name. Please try again.")
            continue
        try:
            size_input = int(input(f"Enter new sample size for '{normalized}': ").strip())
            if size_input <= 0:
                print("Sample size must be a positive integer.")
                continue
            SAMPLE_SIZE_MAP[normalized] = size_input
            print(f"Sample size for '{normalized}' updated to {size_input}.")
        except ValueError:
            print("Invalid input. Please enter a valid integer.")


def get_coordinates(prompt: str) -> Tuple[float, float]:
    """
    Prompt the user to enter coordinates in 'lon, lat' format.
    """
    while True:
        coord_input = input(f"{prompt} (format: lon, lat): ").strip()
        try:
            lon, lat = map(float, coord_input.split(','))
            if not (-180 <= lon <= 180 and -90 <= lat <= 90):
                print("Coordinates out of bounds. Longitude must be between -180 and 180, "
                      "and latitude between -90 and 90.")
                continue
            return lon, lat
        except ValueError:
            print("Invalid format. Please enter coordinates as 'lon, lat'.")


def get_country(prompt: str) -> str:
    """
    Prompt the user to enter a country name and normalize it.
    """
    while True:
        country_input = input(f"{prompt}: ").strip().lower()
        normalized = normalize_name(country_input)
        if normalized in VALID_COUNTRIES:
            return normalized
        else:
            print("Invalid or unrecognized country name. Please try again.")


def save_json(data: dict, filename: str):
    """
    Save the given data dictionary to a JSON file with the specified filename
    in the output directory.
    """
    filepath = os.path.join(OUTPUT_DIR, filename)
    try:
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        print(f"Results saved to '{filepath}'.")
    except Exception as e:
        print(f"Failed to save JSON file: {e}")


def coordinates_to_coordinates():
    """
    Option 1: Calculate distance between two coordinates.
    """
    print("\n--- Coordinates to Coordinates ---")
    lon1, lat1 = get_coordinates("Enter first coordinate")
    lon2, lat2 = get_coordinates("Enter second coordinate")
    distance = geodesic_distance(lon1, lat1, lon2, lat2)
    print(f"Distance: {distance:.1f} km")


def coordinates_to_single_country(final_polygons: dict):
    """
    Option 2: Calculate distance from coordinates to a single country.
    """
    print("\n--- Coordinates to Single Country ---")
    lon, lat = get_coordinates("Enter coordinate")
    country = get_country("Enter target country name")
    poly = final_polygons.get(country)
    if not poly:
        print(f"No geometry found for '{country}'.")
        return
    sample_size = get_sample_size(country)
    point = (lon, lat)
    if country in MICRONATION_COORDS:
        # If the country is a micronation with coordinates
        points = MICRONATION_COORDS[country]
        distance = distance_multiple_points_to_polygon(points, poly, sample_size)
    else:
        distance = distance_multiple_points_to_polygon([point], poly, sample_size)
    if distance is not None:
        print(f"Distance from point to '{country}': {distance:.1f} km")
    else:
        print("Unable to calculate distance.")


def coordinates_to_all_countries(final_polygons: dict):
    """
    Option 3: Calculate distances from coordinates to all countries and save to JSON.
    """
    print("\n--- Coordinates to All Countries ---")
    lon, lat = get_coordinates("Enter coordinate")
    point = (lon, lat)
    distances = {}
    for country in sorted(VALID_COUNTRIES):
        poly = final_polygons.get(country)
        sample_size = get_sample_size(country)
        if poly:
            distance = distance_multiple_points_to_polygon([point], poly, sample_size)
        elif country in MICRONATION_COORDS:
            points = MICRONATION_COORDS[country]
            distance = distance_multiple_points_to_multiple_points([point], points)
        else:
            distance = None
        distances[country] = round(distance, 1) if distance is not None else 9999999

    print("Distances calculated.")
    output_name = input("Enter name for the output JSON file (e.g., 'distances.json'): ").strip()
    if not output_name.endswith('.json'):
        output_name += '.json'
    save_json(distances, output_name)


def single_country_to_single_country(final_polygons: dict):
    """
    Option 4: Calculate distance between two countries.
    """
    print("\n--- Single Country to Single Country ---")
    country1 = get_country("Enter first country name")
    country2 = get_country("Enter second country name")
    if country1 == country2:
        print("Both countries are the same. Distance is 0 km.")
        return
    poly1 = final_polygons.get(country1)
    poly2 = final_polygons.get(country2)
    if not poly1:
        print(f"No geometry found for '{country1}'.")
        return
    if not poly2:
        print(f"No geometry found for '{country2}'.")
        return
    sample_size1 = get_sample_size(country1)
    sample_size2 = get_sample_size(country2)
    if poly1 and poly2:
        distance = distance_polygon_to_polygon(poly1, poly2, sample_size1, sample_size2)
    elif (country1 in MICRONATION_COORDS) and (country2 in MICRONATION_COORDS):
        points1 = MICRONATION_COORDS[country1]
        points2 = MICRONATION_COORDS[country2]
        distance = distance_multiple_points_to_multiple_points(points1, points2)
    elif poly1 and (country2 in MICRONATION_COORDS):
        points2 = MICRONATION_COORDS[country2]
        distance = distance_multiple_points_to_polygon(points2, poly1, sample_size1)
    elif poly2 and (country1 in MICRONATION_COORDS):
        points1 = MICRONATION_COORDS[country1]
        distance = distance_multiple_points_to_polygon(points1, poly2, sample_size2)
    else:
        distance = None

    if distance is not None:
        print(f"Distance between '{country1}' and '{country2}': {distance:.1f} km")
    else:
        print("Unable to calculate distance.")


def single_country_to_all_countries(final_polygons: dict):
    """
    Option 5: Calculate distances from a single country to all countries and save to JSON.
    """
    print("\n--- Single Country to All Countries ---")
    source_country = get_country("Enter source country name")
    poly_source = final_polygons.get(source_country)
    sample_size_source = get_sample_size(source_country)
    if not poly_source and source_country not in MICRONATION_COORDS:
        print(f"No geometry or coordinates found for '{source_country}'.")
        return

    distances = {}
    for target_country in sorted(VALID_COUNTRIES):
        if target_country == source_country:
            distances[target_country] = 0.0
            continue
        poly_target = final_polygons.get(target_country)
        sample_size_target = get_sample_size(target_country)
        if poly_source and poly_target:
            distance = distance_polygon_to_polygon(poly_source, poly_target, sample_size_source, sample_size_target)
        elif (source_country in MICRONATION_COORDS) and (target_country in MICRONATION_COORDS):
            points1 = MICRONATION_COORDS[source_country]
            points2 = MICRONATION_COORDS[target_country]
            distance = distance_multiple_points_to_multiple_points(points1, points2)
        elif poly_source and (target_country in MICRONATION_COORDS):
            points2 = MICRONATION_COORDS[target_country]
            distance = distance_multiple_points_to_polygon(points2, poly_source, sample_size_source)
        elif poly_target and (source_country in MICRONATION_COORDS):
            points1 = MICRONATION_COORDS[source_country]
            distance = distance_multiple_points_to_polygon(points1, poly_target, sample_size_target)
        else:
            distance = None
        distances[target_country] = round(distance, 1) if distance is not None else 9999999

    print("Distances calculated.")
    output_name = input("Enter name for the output JSON file (e.g., 'country_distances.json'): ").strip()
    if not output_name.endswith('.json'):
        output_name += '.json'
    save_json(distances, output_name)


def main():
    # Load country geometries
    shapefile_path = os.path.join(OUTPUT_DIR, SHAPEFILE_NAME)
    if not os.path.exists(shapefile_path):
        print(f"Shapefile '{SHAPEFILE_NAME}' not found in '{OUTPUT_DIR}'. Exiting.")
        sys.exit(1)

    final_polygons = load_country_geometries(shapefile_path)

    # Allow user to update sample sizes
    update_sample_sizes()

    # Menu loop
    while True:
        print("\n=== Distance Calculation Tool ===")
        print("Select an option:")
        print("1. Coordinates to Coordinates")
        print("2. Coordinates to Single Country")
        print("3. Coordinates to All Countries")
        print("4. Single Country to Single Country")
        print("5. Single Country to All Countries")
        print("6. Exit")

        choice = input("Enter your choice (1-6): ").strip()
        if choice == '1':
            coordinates_to_coordinates()
        elif choice == '2':
            coordinates_to_single_country(final_polygons)
        elif choice == '3':
            coordinates_to_all_countries(final_polygons)
        elif choice == '4':
            single_country_to_single_country(final_polygons)
        elif choice == '5':
            single_country_to_all_countries(final_polygons)
        elif choice == '6':
            print("Exiting the tool. Goodbye!")
            break
        else:
            print("Invalid choice. Please enter a number between 1 and 6.")


if __name__ == "__main__":
    main()
