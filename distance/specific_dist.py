import os
import json
import sys
from dist_updated import (
    geodesic_distance,
    distance_point_to_polygon,
    distance_polygon_to_polygon,
    normalize_name,
    VALID_COUNTRIES,
    MICRONATION_COORDS,
    SAMPLE_SIZE_MAP,
    DEFAULT_SAMPLE_SIZE,
    get_sample_size,
    filter_polygons,
)
import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import unary_union

# ------------------- Configuration -------------------
# Directory paths
BASE_DIR = r"C:\Users\raksh\My Drive\Tech\deploy\distance"
SHAPEFILE_NAME = "ne_110m_admin_0_countries.shp"
SHAPEFILE_PATH = os.path.join(BASE_DIR, SHAPEFILE_NAME)

# Optional: Override sample sizes here if desired
# SAMPLE_SIZE_OVERRIDE = {
#     'united states': 500,
#     'canada': 500,
# }
SAMPLE_SIZE_OVERRIDE = {}  # Empty by default

# -----------------------------------------------------

def load_country_geometries():
    """
    Load and process country geometries from the shapefile.
    Returns a dictionary mapping normalized country names to their geometries.
    """
    print("Loading shapefile...")
    try:
        gdf = gpd.read_file(SHAPEFILE_PATH)
    except FileNotFoundError:
        print(f"Shapefile not found at '{SHAPEFILE_PATH}'. Please check the path.")
        sys.exit(1)

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
        print("No known name columns found in the shapefile!")
        sys.exit(1)

    def extract_valid_country(row_series):
        for col in found_cols:
            raw_val = str(row_series[col]).strip()
            normed = normalize_name(raw_val)
            if normed in VALID_COUNTRIES:
                return normed
        return ""

    print("Processing country geometries...")
    country_geoms = {}
    total_rows = len(gdf)
    for idx, (_, row) in enumerate(gdf.iterrows(), start=1):
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

    # Ensure all valid countries are present
    final_polygons = {}
    for country in VALID_COUNTRIES:
        final_polygons[country] = country_geoms.get(country, None)

    print("Country geometries loaded successfully.\n")
    return final_polygons

def get_user_choice():
    """
    Display the menu and get the user's choice.
    """
    print("What would you like to find?")
    print("1. Coordinates to Coordinates Distance")
    print("2. Coordinates to Country Distance")
    print("3. Coordinates to All Countries Distance")
    print("4. Country to Country Distance")
    print("5. Country to All Countries Distance")
    choice = input("Enter the number corresponding to your choice (1-5): ").strip()
    return choice

def get_coordinates(prompt="Enter coordinates"):
    """
    Prompt the user to input coordinates in the format 'longitude, latitude'.
    Returns a tuple of floats: (lon, lat).
    """
    while True:
        coord_input = input(f"{prompt} (format: lon, lat): ").strip()
        try:
            lon_str, lat_str = coord_input.split(",")
            lon = float(lon_str.strip())
            lat = float(lat_str.strip())
            if not (-180 <= lon <= 180 and -90 <= lat <= 90):
                raise ValueError
            return (lon, lat)
        except ValueError:
            print("Invalid format or out of bounds. Please enter as 'longitude, latitude'.")

def get_country(prompt="Enter country name"):
    """
    Prompt the user to input a country name and normalize it.
    Returns the normalized country name if valid.
    """
    while True:
        country_input = input(f"{prompt}: ").strip()
        normed = normalize_name(country_input)
        if normed in VALID_COUNTRIES:
            return normed
        else:
            print("Invalid country name. Please try again.")

def compute_coordinates_distance():
    """
    Compute distance between two sets of coordinates.
    """
    print("\n--- Coordinates to Coordinates Distance ---")
    lon1, lat1 = get_coordinates("Enter first set of coordinates")
    lon2, lat2 = get_coordinates("Enter second set of coordinates")
    distance = geodesic_distance(lon1, lat1, lon2, lat2)
    print(f"Distance between points: {distance:.2f} km\n")

def compute_coordinates_to_country_distance(country_geoms):
    """
    Compute distance between a coordinate and a country.
    """
    print("\n--- Coordinates to Country Distance ---")
    lon, lat = get_coordinates("Enter coordinates")
    country = get_country("Enter country name")
    geom = country_geoms.get(country)
    if not geom:
        print(f"Geometry for '{country}' not found.\n")
        return
    distance = distance_point_to_polygon(lon, lat, geom, get_sample_size_override(country))
    if distance is not None:
        print(f"Distance from point to '{country}': {distance:.2f} km\n")
    else:
        print(f"Could not compute distance to '{country}'.\n")

def compute_coordinates_to_all_countries_distance(country_geoms):
    """
    Compute distances between a coordinate and all countries.
    Save the results to a JSON file.
    """
    print("\n--- Coordinates to All Countries Distance ---")
    lon, lat = get_coordinates("Enter coordinates")
    output_filename = input("Enter a name for the output JSON file (without extension): ").strip()
    output_path = os.path.join(BASE_DIR, f"{output_filename}.json")

    distances = {}
    for country, geom in country_geoms.items():
        if geom:
            distance = distance_point_to_polygon(lon, lat, geom, get_sample_size_override(country))
            distances[country] = round(distance, 1) if distance is not None else None
        elif country in MICRONATION_COORDS:
            c_lon, c_lat = MICRONATION_COORDS[country]
            distance = geodesic_distance(lon, lat, c_lon, c_lat)
            distances[country] = round(distance, 1)
        else:
            distances[country] = None

    # Save to JSON
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(distances, f, indent=2)
    print(f"Distances saved to '{output_path}'.\n")

def compute_country_to_country_distance(country_geoms):
    """
    Compute distance between two countries.
    """
    print("\n--- Country to Country Distance ---")
    country1 = get_country("Enter the first country name")
    country2 = get_country("Enter the second country name")

    if country1 == country2:
        print("Both countries are the same. Distance is 0 km.\n")
        return

    geom1 = country_geoms.get(country1)
    geom2 = country_geoms.get(country2)

    if geom1 and geom2:
        distance = distance_polygon_to_polygon(geom1, geom2, get_sample_size_override(country1))
        distance = round(distance, 1) if distance is not None else None
    elif country1 in MICRONATION_COORDS and country2 in MICRONATION_COORDS:
        lon1, lat1 = MICRONATION_COORDS[country1]
        lon2, lat2 = MICRONATION_COORDS[country2]
        distance = geodesic_distance(lon1, lat1, lon2, lat2)
        distance = round(distance, 1)
    elif geom1 and country2 in MICRONATION_COORDS:
        lon2, lat2 = MICRONATION_COORDS[country2]
        distance = distance_point_to_polygon(lon2, lat2, geom1, get_sample_size_override(country1))
        distance = round(distance, 1) if distance is not None else None
    elif geom2 and country1 in MICRONATION_COORDS:
        lon1, lat1 = MICRONATION_COORDS[country1]
        distance = distance_point_to_polygon(lon1, lat1, geom2, get_sample_size_override(country2))
        distance = round(distance, 1) if distance is not None else None
    else:
        distance = None

    if distance is not None:
        print(f"Distance between '{country1}' and '{country2}': {distance:.2f} km\n")
    else:
        print(f"Could not compute distance between '{country1}' and '{country2}'.\n")

def compute_country_to_all_countries_distance(country_geoms):
    """
    Compute distances between a country and all other countries.
    Save the results to a JSON file.
    """
    print("\n--- Country to All Countries Distance ---")
    country = get_country("Enter country name")
    output_filename = input("Enter a name for the output JSON file (without extension): ").strip()
    output_path = os.path.join(BASE_DIR, f"{output_filename}.json")

    geom = country_geoms.get(country)
    distances = {}

    for target_country, target_geom in country_geoms.items():
        if target_country == country:
            distances[target_country] = 0.0
            continue

        if geom and target_geom:
            distance = distance_polygon_to_polygon(geom, target_geom, get_sample_size_override(country))
            distances[target_country] = round(distance, 1) if distance is not None else None
        elif geom and target_country in MICRONATION_COORDS:
            lon, lat = MICRONATION_COORDS[target_country]
            distance = distance_point_to_polygon(lon, lat, geom, get_sample_size_override(country))
            distances[target_country] = round(distance, 1) if distance is not None else None
        elif country in MICRONATION_COORDS and target_geom:
            lon, lat = MICRONATION_COORDS[country]
            distance = distance_point_to_polygon(lon, lat, target_geom, get_sample_size_override(target_country))
            distances[target_country] = round(distance, 1) if distance is not None else None
        elif country in MICRONATION_COORDS and target_country in MICRONATION_COORDS:
            lon1, lat1 = MICRONATION_COORDS[country]
            lon2, lat2 = MICRONATION_COORDS[target_country]
            distance = geodesic_distance(lon1, lat1, lon2, lat2)
            distances[target_country] = round(distance, 1)
        else:
            distances[target_country] = None

    # Save to JSON
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(distances, f, indent=2)
    print(f"Distances saved to '{output_path}'.\n")

def get_sample_size_override(country):
    """
    Get the sample size for a given country, considering any overrides.
    """
    return SAMPLE_SIZE_OVERRIDE.get(country, get_sample_size(country))

def main():
    # Load country geometries
    country_geoms = load_country_geometries()

    while True:
        choice = get_user_choice()

        if choice == '1':
            compute_coordinates_distance()
        elif choice == '2':
            compute_coordinates_to_country_distance(country_geoms)
        elif choice == '3':
            compute_coordinates_to_all_countries_distance(country_geoms)
        elif choice == '4':
            compute_country_to_country_distance(country_geoms)
        elif choice == '5':
            compute_country_to_all_countries_distance(country_geoms)
        else:
            print("Invalid choice. Please enter a number between 1 and 5.\n")
            continue

        # Ask if the user wants to perform another operation
        again = input("Do you want to perform another operation? (y/n): ").strip().lower()
        if again != 'y':
            print("Exiting the program. Goodbye!")
            break
        print("\n----------------------------------------\n")

if __name__ == "__main__":
    main()
