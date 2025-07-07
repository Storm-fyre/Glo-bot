import os
import json
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon

def normalize_name(raw_name: str) -> str:
    """
    Normalize country names by converting to lowercase and handling common synonyms.
    """
    name = raw_name.lower().strip()
    synonym_map = {
        'united states of america': 'united states',
        'democratic republic of the congo': 'democratic republic of the congo',
        'republic of the congo': 'republic of the congo',
        'côte d’ivoire': 'ivory coast',
        'cote d’ivoire': 'ivory coast',
        'ivory coast': 'ivory coast',
        'palestine': 'palestine',
        # Add more synonyms or normalization rules as needed
    }
    return synonym_map.get(name, name)

def main():
    # Define the directory and shapefile name
    directory = r"C:\Users\raksh\My Drive\Tech\deploy\distance"
    shapefile_name = "ne_110m_admin_0_countries.shp"
    shapefile_path = os.path.join(directory, shapefile_name)
    
    # Output file name
    output_file = os.path.join(directory, "country_polygon_centroids.json")
    
    # Load the shapefile
    print("Loading shapefile...")
    try:
        gdf = gpd.read_file(shapefile_path)
    except Exception as e:
        print(f"Error loading shapefile: {e}")
        return
    
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
    
    # Build a mapping from country name to list of geometries
    country_geometries = {}
    total_rows = len(gdf)
    print(f"Found {total_rows} rows in the shapefile. Processing...")
    
    for idx, (index, row) in enumerate(gdf.iterrows(), start=1):
        # Extract country name from available columns
        country_name = None
        for col in found_cols:
            raw_name = str(row[col]).strip()
            if raw_name and raw_name.upper() != "NULL":
                country_name = normalize_name(raw_name)
                break  # Use the first non-empty name found
        
        if not country_name:
            print(f"Row {idx}: Country name not found or is NULL. Skipping.")
            continue
        
        # Get the geometry
        geometry = row.geometry
        if geometry is None or geometry.is_empty:
            print(f"Row {idx}: Geometry is empty for country '{country_name}'. Skipping.")
            continue
        
        # Initialize the list for the country if not already
        if country_name not in country_geometries:
            country_geometries[country_name] = []
        
        # Append the geometry to the country's list
        country_geometries[country_name].append(geometry)
        
        print(f"Row {idx}: Country '{country_name}' processed.")
    
    # Prepare the final mapping
    country_polygon_info = {}
    
    print("\nCalculating polygon counts and centroids...")
    for country, geometries in country_geometries.items():
        total_polygons = 0
        centroids = []
        
        for geom in geometries:
            if isinstance(geom, Polygon):
                total_polygons += 1
                centroid = geom.centroid
                centroids.append({"lon": round(centroid.x, 6), "lat": round(centroid.y, 6)})
            elif isinstance(geom, MultiPolygon):
                num_polys = len(geom.geoms)
                total_polygons += num_polys
                for poly in geom.geoms:
                    centroid = poly.centroid
                    centroids.append({"lon": round(centroid.x, 6), "lat": round(centroid.y, 6)})
            else:
                print(f"Unsupported geometry type for country '{country}': {type(geom)}. Skipping.")
        
        country_polygon_info[country] = {
            "num_polygons": total_polygons,
            "centroids": centroids
        }
        print(f"Country '{country}': {total_polygons} polygon(s), {len(centroids)} centroid(s).")
    
    # Handle countries with no polygons (if any)
    all_countries = set([normalize_name(name) for name in gdf[found_cols[0]].dropna()])
    missing_countries = all_countries - set(country_geometries.keys())
    for missing in missing_countries:
        country_polygon_info[missing] = {
            "num_polygons": 0,
            "centroids": []
        }
        print(f"Country '{missing}' has no polygons.")
    
    # Save the mapping to a JSON file
    try:
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(country_polygon_info, f, indent=2)
        print(f"\nPolygon counts and centroids saved to '{output_file}'.")
    except Exception as e:
        print(f"Error writing to JSON file: {e}")

if __name__ == "__main__":
    main()
