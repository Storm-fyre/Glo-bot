import json
import time
import os
import math
import difflib
import geopandas as gpd
from shapely.geometry import MultiPolygon, Polygon, Point
from shapely.ops import unary_union
from pyproj import Geod

##############################################################################
# 1. CONFIG & GLOBALS
##############################################################################

# ------------------- 1.1 The 196 Countries -------------------
VALID_COUNTRIES = {
    # 193 UN members + Taiwan + Vatican City + Palestine
    'afghanistan', 'albania', 'algeria', 'andorra', 'angola', 'antigua and barbuda', 'argentina',
    'armenia', 'australia', 'austria', 'azerbaijan', 'bahamas', 'bahrain', 'bangladesh',
    'barbados', 'belarus', 'belgium', 'belize', 'benin', 'bhutan', 'bolivia',
    'bosnia and herzegovina', 'botswana', 'brazil', 'brunei', 'bulgaria', 'burkina faso',
    'burundi', 'cambodia', 'cameroon', 'canada', 'cape verde', 'central african republic',
    'chad', 'chile', 'china', 'colombia', 'comoros', 'congo', 'costa rica', 'croatia',
    'cuba', 'cyprus', 'czech republic', 'democratic republic of the congo', 'denmark',
    'djibouti', 'dominica', 'dominican republic', 'ecuador', 'egypt', 'el salvador',
    'equatorial guinea', 'eritrea', 'estonia', 'eswatini', 'ethiopia', 'fiji', 'finland',
    'france', 'gabon', 'gambia', 'georgia', 'germany', 'ghana', 'greece', 'grenada',
    'guatemala', 'guinea', 'guinea-bissau', 'guyana', 'haiti', 'honduras', 'hungary',
    'iceland', 'india', 'indonesia', 'iran', 'iraq', 'ireland', 'israel', 'italy',
    'ivory coast', 'jamaica', 'japan', 'jordan', 'kazakhstan', 'kenya', 'kiribati',
    'kuwait', 'kyrgyzstan', 'laos', 'latvia', 'lebanon', 'lesotho', 'liberia', 'libya',
    'liechtenstein', 'lithuania', 'luxembourg', 'madagascar', 'malawi', 'malaysia',
    'maldives', 'mali', 'malta', 'marshall islands', 'mauritania', 'mauritius', 'mexico',
    'micronesia', 'moldova', 'monaco', 'mongolia', 'montenegro', 'morocco', 'mozambique',
    'myanmar', 'namibia', 'nauru', 'nepal', 'netherlands', 'new zealand', 'nicaragua',
    'niger', 'nigeria', 'north korea', 'north macedonia', 'norway', 'oman', 'pakistan',
    'palau', 'palestine', 'panama', 'papua new guinea', 'paraguay', 'peru', 'philippines',
    'poland', 'portugal', 'qatar', 'romania', 'russia', 'rwanda', 'saint kitts and nevis',
    'saint lucia', 'saint vincent and the grenadines', 'samoa', 'san marino',
    'sao tome and principe', 'saudi arabia', 'senegal', 'serbia', 'seychelles',
    'sierra leone', 'singapore', 'slovakia', 'slovenia', 'solomon islands', 'somalia',
    'south africa', 'south korea', 'south sudan', 'spain', 'sri lanka', 'sudan', 'suriname',
    'sweden', 'switzerland', 'syria', 'taiwan', 'tajikistan', 'tanzania', 'thailand',
    'timor-leste', 'togo', 'tonga', 'trinidad and tobago', 'tunisia', 'turkey',
    'turkmenistan', 'tuvalu', 'uganda', 'ukraine', 'united arab emirates',
    'united kingdom', 'united states', 'uruguay', 'uzbekistan', 'vanuatu',
    'vatican city', 'venezuela', 'vietnam', 'yemen', 'zambia', 'zimbabwe'
}

# ------------------- 1.2 Synonym/Alternate Name Map -------------------
SYNONYM_MAP = {
    'republic of korea': 'south korea',
    'korea, republic of': 'south korea',
    'korea (rep.)': 'south korea',
    'korea (south)': 'south korea',

    'korea, dem. people’s rep. (north korea)': 'north korea',
    'democratic people\'s republic of korea': 'north korea',
    'korea, dem. rep. of': 'north korea',

    'ivory coast': 'ivory coast',
    'côte d’ivoire': 'ivory coast',
    'cote d’ivoire': 'ivory coast',

    'bolivia (plurinational state of)': 'bolivia',
    'lao people’s democratic republic': 'laos',
    'lao pdr': 'laos',
    'brunei darussalam': 'brunei',
    'union of the comoros': 'comoros',
    'tibet': 'china',
    'hong kong': 'china',

    'syrian arab republic': 'syria',
    'russian federation': 'russia',
    'viet nam': 'vietnam',
    'timor leste': 'timor-leste',
    'swaziland': 'eswatini',

    'bosnia & herzegovina': 'bosnia and herzegovina',
    'bosnia-herzegovina': 'bosnia and herzegovina',

    'united states of america': 'united states',
    'u.s.a.': 'united states',
    'united states (usa)': 'united states',

    'republic of the congo': 'congo',
    'congo, rep. of the': 'congo',
    'republic of congo': 'congo',

    'democratic republic of congo': 'democratic republic of the congo',
    'congo, the democratic republic of': 'democratic republic of the congo',
    'drc': 'democratic republic of the congo',

    'iran (islamic republic of)': 'iran',
    'iran, islamic republic of': 'iran',

    'gambia, the': 'gambia',
    'the gambia': 'gambia',

    'moldova (republic of)': 'moldova',
    'republic of moldova': 'moldova',

    'tanzania, united republic of': 'tanzania',
    'united republic of tanzania': 'tanzania',

    'venezuela (bolivarian republic of)': 'venezuela',
    'venezuela, bolivarian republic of': 'venezuela',

    'south sudan (republic of)': 'south sudan',
}

# ------------------- 1.3 Sample Size Mapping -------------------
SAMPLE_SIZE_MAP = {
    'russia': 1700,
    'canada': 900,
    'china': 900,
    'united states': 900,
    'brazil': 900,
    'australia': 900,
    'india': 900,
}

DEFAULT_SAMPLE_SIZE = 300

def get_sample_size(country: str, custom_sample: int = None) -> int:
    """
    Returns the sample size based on the country.
    If custom_sample is provided, it overrides the default/sample size.
    """
    if custom_sample is not None:
        return custom_sample
    return SAMPLE_SIZE_MAP.get(country, DEFAULT_SAMPLE_SIZE)

# ------------------- 1.4 Mainland (Polygon) Filtering Rules -------------------
POLYGON_SELECTION_RULES = {
    'angola':            {'count': 2},
    'azerbaijan':        {'count': 2},
    'chile':             {'distance_threshold': 500},   # exclude Easter Island
    'china':             {'count': 2},                  # Mainland + Hainan (exclude Taiwan)
    'denmark':           {'count': 3},                  # exclude Faroe, Greenland
    'equatorial guinea': {'count': 2},
    'france':            {'count': 1},                  # exclude Corsica
    'greece':            {'count': 3},                  # Mainland, Crete, Lesbos
    'indonesia':         {'all': True},
    'italy':             {'count': 3},                  # Mainland, Sicily, Sardinia
    'japan':             {'custom_japan': True},        # special logic for Tsushima
    'malaysia':          {'count': 2},
    'morocco':           {'count': 1},                  # exclude Western Sahara
    'netherlands':       {'count': 2},
    'new zealand':       {'count': 3},                  # North, South, Stewart
    'norway':            {'count': 2},
    'oman':              {'count': 2},
    'papua new guinea':  {'count': 5},
    'philippines':       {'all': True},
    'portugal':          {'count': 1},   # exclude Azores, Madeira
    'russia':            {'all': True},  # but filter out Kaliningrad below
    'south korea':       {'count': 2},
    'spain':             {'count': 2},   # Mainland + Balearics; exclude Canary
    'timor-leste':       {'count': 2},
    'turkey':            {'count': 1},
    'united kingdom':    {'count': 3},   # GB, Northern Ireland, Lewis/Harris
    'united states':     {'all': True},  # but filter out non-states below
    'venezuela':         {'count': 1},
    'yemen':             {'count': 2},
}

# ------------------- 1.5 Micronation Coordinates (no polygons) ---------------
MICRONATION_COORDS = {
    "andorra": (1.6016, 42.5424),
    "antigua and barbuda": (-61.8456, 17.0747),
    "bahrain": (50.5577, 26.0667),
    "barbados": (-59.5432, 13.1939),
    "cape verde": (-23.5087, 14.9300),
    "comoros": (43.3333, -11.6455),
    "dominica": (-61.3710, 15.4239),
    "grenada": (-61.6792, 12.1165),
    "kiribati": (173.0314, 1.3382),
    "liechtenstein": (9.5537, 47.1660),
    "maldives": (73.5093, 4.1755),
    "malta": (14.3754, 35.9375),
    "marshall islands": (171.1854, 7.1315),
    "mauritius": (57.5522, -20.3484),
    "micronesia": (158.2239, 6.9248),
    "monaco": (7.4128, 43.7306),
    "nauru": (166.9315, -0.5338),
    "palau": (134.4795, 7.3419),
    "saint kitts and nevis": (-62.7830, 17.3578),
    "saint lucia": (-60.9789, 13.9094),
    "saint vincent and the grenadines": (-61.2872, 13.2528),
    "samoa": (-172.1046, -13.7590),
    "san marino": (12.4578, 43.9424),
    "sao tome and principe": (6.7273, 0.3302),
    "seychelles": (55.4915, -4.6796),
    "singapore": (103.8198, 1.3521),
    "tonga": (-175.1982, -21.1789),
    "tuvalu": (179.2168, -8.5199),
    "vatican city": (12.4534, 41.9029)
}

# ------------------- 1.6 Boundary Sample Settings -------------------
# Removed the global BOUNDARY_SAMPLES as we are using dynamic sample sizes.

##############################################################################
# 2. HELPER FUNCTIONS
##############################################################################

# A global geodesic object (WGS84 ellipsoid)
geod = Geod(ellps='WGS84')

def geodesic_distance(lon1, lat1, lon2, lat2):
    """
    Return geodesic distance (km) between two (lon, lat) points on WGS84.
    """
    _, _, dist_m = geod.inv(lon1, lat1, lon2, lat2)
    return dist_m / 1000.0

def distance_polygon_to_polygon(polyA, polyB, samples):
    """
    Approx. minimal geodesic distance between boundaries of two polygons
    by sampling their boundaries. Returns distance in km or None if empty.
    """
    if not polyA or polyA.is_empty or not polyB or polyB.is_empty:
        return None

    boundaryA = polyA.boundary
    boundaryB = polyB.boundary

    lenA = boundaryA.length
    lenB = boundaryB.length

    stepA = lenA / samples
    stepB = lenB / samples

    min_dist = float('inf')

    # Sample boundary of polyA
    ptsA = []
    for i in range(samples + 1):
        dA = i * stepA
        ptA = boundaryA.interpolate(dA)
        ptsA.append((ptA.x, ptA.y))

    # Sample boundary of polyB
    ptsB = []
    for j in range(samples + 1):
        dB = j * stepB
        ptB = boundaryB.interpolate(dB)
        ptsB.append((ptB.x, ptB.y))

    # Compare each point in A to each in B
    for (xA, yA) in ptsA:
        for (xB, yB) in ptsB:
            d_km = geodesic_distance(xA, yA, xB, yB)
            if d_km < min_dist:
                min_dist = d_km

    return min_dist if min_dist != float('inf') else None

def centroid_distance(polyA, polyB):
    """
    Return geodesic distance (km) between polygon centroids.
    """
    cA = polyA.centroid
    cB = polyB.centroid
    return geodesic_distance(cA.x, cA.y, cB.x, cB.y)

def normalize_name(raw_name: str) -> str:
    """
    Convert various official or alternate forms (like "Republic of Korea")
    to our canonical forms (like "south korea"), using the synonyms above.
    """
    import re

    candidate = raw_name.lower().strip()

    # Direct check
    if candidate in VALID_COUNTRIES:
        return candidate

    # Synonym check
    if candidate in SYNONYM_MAP:
        mapped = SYNONYM_MAP[candidate]
        if mapped in VALID_COUNTRIES:
            return mapped

    # Extra cleaning of punctuation (parentheses, commas, etc.)
    fallback = re.sub(r'[\(\),\'’\.]', '', candidate)
    fallback = fallback.replace("  ", " ").strip()

    if fallback in SYNONYM_MAP:
        mapped = SYNONYM_MAP[fallback]
        if mapped in VALID_COUNTRIES:
            return mapped
    if fallback in VALID_COUNTRIES:
        return fallback

    return ""  # Not recognized

def distance_point_to_point(lon1, lat1, lon2, lat2):
    """
    Geodesic distance between two (lon, lat) points.
    """
    return geodesic_distance(lon1, lat1, lon2, lat2)

def distance_point_to_polygon(point_lon, point_lat, poly, samples):
    """
    Approx. minimal geodesic distance between a single point and a polygon
    by sampling the polygon boundary.
    """
    if not poly or poly.is_empty:
        return None

    boundary = poly.boundary
    length = boundary.length
    step = length / samples

    min_dist = float('inf')
    for i in range(samples + 1):
        d = i * step
        pt = boundary.interpolate(d)
        d_km = geodesic_distance(pt.x, pt.y, point_lon, point_lat)
        if d_km < min_dist:
            min_dist = d_km

    return min_dist if min_dist != float('inf') else None

##############################################################################
# 3. MAINLAND FILTERING LOGIC
##############################################################################

from shapely.geometry import MultiPolygon, Polygon, Point
from shapely.ops import unary_union

def filter_polygons(multi_or_poly, country_name):
    """
    Returns a (Multi)Polygon for the “mainland” of a country based on your rules,
    including special exclusions for Russia, US, Japan, etc.
    """
    # If geometry is a single Polygon, wrap it in a list
    if isinstance(multi_or_poly, Polygon):
        polygons = [multi_or_poly]
    # If geometry is a MultiPolygon, get its .geoms list
    elif isinstance(multi_or_poly, MultiPolygon):
        polygons = list(multi_or_poly.geoms)
    else:
        # e.g. None or invalid geometry
        return None

    # Remove empty
    polygons = [p for p in polygons if p and not p.is_empty]

    # 1) Modifications for Russia
    if country_name == 'russia':

        # Exclude Kaliningrad
        kal_min_lon, kal_max_lon = (19.0, 23.0)
        kal_min_lat, kal_max_lat = (54.0, 55.5)

        # Exclude Crimea
        cr_min_lon, cr_max_lon = (32.0, 36.5)
        cr_min_lat, cr_max_lat = (44.0, 46.5)

        def is_kaliningrad(p):
            c = p.centroid
            return (kal_min_lon <= c.x <= kal_max_lon) and (kal_min_lat <= c.y <= kal_max_lat)

        def is_crimea(p):
            c = p.centroid
            return (cr_min_lon <= c.x <= cr_max_lon) and (cr_min_lat <= c.y <= cr_max_lat)

        polygons = [p for p in polygons if not is_kaliningrad(p) and not is_crimea(p)]

    # 2) US => only 50 states (exclude territories)
    if country_name == 'united states':
        def is_in_50_states(p):
            c = p.centroid
            lat, lon = c.y, c.x
            # Mainland bounding box
            in_mainland = (-125.0 <= lon <= -66.0) and (24.4 <= lat <= 49.5)
            # Alaska bounding box
            in_alaska = (-172.0 <= lon <= -130.0) and (51.2 <= lat <= 72.0)
            # Hawaii bounding box
            in_hawaii = (-161.0 <= lon <= -154.0) and (18.5 <= lat <= 23.0)
            return in_mainland or in_alaska or in_hawaii

        polygons = [p for p in polygons if is_in_50_states(p)]

    # 3) Japan => 4 main islands + Tsushima coordinate
    if country_name == 'japan':
        # Step 1: exclude polygons south of 30N (Okinawa, etc.)
        polygons = [p for p in polygons if p.centroid.y >= 30.0]

        # Step 2: sort polygons by area and select the top 4
        polygons_sorted = sorted(polygons, key=lambda p: p.area, reverse=True)
        main_4 = polygons_sorted[:4]

        # Step 3: Create a small buffer polygon around (129.3, 34.4) to represent Tsushima
        tsushima_point = Point(129.3, 34.4)
        # Buffer size is in degrees; adjust as needed. Here, 0.05 degrees ~5.5 km
        tsushima_buffer = tsushima_point.buffer(0.05)

        # Combine the four main islands with the Tsushima buffer
        out_polys = main_4 + [tsushima_buffer]

        if out_polys:
            return MultiPolygon(out_polys)
        else:
            return None

    # 4) Standard selection
    polygons_sorted = sorted(polygons, key=lambda p: p.area, reverse=True)
    if not polygons_sorted:
        return None

    rule = POLYGON_SELECTION_RULES.get(country_name, {'count': 1})

    if 'all' in rule and rule['all'] == True:
        return MultiPolygon(polygons_sorted)

    if 'count' in rule:
        n = rule['count']
        return MultiPolygon(polygons_sorted[:n])

    if 'distance_threshold' in rule:
        dist_thr = rule['distance_threshold']
        main_poly = polygons_sorted[0]
        keep = [main_poly]
        for p in polygons_sorted[1:]:
            if centroid_distance(main_poly, p) <= dist_thr:
                keep.append(p)
        return MultiPolygon(keep)

    # Default => largest polygon
    return polygons_sorted[0]

##############################################################################
# 4. MAIN SCRIPT WITH USER INTERFACE
##############################################################################

def load_shapefile(shapefile_path):
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

    return gdf, found_cols

def extract_valid_country(row_series, found_cols):
    """
    Try each column in found_cols; return first valid country name or "".
    row_series is a pandas Series, so we do row_series[col].
    """
    for col in found_cols:
        raw_val = str(row_series[col]).strip()
        normed = normalize_name(raw_val)
        if normed in VALID_COUNTRIES:
            return normed
    return ""

def build_country_geometries(gdf, found_cols):
    """
    Builds a dictionary mapping each valid country to its filtered geometry.
    """
    total_rows = len(gdf)
    print(f"Found {total_rows} rows in the shapefile. Processing...")

    country_geoms = {}
    for idx, (index, row) in enumerate(gdf.iterrows(), start=1):
        print(f"[Polygon Extraction] {idx}/{total_rows} => ", end="")
        country_name = extract_valid_country(row, found_cols)
        if not country_name:
            print("Skipped (not recognized).")
            continue

        print(f"Recognized as '{country_name}'.")
        raw_geom = row.geometry
        filtered = filter_polygons(raw_geom, country_name)
        if filtered and not filtered.is_empty:
            # If repeated: union
            if country_name not in country_geoms:
                country_geoms[country_name] = filtered
            else:
                combined = unary_union([country_geoms[country_name], filtered])
                country_geoms[country_name] = combined

    # We want all 196 countries. If no polygon => None (we might do micronation logic).
    final_polygons = {}
    for ctry in VALID_COUNTRIES:
        if ctry in country_geoms:
            final_polygons[ctry] = country_geoms[ctry]
        else:
            final_polygons[ctry] = None

    return final_polygons

def compute_distance_map(final_polygons, option, **kwargs):
    """
    Computes distances based on the selected option.
    """
    distance_map = {}

    if option == 1:
        # Coordinates to Coordinates
        lon1, lat1, lon2, lat2 = kwargs.get('coords1'), kwargs.get('coords2')
        distance = distance_point_to_point(lon1, lat1, lon2, lat2)
        return distance

    elif option == 2:
        # Coordinates to Single Country
        lon, lat = kwargs.get('coords')
        country = kwargs.get('country')
        sample_size = kwargs.get('sample_size')

        poly = final_polygons.get(country, None)
        if poly:
            dist_km = distance_point_to_polygon(lon, lat, poly, samples=sample_size)
            return round(dist_km, 1) if dist_km is not None else "No valid distance"
        else:
            # If the country is a micronation with coordinates
            if country in MICRONATION_COORDS:
                c_lon, c_lat = MICRONATION_COORDS[country]
                return round(distance_point_to_point(lon, lat, c_lon, c_lat), 1)
            else:
                return "Country geometry not available."

    elif option == 3:
        # Coordinates to All Countries
        lon, lat = kwargs.get('coords')
        output_filename = kwargs.get('output_filename')
        sample_size = kwargs.get('sample_size')

        distance_map = {}
        for ctry in sorted(VALID_COUNTRIES):
            if ctry in MICRONATION_COORDS:
                c_lon, c_lat = MICRONATION_COORDS[ctry]
                dist = distance_point_to_point(lon, lat, c_lon, c_lat)
                distance_map[ctry] = round(dist, 1)
            else:
                poly = final_polygons.get(ctry, None)
                if poly:
                    dist = distance_point_to_polygon(lon, lat, poly, samples=sample_size)
                    distance_map[ctry] = round(dist, 1) if dist is not None else 9999999
                else:
                    distance_map[ctry] = "No geometry available."

        # Save to JSON
        with open(output_filename, "w", encoding="utf-8") as f:
            json.dump(distance_map, f, indent=2)

        return distance_map

    elif option == 4:
        # Single Country to Single Country
        country1 = kwargs.get('country1')
        country2 = kwargs.get('country2')
        sample_size1 = kwargs.get('sample_size1')
        sample_size2 = kwargs.get('sample_size2')

        if country1 == country2:
            return 0.0

        poly1 = final_polygons.get(country1, None)
        poly2 = final_polygons.get(country2, None)

        if poly1 and poly2:
            dist_km = distance_polygon_to_polygon(poly1, poly2, samples=max(sample_size1, sample_size2))
            return round(dist_km, 1) if dist_km is not None else "No valid distance"
        elif (country1 in MICRONATION_COORDS) and (country2 in MICRONATION_COORDS):
            lon1, lat1 = MICRONATION_COORDS[country1]
            lon2, lat2 = MICRONATION_COORDS[country2]
            dist = distance_point_to_point(lon1, lat1, lon2, lat2)
            return round(dist, 1)
        elif poly1 and (country2 in MICRONATION_COORDS):
            lon2, lat2 = MICRONATION_COORDS[country2]
            dist_km = distance_point_to_polygon(lon2, lat2, poly1, samples=sample_size1)
            return round(dist_km, 1) if dist_km is not None else "No valid distance"
        elif poly2 and (country1 in MICRONATION_COORDS):
            lon1, lat1 = MICRONATION_COORDS[country1]
            dist_km = distance_point_to_polygon(lon1, lat1, poly2, samples=sample_size2)
            return round(dist_km, 1) if dist_km is not None else "No valid distance"
        else:
            return "Distance cannot be computed due to missing geometries."

    elif option == 5:
        # Single Country to All Countries
        country = kwargs.get('country')
        sample_size = kwargs.get('sample_size')
        output_filename = kwargs.get('output_filename')

        distance_map = {}
        poly1 = final_polygons.get(country, None)

        for c2 in sorted(VALID_COUNTRIES):
            if c2 == country:
                distance_map[c2] = 0.0
                continue

            if c2 in MICRONATION_COORDS:
                if country in MICRONATION_COORDS:
                    # Both are micronations
                    lon1, lat1 = MICRONATION_COORDS[country]
                    lon2, lat2 = MICRONATION_COORDS[c2]
                    dist = distance_point_to_point(lon1, lat1, lon2, lat2)
                    distance_map[c2] = round(dist, 1)
                else:
                    # Country to micronation
                    lon2, lat2 = MICRONATION_COORDS[c2]
                    if poly1:
                        dist_km = distance_point_to_polygon(lon2, lat2, poly1, samples=sample_size)
                        distance_map[c2] = round(dist_km, 1) if dist_km is not None else 9999999
                    else:
                        distance_map[c2] = "No geometry available."
            else:
                if poly1 and final_polygons.get(c2, None):
                    dist_km = distance_polygon_to_polygon(poly1, final_polygons[c2], samples=sample_size)
                    distance_map[c2] = round(dist_km, 1) if dist_km is not None else 9999999
                elif poly1 and (c2 in MICRONATION_COORDS):
                    lon2, lat2 = MICRONATION_COORDS[c2]
                    dist_km = distance_point_to_polygon(lon2, lat2, poly1, samples=sample_size)
                    distance_map[c2] = round(dist_km, 1) if dist_km is not None else 9999999
                elif (country in MICRONATION_COORDS) and (c2 not in MICRONATION_COORDS):
                    lon1, lat1 = MICRONATION_COORDS[country]
                    poly2 = final_polygons.get(c2, None)
                    if poly2:
                        dist_km = distance_point_to_polygon(lon1, lat1, poly2, samples=sample_size)
                        distance_map[c2] = round(dist_km, 1) if dist_km is not None else 9999999
                    else:
                        distance_map[c2] = "No geometry available."
                else:
                    distance_map[c2] = "No valid distance."

        # Save to JSON
        with open(output_filename, "w", encoding="utf-8") as f:
            json.dump(distance_map, f, indent=2)

        return distance_map

    else:
        return "Invalid option."

def main():
    shapefile_path = r"C:\Users\raksh\My Drive\Tech\deploy\ne_110m_admin_0_countries.shp"

    # Load shapefile and build geometries
    try:
        gdf, found_cols = load_shapefile(shapefile_path)
        final_polygons = build_country_geometries(gdf, found_cols)
    except Exception as e:
        print(f"Error loading shapefile: {e}")
        return

    # User Interface Loop
    while True:
        print("\nDistance Calculation Options:")
        print("1. Coordinates to Coordinates")
        print("2. Coordinates to Single Country")
        print("3. Coordinates to All Countries")
        print("4. Single Country to Single Country")
        print("5. Single Country to All Countries")
        print("6. Exit")

        try:
            choice = int(input("Select an option (1-6): ").strip())
        except ValueError:
            print("Invalid input. Please enter a number between 1 and 6.")
            continue

        if choice == 6:
            print("Exiting program.")
            break

        if choice == 1:
            # Coordinates to Coordinates
            try:
                lon1 = float(input("Enter Longitude 1: ").strip())
                lat1 = float(input("Enter Latitude 1: ").strip())
                lon2 = float(input("Enter Longitude 2: ").strip())
                lat2 = float(input("Enter Latitude 2: ").strip())
                distance = compute_distance_map(
                    final_polygons,
                    option=1,
                    coords1=(lon1, lat1),
                    coords2=(lon2, lat2)
                )
                print(f"Distance: {distance} km")
            except ValueError:
                print("Invalid coordinates entered.")

        elif choice == 2:
            # Coordinates to Single Country
            try:
                lon = float(input("Enter Longitude: ").strip())
                lat = float(input("Enter Latitude: ").strip())
                country_input = input("Enter the target country: ").strip().lower()
                country = normalize_name(country_input)
                if not country:
                    print("Country not recognized.")
                    continue

                sample_size = input(f"Enter sample size for '{country}' (Press Enter for default {get_sample_size(country)}): ").strip()
                sample_size = int(sample_size) if sample_size else None
                sample_size = get_sample_size(country, custom_sample=sample_size)

                distance = compute_distance_map(
                    final_polygons,
                    option=2,
                    coords=(lon, lat),
                    country=country,
                    sample_size=sample_size
                )
                print(f"Distance to '{country}': {distance} km")
            except ValueError:
                print("Invalid input.")

        elif choice == 3:
            # Coordinates to All Countries
            try:
                lon = float(input("Enter Longitude: ").strip())
                lat = float(input("Enter Latitude: ").strip())
                sample_size = input(f"Enter sample size for distance calculations (Press Enter for default {DEFAULT_SAMPLE_SIZE}): ").strip()
                sample_size = int(sample_size) if sample_size else DEFAULT_SAMPLE_SIZE

                output_filename = input("Enter the output JSON filename (e.g., distances.json): ").strip()
                if not output_filename.endswith('.json'):
                    output_filename += '.json'

                distance_map = compute_distance_map(
                    final_polygons,
                    option=3,
                    coords=(lon, lat),
                    output_filename=output_filename,
                    sample_size=sample_size
                )
                print(f"Distances saved to '{output_filename}'.")
            except ValueError:
                print("Invalid input.")

        elif choice == 4:
            # Single Country to Single Country
            try:
                country1_input = input("Enter the first country: ").strip().lower()
                country1 = normalize_name(country1_input)
                if not country1:
                    print("First country not recognized.")
                    continue

                country2_input = input("Enter the second country: ").strip().lower()
                country2 = normalize_name(country2_input)
                if not country2:
                    print("Second country not recognized.")
                    continue

                sample_size1 = input(f"Enter sample size for '{country1}' (Press Enter for default {get_sample_size(country1)}): ").strip()
                sample_size1 = int(sample_size1) if sample_size1 else None
                sample_size1 = get_sample_size(country1, custom_sample=sample_size1)

                sample_size2 = input(f"Enter sample size for '{country2}' (Press Enter for default {get_sample_size(country2)}): ").strip()
                sample_size2 = int(sample_size2) if sample_size2 else None
                sample_size2 = get_sample_size(country2, custom_sample=sample_size2)

                distance = compute_distance_map(
                    final_polygons,
                    option=4,
                    country1=country1,
                    country2=country2,
                    sample_size1=sample_size1,
                    sample_size2=sample_size2
                )
                print(f"Distance between '{country1}' and '{country2}': {distance} km")
            except ValueError:
                print("Invalid input.")

        elif choice == 5:
            # Single Country to All Countries
            try:
                country_input = input("Enter the source country: ").strip().lower()
                country = normalize_name(country_input)
                if not country:
                    print("Country not recognized.")
                    continue

                sample_size = input(f"Enter sample size for '{country}' (Press Enter for default {get_sample_size(country)}): ").strip()
                sample_size = int(sample_size) if sample_size else None
                sample_size = get_sample_size(country, custom_sample=sample_size)

                output_filename = input("Enter the output JSON filename (e.g., distances.json): ").strip()
                if not output_filename.endswith('.json'):
                    output_filename += '.json'

                distance_map = compute_distance_map(
                    final_polygons,
                    option=5,
                    country=country,
                    sample_size=sample_size,
                    output_filename=output_filename
                )
                print(f"Distances saved to '{output_filename}'.")
            except ValueError:
                print("Invalid input.")

        else:
            print("Invalid option. Please select a number between 1 and 6.")

    print("Program terminated.")

if __name__ == "__main__":
    main()
