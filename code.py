import json
import os

##############################################################################
# 1) CONFIG
##############################################################################

JSON_FILENAME = "country_distances.json"
DISTANCE_UNKNOWN = 9999999

# List of mischievous countries with their custom tolerance settings
MISCHIEVOUS_COUNTRIES = {
    "canada": {
        "tolerance": 500,
        "distance_threshold": 2500
    },
    # Add more mischievous countries here as needed
}

def get_tolerance(guess_country, dist_value):
    """
    Returns the tolerance based on the guess_country and the reported distance.
    For mischievous countries, applies custom tolerance rules.
    Otherwise, returns the default tolerance.
    """
    if guess_country in MISCHIEVOUS_COUNTRIES:
        config = MISCHIEVOUS_COUNTRIES[guess_country]
        if dist_value > config["distance_threshold"]:
            return config["tolerance"]
    return 50  # Default tolerance

##############################################################################
# 2) LOADING & FILTERING DATA
##############################################################################

def load_distance_data(json_file):
    with open(json_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    return data

def get_supported_countries(distance_map):
    """
    Return countries that do NOT have only unknown distances (9999999).
    """
    supported = []
    for c1 in distance_map:
        distances_c1 = distance_map[c1].values()
        # If absolutely everything is 9999999, skip
        if all(d == DISTANCE_UNKNOWN for d in distances_c1):
            continue
        supported.append(c1)
    return set(supported)

##############################################################################
# 3) REFINEMENT LOGIC
##############################################################################

def match_distance_with_tolerance(actual_dist, reported_dist, guess_country):
    """
    Check if 'actual_dist' is within ±some tolerance of 'reported_dist',
    where that tolerance depends on 'reported_dist' and the guess_country.
    """
    tol = get_tolerance(guess_country, reported_dist)
    low = reported_dist - tol
    high = reported_dist + tol
    return (low <= actual_dist <= high)

def refine_distance_with_tolerance(distance_map, possible_targets, guess_country, reported_dist):
    """
    Filter possible_targets to those for which distance_map[guess_country][country]
    is consistent with reported_dist ± tolerance.
    """
    new_possible = set()
    for c in possible_targets:
        dist_c = distance_map[guess_country].get(c, DISTANCE_UNKNOWN)
        if dist_c == DISTANCE_UNKNOWN:
            # skip
            continue
        if match_distance_with_tolerance(dist_c, reported_dist, guess_country):
            new_possible.add(c)
    return new_possible

def refine_further(distance_map, possible_targets, old_guess, new_guess):
    """
    If the game said "new_guess" is COOLER than "old_guess" (further away from the target),
    that means dist(new_guess, target) > dist(old_guess, target).
    Keep only countries C where distance_map[new_guess][C] > distance_map[old_guess][C].
    """
    new_possible = set()
    for c in possible_targets:
        d_old = distance_map[old_guess].get(c, DISTANCE_UNKNOWN)
        d_new = distance_map[new_guess].get(c, DISTANCE_UNKNOWN)
        if d_old == DISTANCE_UNKNOWN or d_new == DISTANCE_UNKNOWN:
            # skip
            continue
        if d_new > d_old:
            new_possible.add(c)
    return new_possible

def refine_warmer(distance_map, possible_targets, old_guess, new_guess):
    """
    The opposite of 'cooler'. If the game said 'new_guess' is WARMER than 'old_guess',
    that means dist(new_guess, target) < dist(old_guess, target).
    Keep only countries C where distance_map[new_guess][C] < distance_map[old_guess][C].
    """
    new_possible = set()
    for c in possible_targets:
        d_old = distance_map[old_guess].get(c, DISTANCE_UNKNOWN)
        d_new = distance_map[new_guess].get(c, DISTANCE_UNKNOWN)
        if d_old == DISTANCE_UNKNOWN or d_new == DISTANCE_UNKNOWN:
            # skip
            continue
        if d_new < d_old:
            new_possible.add(c)
    return new_possible

def refine_adjacent(distance_map, possible_targets, guess_country, threshold=10):
    """
    If the user types 'adjacent', it means the guessed country shares a border
    with the target. We interpret that as countries whose distance from 'guess_country'
    is below some threshold (default=10 km). Adjust threshold to taste.
    Additionally, remove the guessed country from "Remaining possibilities".
    """
    new_possible = set()
    for c in possible_targets:
        dist_c = distance_map[guess_country].get(c, DISTANCE_UNKNOWN)
        if dist_c == DISTANCE_UNKNOWN:
            continue
        if dist_c < threshold:
            new_possible.add(c)
    # Remove the guessed country itself from possible targets
    new_possible.discard(guess_country)
    return new_possible

##############################################################################
# 4) MAIN SOLVER LOGIC
##############################################################################

def main():
    # Load data
    if not os.path.exists(JSON_FILENAME):
        print(f"Error: {JSON_FILENAME} not found in the current directory.")
        return

    distance_map = load_distance_data(JSON_FILENAME)
    supported = get_supported_countries(distance_map)

    # Start with all supported
    possible_targets = set(supported)

    # A small helper function to pick a "next guess" from possible_targets
    # Mischievous countries are picked last
    def pick_next_guess(possible_set, guesses_so_far):
        non_mischief = [c for c in possible_set if c not in MISCHIEVOUS_COUNTRIES]
        for c in non_mischief:
            if c not in guesses_so_far:
                return c
        # If no non-mischief left, pick from mischief
        for c in possible_set:
            if c not in guesses_so_far and c in MISCHIEVOUS_COUNTRIES:
                return c
        # Fallback if everything is guessed
        return next(iter(possible_set)) if possible_set else None

    # 1) Ask user for the first guess
    first_guess = input("Enter your FIRST guess country (e.g., 'france'): ").strip().lower()
    if first_guess not in supported:
        print(f"'{first_guess}' is not recognized or not supported. Exiting.")
        return

    # 2) Ask user for distance => For the FIRST guess, accept only standalone numerical distance
    first_dist_str = input(
        f"Distance from game for '{first_guess}' (e.g., '1500'): "
    ).strip().lower()

    guesses = []
    last_guess = None

    # Handle first guess input
    try:
        # Attempt to parse as standalone numerical distance
        first_dist = float(first_dist_str)
        # Refine based on distance tolerance
        possible_targets = refine_distance_with_tolerance(distance_map, possible_targets, first_guess, first_dist)
        guesses.append(first_guess)
        last_guess = first_guess
    except ValueError:
        # Invalid input since only numerical distance is allowed for the first guess
        print("Invalid input for the first guess. Must be a numerical distance (e.g., '1500'). Exiting.")
        return

    # Main loop: pick next guess => user inputs distance or 'warmer <distance>' or 'cooler' or 'adjacent' or 'done' => refine
    while True:
        if len(possible_targets) == 1:
            # We have exactly 1 solution
            sol = next(iter(possible_targets))
            print(f"** We have EXACTLY ONE possible target: {sol.upper()}! **")
            break
        elif len(possible_targets) == 0:
            print("No possible targets remain! Possibly contradictory data.")
            break

        # 3) Suggest next guess
        next_guess = pick_next_guess(possible_targets, guesses)
        if not next_guess:
            print("No next guess found. Possibly we've guessed everything. Stopping.")
            break

        print(f"** SUGGESTED NEXT GUESS: {next_guess.upper()} **")
        guesses.append(next_guess)

        # 4) Ask user for distance or 'warmer <distance>' or 'cooler' or 'adjacent' or 'done'
        dist_str = input(
            f"Distance from game for '{next_guess}' (e.g., 'warmer 1500', 'cooler', 'adjacent', 'done'): "
        ).strip().lower()
        if dist_str == "done":
            break

        if dist_str.startswith("warmer"):
            # Expecting 'warmer <distance>'
            parts = dist_str.split()
            if len(parts) != 2:
                print("Invalid input for 'warmer'. Please enter as 'warmer <distance>'. Exiting.")
                break
            try:
                warmer_dist = float(parts[1])
            except ValueError:
                print("Invalid distance value for 'warmer'. Exiting.")
                break
            # Refine based on warmer
            possible_targets = refine_warmer(distance_map, possible_targets, last_guess, next_guess)
            # Refine based on distance tolerance
            possible_targets = refine_distance_with_tolerance(distance_map, possible_targets, next_guess, warmer_dist)
        elif dist_str == "cooler":
            # new_guess is further away from target than last_guess
            possible_targets = refine_further(distance_map, possible_targets, last_guess, next_guess)
        elif dist_str == "adjacent":
            # new_guess shares a border with the target
            possible_targets = refine_adjacent(distance_map, possible_targets, next_guess)
        else:
            # Invalid input since standalone numerical inputs are not allowed after the first guess
            print("Invalid input. Must be 'warmer <distance>', 'cooler', 'adjacent', or 'done'. Exiting.")
            break

        # Remove the guessed country from possible_targets
        possible_targets.discard(next_guess)

        # Update last_guess
        last_guess = next_guess

    # When 'done' is entered or loop is broken, do nothing else

##############################################################################
# 5) RUN MAIN
##############################################################################

if __name__ == "__main__":
    main()
