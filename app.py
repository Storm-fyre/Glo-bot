from flask import Flask, render_template, request, session, redirect, url_for
import code  # Your logic in 'code.py'

app = Flask(__name__)
app.secret_key = "SOME_RANDOM_SECRET_KEY_HERE"

# Load global data
distance_map = code.load_distance_data("country_distances.json")
supported = code.get_supported_countries(distance_map)

def pick_next_guess(possible_set, guesses_so_far):
    """Pick the next guess from possible_set that hasn't been guessed before."""
    for c in possible_set:
        if c not in guesses_so_far:
            return c
    return None

@app.route("/", methods=["GET", "POST"])
def index():
    # Initialize session if not present
    if "initialized" not in session:
        session["initialized"] = True
        session["possible_targets"] = list(supported)
        session["guesses"] = []
        session["last_guess"] = None
        session["message"] = ""

    possible_targets = set(session["possible_targets"])
    guesses = session["guesses"]
    last_guess = session["last_guess"]
    message = session["message"]

    if request.method == "POST":
        form_type = request.form.get("form_type", "")

        # =========== FIRST GUESS FORM =============
        if form_type == "first_guess":
            first_guess = request.form.get("first_guess_input", "").strip().lower()
            first_distance_str = request.form.get("first_distance_input", "").strip()

            if first_guess not in supported:
                session["message"] = f"'{first_guess}' is not recognized or not supported."
                return redirect(url_for("index"))

            # Parse distance
            try:
                first_dist = float(first_distance_str)
            except ValueError:
                session["message"] = "Invalid first guess distance."
                return redirect(url_for("index"))

            # Refine possible targets
            possible_targets = code.refine_distance_with_tolerance(
                distance_map, 
                possible_targets, 
                first_guess, 
                first_dist
            )
            guesses.append(first_guess)
            last_guess = first_guess

            # Pick next guess automatically if multiple remain
            if len(possible_targets) > 1:
                new_guess = pick_next_guess(possible_targets, guesses)
                if new_guess:
                    guesses.append(new_guess)
                    last_guess = new_guess

            # Store back in session
            session["possible_targets"] = list(possible_targets)
            session["guesses"] = guesses
            session["last_guess"] = last_guess
            session["message"] = ""
            return redirect(url_for("index"))

        # =========== FEEDBACK: COOLER / WARMER / ADJACENT / DONE ============
        elif form_type == "feedback":
            action = request.form.get("action", "")
            current_guess = guesses[-1] if guesses else None

            if not current_guess:
                session["message"] = "No guess to evaluate feedback for. Please restart."
                return redirect(url_for("index"))

            if action == "done":
                session["message"] = "Game ended by user."
                return redirect(url_for("index"))

            if action == "warmer":
                warmer_distance_str = request.form.get("warmer_distance_input", "").strip()
                try:
                    warmer_dist = float(warmer_distance_str)
                except ValueError:
                    session["message"] = "Invalid distance value for 'warmer'."
                    return redirect(url_for("index"))

                if len(guesses) < 2:
                    session["message"] = "Not enough guesses to compare cooler/warmer."
                    return redirect(url_for("index"))
                old_guess = guesses[-2]
                possible_targets = code.refine_warmer(distance_map, possible_targets, old_guess, current_guess)
                possible_targets = code.refine_distance_with_tolerance(distance_map, possible_targets, current_guess, warmer_dist)

            elif action == "cooler":
                if len(guesses) < 2:
                    session["message"] = "Not enough guesses to compare cooler/warmer."
                    return redirect(url_for("index"))
                old_guess = guesses[-2]
                possible_targets = code.refine_further(distance_map, possible_targets, old_guess, current_guess)

            elif action == "adjacent":
                possible_targets = code.refine_adjacent(distance_map, possible_targets, current_guess)

            # Remove the current guess from possibilities
            if current_guess in possible_targets:
                possible_targets.discard(current_guess)

            # Then pick next guess if multiple remain
            if len(possible_targets) > 1:
                new_guess = pick_next_guess(possible_targets, guesses)
                if new_guess:
                    guesses.append(new_guess)
                    last_guess = new_guess

            # Store updated session data
            session["possible_targets"] = list(possible_targets)
            session["guesses"] = guesses
            session["last_guess"] = last_guess
            session["message"] = ""
            return redirect(url_for("index"))

        # =========== RESET/RESTART ===========
        elif form_type == "reset":
            session.clear()
            return redirect(url_for("index"))

    # For GET or after POST -> show updated page
    possible_count = len(possible_targets)
    final_solution = None
    contradiction = False

    if possible_count == 0:
        contradiction = True
    elif possible_count == 1:
        # "We have EXACTLY ONE possible target" => we say "Your mystery country is:"
        final_solution = list(possible_targets)[0]

    # The current guess is always the last guessed item
    current_guess = guesses[-1] if guesses else None

    return render_template(
        "index.html",
        message=message,
        current_guess=current_guess,
        possible_count=possible_count,
        final_solution=final_solution,
        contradiction=contradiction
    )

if __name__ == "__main__":
    # Run locally: http://127.0.0.1:5000
    app.run()
