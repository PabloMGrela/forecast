import sys
import json
from pathlib import Path


def load_sites():
    """Load site definitions from dive_sites.json relative to this script."""
    sites_file = Path(__file__).parent / "dive_sites.json"
    with open(sites_file, "r") as f:
        return json.load(f)["sites"]


# Compass direction → degrees
COMPASS_DEG = {
    "N": 0, "NE": 45, "E": 90, "SE": 135,
    "S": 180, "SW": 225, "W": 270, "NW": 315,
}


def is_offshore_wind(wind_dir, site_exposure):
    """Check if wind blows away from the site (±90° from the opposite of exposure)."""
    if site_exposure not in COMPASS_DEG:
        return False
    offshore_center = (COMPASS_DEG[site_exposure] + 180) % 360
    diff = abs(wind_dir - offshore_center)
    if diff > 180:
        diff = 360 - diff
    return diff <= 90


def _normalize_input(raw):
    """Accept both flat dicts and the nested output of oceanographic_engine.py."""
    if "telemetry" in raw:
        t = raw["telemetry"]
        tide = raw.get("tide") or {}
        return {
            "wave_height":      t.get("wave_height_m", 0),
            "wave_period":      t.get("wave_period_s", 0),
            "wind_wave_height": t.get("wind_wave_height_m", 0),
            "swell_height":     t.get("swell_height_m", 0),
            "swell_period":     t.get("swell_period_s", 0),
            "wave_direction":   t.get("wave_dir_deg", 0),
            "wind_direction":   t.get("wave_dir_deg", 0),  # approx: use wave dir
            "water_temp":       t.get("water_temp_c", 15),
            "current_speed":    t.get("current_speed_mps", 0),
            "current_dir":      t.get("current_dir_deg", 0),
            "salinity":         t.get("salinity_psu") or 35.5,
            "chlorophyll":      t.get("chlorophyll_mg_m3", 0),
            "visibility_m":     raw.get("dive", {}).get("visibility_m", 5),
            "tide_state":       tide.get("state", ""),
        }
    return raw  # already flat


def score_surf_site(site, data):
    """Score a surf site based on wave and wind conditions."""
    swell_h = data.get("swell_height", 0)
    swell_p = data.get("swell_period", 0)
    wind_dir = data.get("wind_direction", 0)
    tide_state = data.get("tide_state", "")

    score = 0
    notes = []

    exposure = site.get("exposure", "")
    if exposure in ("W", "NW"):
        if 1.5 <= swell_h <= 2.5:
            score += 85
            notes.append("perfect swell height")
        elif 1.0 <= swell_h < 1.5:
            score += 70
            notes.append("good swell")
        elif swell_h < 1.0:
            score += 40
            notes.append("small swell")
        elif swell_h > 2.5:
            score += 60
            notes.append("big swell, challenging")
    else:
        if 0.7 <= swell_h <= 1.5:
            score += 80
            notes.append("ideal swell range")
        elif 1.5 < swell_h <= 2.0:
            score += 60
            notes.append("good but big")
        elif swell_h < 0.7:
            score += 50
            notes.append("small but rideable")

    if 10 <= swell_p <= 14:
        score += 15
        notes.append("excellent period")
    elif 8 <= swell_p < 10 or 14 < swell_p <= 16:
        score += 8
        notes.append("decent period")

    if is_offshore_wind(wind_dir, exposure):
        score += 10
        notes.append("offshore wind")

    if site["id"] == "matadero" and tide_state not in ("low", "falling"):
        score -= 40
        notes.append("NOT at ideal tide (needs low/falling)")

    return max(0, min(100, score)), " + ".join(notes) or "Conditions present"


def _wave_exposure_penalty(wave_dir, site_exposure):
    """How directly the swell hits this site. Returns 0 (sheltered) to 1 (head-on).

    Compares the wave direction (where waves come FROM) with the site exposure
    (the direction the site faces / is open to). A site facing NW hit by NW
    swell gets full penalty; a NE-facing site is sheltered from NW swell.
    """
    if site_exposure not in COMPASS_DEG or wave_dir is None:
        return 0.5  # unknown → neutral
    exp_deg = COMPASS_DEG[site_exposure]
    diff = abs(wave_dir - exp_deg)
    if diff > 180:
        diff = 360 - diff
    # 0° diff = head-on (1.0), 90° = glancing (0.3), 180° = fully sheltered (0.0)
    return max(0.0, 1.0 - diff / 180.0)


def score_dive_site(site, data):
    """Score a dive site based on water conditions and site characteristics."""
    visibility = data.get("visibility_m", 5)
    water_temp = data.get("water_temp", 15)
    current_speed = data.get("current_speed", 0)
    salinity = data.get("salinity", 35.5)
    wave_dir = data.get("wave_direction", data.get("wind_direction", 0))
    wave_height = data.get("wave_height", 0)
    shelter = site.get("shelter", "exposed")
    exposure = site.get("exposure", "")

    notes = []
    vis_adj = visibility

    # --- Shelter visibility bonus ---
    if shelter == "sheltered":
        vis_adj *= 1.30
        notes.append("sheltered +30% vis")
    elif shelter == "semi-sheltered":
        vis_adj *= 1.15
        notes.append("semi-sheltered +15% vis")
    vis_adj = min(vis_adj, 15.0)

    # Visibility → base score
    if vis_adj >= 8:
        score = 90
    elif vis_adj >= 5:
        score = 75
    elif vis_adj >= 3:
        score = 55
    else:
        score = 35

    # --- WAVE DIRECTION vs SITE EXPOSURE (biggest factor) ---
    exp_factor = _wave_exposure_penalty(wave_dir, exposure)
    # Scale penalty by wave height: small waves don't matter even head-on
    swell_penalty = int(exp_factor * wave_height * 15)
    if swell_penalty > 5:
        score -= swell_penalty
        notes.append(f"swell exposure -{swell_penalty}")
    elif exp_factor < 0.3:
        bonus = int((1 - exp_factor) * 8)
        score += bonus
        notes.append(f"sheltered from swell +{bonus}")

    # --- Temperature ---
    if water_temp < 13:
        penalty = int((13 - water_temp) * 3)
        score -= penalty
        notes.append(f"cold water -{penalty}")
    elif water_temp < 15:
        score -= 5
        notes.append("cool water -5")

    # --- Currents (threshold scaled by shelter) ---
    threshold = {"sheltered": 0.5, "semi-sheltered": 0.4}.get(shelter, 0.3)
    if current_speed > threshold:
        mult = {"sheltered": 5, "semi-sheltered": 10}.get(shelter, 20)
        penalty = int(current_speed * mult)
        score -= penalty
        notes.append(f"current penalty ({shelter}) -{penalty}")

    # --- Salinity — ría sites hardest hit by runoff ---
    if salinity < 34 and site["id"] in ("seselle", "ares_ria"):
        score -= 20
        notes.append("low salinity (river discharge) -20")

    score = max(0, min(100, score))
    return score, round(vis_adj, 1), " + ".join(notes) or "Good conditions"


def rating(score):
    if score >= 80: return "excellent"
    if score >= 60: return "good"
    if score >= 40: return "fair"
    return "poor"


def analyze(raw_data):
    """Analyze both surf and dive conditions for all sites."""
    data = _normalize_input(raw_data)
    sites = load_sites()

    surf_results = []
    dive_results = []

    for site in sites:
        st = site.get("type", "")

        if st in ("surf", "both"):
            sc, note = score_surf_site(site, data)
            surf_results.append({
                "name": site["name"], "score": sc,
                "rating": rating(sc), "note": note,
            })

        if st in ("dive", "both"):
            sc, vis, note = score_dive_site(site, data)
            dive_results.append({
                "name": site["name"], "score": sc,
                "rating": rating(sc), "visibility_m": vis, "note": note,
            })

    surf_results.sort(key=lambda x: x["score"], reverse=True)
    dive_results.sort(key=lambda x: x["score"], reverse=True)

    return {
        "surf_spots": surf_results,
        "dive_spots": dive_results,
        "best_surf": surf_results[0]["name"] if surf_results else None,
        "best_dive": dive_results[0]["name"] if dive_results else None,
    }


if __name__ == "__main__":
    try:
        result = analyze(json.load(sys.stdin))
        print(json.dumps(result, indent=2, ensure_ascii=False))
    except Exception as e:
        print(json.dumps({"error": str(e)}), file=sys.stderr)
        sys.exit(1)
