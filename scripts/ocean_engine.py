import sys
import json
from datetime import datetime


def analyze(data):
    h = data.get('wave_height', 0)
    p = data.get('wave_period', 0)
    wd_speed = data.get('wind_speed', 0)
    wd_dir = data.get('wind_direction', 0)
    curr = data.get('current_speed', 0)
    chl = data.get('chlorophyll', 0)

    # Offshore wind for these beaches: southerly component (135° to 225°)
    is_offshore = 135 <= wd_dir <= 225

    # --- PER-BEACH EVALUATION (Beginner/Returning level) ---
    beaches = {}

    # 1. BASTIAGUEIRO (safe zone when there's big swell)
    basti_score = 0
    if h > 1.8: basti_score = 80
    elif 1.0 <= h <= 1.8: basti_score = 60
    if is_offshore: basti_score += 20
    beaches["Bastiagueiro"] = {"score": min(100, basti_score), "note": "Ideal when it's big outside."}

    # 2. MATADERO (low tide only)
    mata_score = 0
    if 0.8 <= h <= 1.6: mata_score = 70
    if is_offshore: mata_score += 20
    beaches["Matadero"] = {"score": mata_score, "note": "LOW TIDE ONLY — doesn't break at high tide."}

    # 3. CAION (versatile)
    caion_score = 0
    if 0.7 <= h <= 1.5: caion_score = 80
    if is_offshore: caion_score += 20
    beaches["Caion"] = {"score": caion_score, "note": "Great for getting back in the water with medium swell."}

    # 4. SABON (small swell days)
    sabon_score = 0
    if h < 1.0: sabon_score = 80
    elif h <= 1.5: sabon_score = 50
    if is_offshore: sabon_score += 20
    beaches["Sabon"] = {"score": sabon_score, "note": "Perfect when there's little swell outside."}

    # --- DIVE SCORE ---
    dive_score = 100 - (h * 20) - (curr * 50) - (chl * 10)
    dive_score = max(0, min(100, int(dive_score)))

    best_beach = max(beaches, key=lambda x: beaches[x]['score'])

    return {
        "level": "Beginner/Returning",
        "measurements": {
            "swell_h": round(h, 2),
            "period_s": round(p, 2),
            "wind": f"{wd_speed}km/h ({wd_dir}°)"
        },
        "surf_analysis": beaches,
        "dive_analysis": {
            "score": dive_score,
            "summary": "good" if dive_score > 60 else "average" if dive_score > 40 else "bad"
        },
        "recommendation": f"Go to {best_beach}. {beaches[best_beach]['note']}"
    }


if __name__ == "__main__":
    try:
        input_data = json.load(sys.stdin)
        result = analyze(input_data)
        print(json.dumps(result, indent=2, ensure_ascii=False))
    except Exception as e:
        print(json.dumps({"error": str(e)}))
