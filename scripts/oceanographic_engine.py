import copernicusmarine
import pandas as pd
import numpy as np
import json
import sys
from datetime import datetime, timedelta, timezone

# Coordinates of the reference point (A Coruña coast, Galicia)
LAT, LON = 43.382167, -8.389000

# Copernicus Marine IBI datasets
DS_WAV = "cmems_mod_ibi_wav_anfc_0.027deg_PT1H-i"
DS_PHY = "cmems_mod_ibi_phy_anfc_0.027deg-3D_PT1H-m"
DS_BGC = "cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1D-m"

# Visibility model constants
KD_BASELINE = 0.12
SEDIMENT_HALFLIFE_H = 24
SECCHI_COEFF = 1.4
DEFAULT_HOUR = 9


def parse_args():
    """Parse arguments: [tomorrow] [HH]
    Examples:
        python3 oceanographic_engine.py             -> today 9:00
        python3 oceanographic_engine.py tomorrow    -> tomorrow 9:00
        python3 oceanographic_engine.py 14          -> today 14:00
        python3 oceanographic_engine.py tomorrow 7  -> tomorrow 7:00
    """
    is_tomorrow = False
    hour = DEFAULT_HOUR
    for arg in sys.argv[1:]:
        if arg == "tomorrow":
            is_tomorrow = True
        else:
            try:
                h = int(arg)
                if 0 <= h <= 23:
                    hour = h
            except ValueError:
                pass
    return is_tomorrow, hour


def fetch_nearest(ds_id, vars, time_start, time_end, R=0.2):
    """Fetch the nearest ocean point to the reference coordinates."""
    try:
        df = copernicusmarine.read_dataframe(
            dataset_id=ds_id, variables=vars,
            minimum_longitude=LON-R, maximum_longitude=LON+R,
            minimum_latitude=LAT-R, maximum_latitude=LAT+R,
            start_datetime=time_start.strftime("%Y-%m-%dT%H:%M:%S"),
            end_datetime=time_end.strftime("%Y-%m-%dT%H:%M:%S")
        ).reset_index().dropna(subset=[vars[0]])
        if df.empty:
            return pd.DataFrame()
        df['dist'] = (df['latitude'] - LAT)**2 + (df['longitude'] - LON)**2
        return df
    except:
        return pd.DataFrame()


def calc_effective_wave_height(wave_df, target):
    """Calculate effective wave height with exponential decay.

    Sediments resuspended by large waves don't disappear instantly.
    We use an exponentially weighted average where recent waves weigh more,
    but those from the last 48-72h still contribute (modelling settling).
    """
    if wave_df.empty:
        return 0.0

    nearest = wave_df.loc[wave_df.groupby('time')['dist'].idxmin()].copy()
    nearest = nearest.sort_values('time')

    nearest["hours_ago"] = (target.replace(tzinfo=None) - nearest["time"]).dt.total_seconds() / 3600
    nearest = nearest[nearest['hours_ago'] >= 0]

    if nearest.empty:
        return 0.0

    weights = np.exp(-nearest['hours_ago'].values / SEDIMENT_HALFLIFE_H)
    effective = np.average(nearest['VHM0'].values, weights=weights)
    return float(effective)


def calc_visibility(chl, effective_wave, curr_speed):
    """Multi-factor underwater visibility model for Case 2 coastal waters."""
    kd_chl = 0.0166 + 0.0773 * max(chl, 0.01) ** 0.6715
    kd_waves = 0.05 * effective_wave
    kd_curr = 0.1 * curr_speed
    kd_total = KD_BASELINE + kd_chl + kd_waves + kd_curr
    vis = SECCHI_COEFF / kd_total
    return round(max(0.5, min(vis, 15)), 1)


def get_data():
    is_tomorrow, hour = parse_args()
    now = datetime.now(timezone.utc)

    target_date = now.date() + timedelta(days=1 if is_tomorrow else 0)
    target = datetime(target_date.year, target_date.month, target_date.day,
                      hour, 0, 0, tzinfo=timezone.utc)

    # Waves: 5 days of history to calculate effective wave height
    wave_df = fetch_nearest(
        DS_WAV, ["VHM0", "VTPK", "VMDR"],
        target - timedelta(days=5),
        target + timedelta(hours=3)
    )

    effective_wave = calc_effective_wave_height(wave_df, target)

    if not wave_df.empty:
        wave_df["time_dist"] = abs((wave_df["time"] - target.replace(tzinfo=None)).dt.total_seconds())
        current_wave = wave_df.sort_values(['time_dist', 'dist']).iloc[0]
        h = float(current_wave.get('VHM0', 0))
        pr = float(current_wave.get('VTPK', 0))
        d = float(current_wave.get('VMDR', 0))
    else:
        h, pr, d = 0, 0, 0

    # Physics: currents and temperature
    phy_df = fetch_nearest(
        DS_PHY, ["uo", "vo", "thetao"],
        target - timedelta(hours=3),
        target + timedelta(hours=3)
    )
    if not phy_df.empty:
        phy_df["time_dist"] = abs((phy_df["time"] - target.replace(tzinfo=None)).dt.total_seconds())
        p = phy_df.sort_values(['time_dist', 'dist']).iloc[0]
        t = float(p.get('thetao', 0))
        uo = float(p.get('uo', 0))
        vo = float(p.get('vo', 0))
    else:
        t, uo, vo = 0, 0, 0

    curr_speed = np.sqrt(uo**2 + vo**2)

    # Biogeochemistry: chlorophyll
    bgc_df = fetch_nearest(
        DS_BGC, ["chl"],
        target - timedelta(hours=24),
        target + timedelta(hours=24)
    )
    if not bgc_df.empty:
        chl = float(bgc_df.sort_values('dist').iloc[0].get('chl', 0.1))
    else:
        chl = 0.1

    vis = calc_visibility(chl, effective_wave, curr_speed)

    surf_score = 0
    if 0.8 <= h <= 1.8: surf_score += 40
    elif 1.8 < h <= 2.5: surf_score += 20
    surf_score += min(pr * 2.5, 30)
    surf_score += 30

    dive_score = max(0, min(100, int(
        100 - (effective_wave * 15) - (curr_speed * 50) - (chl * 10)
    )))

    def rating(score):
        if score >= 80: return "excellent"
        if score >= 60: return "good"
        if score >= 40: return "fair"
        return "poor"

    res = {
        "location": {"lat": LAT, "lon": LON},
        "target": target.strftime("%Y-%m-%d %H:%M UTC"),
        "telemetry": {
            "wave_height_m": round(float(h), 2),
            "wave_period_s": round(float(pr), 2),
            "wave_dir_deg": round(float(d), 1),
            "water_temp_c": round(float(t), 1),
            "chlorophyll": round(float(chl), 2),
            "current_speed_mps": round(float(curr_speed), 2),
            "effective_wave_m": round(effective_wave, 2)
        },
        "surf": {
            "score": int(min(100, surf_score)),
            "rating": rating(int(min(100, surf_score)))
        },
        "dive": {
            "score": dive_score,
            "visibility_m": vis,
            "rating": rating(dive_score)
        }
    }
    print(json.dumps(res, indent=2))

if __name__ == "__main__":
    get_data()
