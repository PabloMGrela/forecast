import copernicusmarine
import pandas as pd
import numpy as np
import json
import sys
import math
from datetime import datetime, timedelta, timezone

# Coordinates of the reference point (A Coruña coast, Galicia)
LAT, LON = 43.382167, -8.389000

# Copernicus Marine IBI datasets
DS_WAV   = "cmems_mod_ibi_wav_anfc_0.027deg_PT1H-i"
DS_PHY   = "cmems_mod_ibi_phy_anfc_0.027deg-3D_PT1H-m"
DS_BGC   = "cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1D-m"
DS_TIDE  = "cmems_mod_ibi_phy_anfc_0.027deg-2D_PT15M-i"

# Visibility / physics constants
SECCHI_COEFF      = 1.7      # Poole-Atkins coefficient
SEDIMENT_HALFLIFE_H = 24     # hours for exponential wave-history weighting
DEFAULT_HOUR      = 9
DEFAULT_DEPTH     = 12.0     # metres — default dive depth for orbital velocity calc
GRAVITY           = 9.81     # m/s²

# Sediment-resuspension model constants
UB_ALPHA          = 0.3      # kd per (m/s) above threshold
UB_THRESHOLD      = 0.15     # m/s — orbital velocity threshold for resuspension


def parse_args():
    """Parse arguments: [YYYY-MM-DD | tomorrow] [HH]
    Examples:
        python3 oceanographic_engine.py                  -> today 9:00
        python3 oceanographic_engine.py tomorrow         -> tomorrow 9:00
        python3 oceanographic_engine.py 2026-05-01       -> 2026-05-01 9:00
        python3 oceanographic_engine.py 2026-05-01 14    -> 2026-05-01 14:00
        python3 oceanographic_engine.py 14               -> today 14:00
    """
    today = datetime.now(timezone.utc).date()
    target_date = today
    hour = DEFAULT_HOUR

    for arg in sys.argv[1:]:
        if arg == "tomorrow":
            target_date = today + timedelta(days=1)
        else:
            try:
                target_date = datetime.strptime(arg, "%Y-%m-%d").date()
                continue
            except ValueError:
                pass
            try:
                h = int(arg)
                if 0 <= h <= 23:
                    hour = h
            except ValueError:
                pass

    return target_date, hour


def fetch_nearest(ds_id, vars, time_start, time_end, R=0.2):
    """Fetch the nearest ocean point to the reference coordinates."""
    try:
        df = copernicusmarine.read_dataframe(
            dataset_id=ds_id, variables=vars,
            minimum_longitude=LON - R, maximum_longitude=LON + R,
            minimum_latitude=LAT - R,  maximum_latitude=LAT + R,
            start_datetime=time_start.strftime("%Y-%m-%dT%H:%M:%S"),
            end_datetime=time_end.strftime("%Y-%m-%dT%H:%M:%S")
        ).reset_index().dropna(subset=[vars[0]])
        if df.empty:
            return pd.DataFrame()
        df['dist'] = (df['latitude'] - LAT)**2 + (df['longitude'] - LON)**2
        return df
    except:
        return pd.DataFrame()


# ---------------------------------------------------------------------------
# Wave history helper
# ---------------------------------------------------------------------------

def calc_effective_wave_height(wave_df, target):
    """Exponentially weighted average of VHM0 over the past 5 days.

    Sediments resuspended by large waves don't disappear instantly; this
    models their settling with a 24-hour half-life.
    """
    if wave_df.empty:
        return 0.0

    nearest = wave_df.loc[wave_df.groupby('time')['dist'].idxmin()].copy()
    nearest = nearest.sort_values('time')

    nearest["hours_ago"] = (
        target.replace(tzinfo=None) - nearest["time"]
    ).dt.total_seconds() / 3600
    nearest = nearest[nearest['hours_ago'] >= 0]

    if nearest.empty:
        return 0.0

    weights = np.exp(-nearest['hours_ago'].values / SEDIMENT_HALFLIFE_H)
    effective = np.average(nearest['VHM0'].values, weights=weights)
    return float(effective)


# ---------------------------------------------------------------------------
# Orbital-velocity bottom sediment resuspension
# ---------------------------------------------------------------------------

def _wavelength_linear(T, depth, iterations=10):
    """Iterative linear dispersion relation to find wavelength at given depth.

    Deep-water approximation as seed, then iterate:
        L = (g * T^2 / 2π) * tanh(2π * depth / L)
    """
    L0 = GRAVITY * T**2 / (2 * math.pi)  # deep-water seed
    L = L0
    for _ in range(iterations):
        L = L0 * math.tanh(2 * math.pi * depth / L)
    return max(L, 0.1)


def calc_bottom_orbital_velocity(H_ww, T_ww, depth=DEFAULT_DEPTH):
    """Near-bottom peak orbital velocity from linear wave theory.

    U_b = π * H / (T * sinh(2π * d / L))
    """
    if T_ww <= 0 or H_ww <= 0:
        return 0.0
    L = _wavelength_linear(T_ww, depth)
    arg = 2 * math.pi * depth / L
    sinh_arg = math.sinh(arg)
    if sinh_arg < 1e-6:
        return 0.0
    U_b = math.pi * H_ww / (T_ww * sinh_arg)
    return float(U_b)


# ---------------------------------------------------------------------------
# Improved visibility model
# ---------------------------------------------------------------------------

def calc_visibility(kd_copernicus, H_ww, T_ww, salinity,
                    depth=DEFAULT_DEPTH, chl=None):
    """Multi-factor underwater visibility model for Case 2 coastal waters.

    Priority:
      1. Use Copernicus kd (BGC) as the baseline when available.
      2. Add sediment resuspension from wind-wave bottom orbital velocity.
      3. Add freshwater-runoff penalty from salinity.

    Falls back to a bio-optical estimate from chlorophyll if kd is missing.
    """
    # --- Baseline kd ---
    if kd_copernicus is not None and kd_copernicus > 0:
        kd_base = float(kd_copernicus)
    else:
        # Bio-optical fallback (Lee et al. / Morel)
        c = max(chl, 0.01) if chl is not None else 0.1
        kd_base = 0.12 + 0.0166 + 0.0773 * c**0.6715

    # --- Sediment resuspension from wind waves ---
    U_b = calc_bottom_orbital_velocity(H_ww, T_ww, depth)
    kd_sediment = UB_ALPHA * max(0.0, U_b - UB_THRESHOLD)

    # --- Runoff / freshwater penalty ---
    if salinity is not None and salinity < 34.5:
        kd_runoff = min(1.0, 0.15 * (35.5 - salinity))
    else:
        kd_runoff = 0.0

    kd_total = kd_base + kd_sediment + kd_runoff
    vis = SECCHI_COEFF / max(kd_total, 0.01)
    return round(max(0.5, min(vis, 15.0)), 1), kd_base, kd_sediment, kd_runoff


# ---------------------------------------------------------------------------
# Tide analysis from 15-min zos time series
# ---------------------------------------------------------------------------

def analyse_tide(tide_df, target):
    """Derive tide state and next event from a 13-hour zos window.

    Returns a dict with height_m, state, and next_event.
    """
    if tide_df.empty:
        return None

    # Keep nearest grid point, sort by time
    nearest_idx = tide_df.groupby('time')['dist'].idxmin()
    ts = tide_df.loc[nearest_idx].copy().sort_values('time').reset_index(drop=True)

    if 'zos' not in ts.columns or ts['zos'].isna().all():
        return None

    # Convert time column to UTC-aware if needed
    if ts['time'].dt.tz is None:
        ts['time'] = ts['time'].dt.tz_localize('UTC')

    # Height at target (nearest 15-min step)
    ts['time_dist'] = (ts['time'] - target).dt.total_seconds().abs()
    idx_now = ts['time_dist'].idxmin()
    height_now = float(ts.loc[idx_now, 'zos'])

    # Rising / falling: compare ±1h neighbours
    t_minus = target - timedelta(hours=1)
    t_plus  = target + timedelta(hours=1)

    def _closest_zos(t_ref):
        d = (ts['time'] - t_ref).dt.total_seconds().abs()
        i = d.idxmin()
        if d[i] > 3600:
            return None
        return float(ts.loc[i, 'zos'])

    z_before = _closest_zos(t_minus)
    z_after  = _closest_zos(t_plus)

    if z_before is not None and z_after is not None:
        delta = z_after - z_before
        if abs(delta) < 0.05:
            state = "high" if height_now > (z_before + z_after) / 2 else "low"
        elif delta > 0:
            state = "rising"
        else:
            state = "falling"
    elif z_after is not None:
        state = "rising" if z_after > height_now else "falling"
    else:
        state = "unknown"

    # Next high/low: sign changes in first derivative of zos
    zos_vals = ts['zos'].values
    times    = ts['time'].values  # numpy datetime64

    next_event = None
    if len(zos_vals) >= 3:
        dz = np.diff(zos_vals.astype(float))
        sign_changes = np.where(np.diff(np.sign(dz)))[0]  # index into dz[:-1]
        for sc in sign_changes:
            # sc+1 is the index in ts where the extremum is located
            event_idx = sc + 1
            event_time_np = times[event_idx]
            # Convert numpy datetime64 to Python datetime
            event_time = pd.Timestamp(event_time_np).to_pydatetime()
            if event_time.tzinfo is None:
                event_time = event_time.replace(tzinfo=timezone.utc)
            if event_time <= target:
                continue
            event_type = "high" if dz[sc] > 0 else "low"
            next_event = {
                "type": event_type,
                "at": event_time.strftime("%H:%M UTC")
            }
            break

    return {
        "height_m": round(height_now, 2),
        "state": state,
        "next_event": next_event
    }


# ---------------------------------------------------------------------------
# Scoring helpers
# ---------------------------------------------------------------------------

def calc_dive_score(vis, temp, curr_speed, effective_wave, salinity):
    """Weighted dive score (0-100) from five environmental factors."""
    vis_score  = min(100.0, vis * 10.0)
    temp_score = min(100.0, max(0.0, (temp - 10.0) * 12.5))
    curr_score = max(0.0, 100.0 - curr_speed * 200.0)
    wave_score = max(0.0, 100.0 - effective_wave * 40.0)
    if salinity is not None and salinity >= 34.5:
        sal_score = 100.0
    elif salinity is not None:
        sal_score = max(0.0, (salinity - 30.0) * 22.0)
    else:
        sal_score = 100.0  # assume normal if missing

    score = (
        0.35 * vis_score
        + 0.10 * temp_score
        + 0.25 * curr_score
        + 0.20 * wave_score
        + 0.10 * sal_score
    )
    return (
        int(round(min(100.0, max(0.0, score)))),
        int(round(vis_score)),
        int(round(temp_score)),
        int(round(curr_score)),
        int(round(wave_score)),
        int(round(sal_score)),
    )


def calc_surf_score(h, pr, swell_h, swell_period):
    """Simple surf score using total and swell wave parameters."""
    score = 0
    # Wave height band
    if 0.8 <= h <= 1.8:
        score += 40
    elif 1.8 < h <= 2.5:
        score += 20
    elif h > 2.5:
        score += 5
    # Period / power
    score += min(pr * 2.5, 30)
    # Swell quality bonus
    if swell_period is not None and swell_period > 0:
        score += min(swell_period * 1.5, 20)
    else:
        score += 10  # neutral if missing
    # Swell height contribution (clean groundswell preferred)
    if swell_h is not None and 0.5 <= swell_h <= 2.0:
        score += 10
    return int(min(100, max(0, score)))


def rating(score):
    if score >= 80: return "excellent"
    if score >= 60: return "good"
    if score >= 40: return "fair"
    return "poor"


# ---------------------------------------------------------------------------
# Safety alerts
# ---------------------------------------------------------------------------

def build_safety_alerts(effective_wave, curr_speed, temp, vis, salinity, H_ww):
    """Return a list of human-readable warning strings."""
    alerts = []
    if effective_wave > 2.5:
        alerts.append("Heavy swell — shore entries dangerous")
    if curr_speed > 0.4:
        alerts.append("Strong currents")
    if temp < 13:
        alerts.append("Cold water — full wetsuit + hood required")
    if vis < 2:
        alerts.append("Very poor visibility")
    if salinity is not None and salinity < 33:
        alerts.append("Heavy river discharge — avoid coastal sites")
    if H_ww is not None and H_ww > 2.0:
        alerts.append("Rough surface conditions")
    return alerts


# ---------------------------------------------------------------------------
# Single time-point data extraction
# ---------------------------------------------------------------------------

def _scalar(row, col, default=0.0):
    """Safely extract a scalar from a DataFrame row or fallback default."""
    if row is None:
        return default
    v = row.get(col, default)
    return float(v) if v is not None and not (isinstance(v, float) and math.isnan(v)) else default


def extract_point(wave_df, phy_df, bgc_df, target, depth=DEFAULT_DEPTH):
    """Extract all variables at a single target datetime and compute scores.

    Returns a dict with all telemetry, scores and metadata.
    """
    # --- Wave variables ---
    effective_wave = calc_effective_wave_height(wave_df, target)

    if not wave_df.empty:
        wave_df = wave_df.copy()
        wave_df["time_dist"] = abs(
            (wave_df["time"] - target.replace(tzinfo=None)).dt.total_seconds()
        )
        cw = wave_df.sort_values(['time_dist', 'dist']).iloc[0]
    else:
        cw = None

    h         = _scalar(cw, 'VHM0')
    pr        = _scalar(cw, 'VTPK')
    d         = _scalar(cw, 'VMDR')
    H_ww      = _scalar(cw, 'VHM0_WW')
    T_ww      = _scalar(cw, 'VTM01_WW')
    swell_h   = _scalar(cw, 'VHM0_SW1')
    swell_per = _scalar(cw, 'VTM01_SW1')

    # --- Physics variables ---
    if not phy_df.empty:
        phy_df = phy_df.copy()
        phy_df["time_dist"] = abs(
            (phy_df["time"] - target.replace(tzinfo=None)).dt.total_seconds()
        )
        p = phy_df.sort_values(['time_dist', 'dist']).iloc[0]
    else:
        p = None

    temp     = _scalar(p, 'thetao')
    uo       = _scalar(p, 'uo')
    vo       = _scalar(p, 'vo')
    salinity = _scalar(p, 'so', None)  # None means "not available"
    if salinity == 0.0 and p is not None and p.get('so') is None:
        salinity = None

    curr_speed = float(np.sqrt(uo**2 + vo**2))
    curr_dir   = float(math.degrees(math.atan2(uo, vo)) % 360) if (uo or vo) else 0.0

    # --- BGC variables ---
    kd_cop  = None
    zeu     = None
    chl     = 0.1

    if not bgc_df.empty:
        bgc_df = bgc_df.copy()
        bgc_df["time_dist"] = abs(
            (bgc_df["time"] - target.replace(tzinfo=None)).dt.total_seconds()
        )
        b = bgc_df.sort_values(['time_dist', 'dist']).iloc[0]
        chl_raw = b.get('chl')
        if chl_raw is not None and not (isinstance(chl_raw, float) and math.isnan(chl_raw)):
            chl = float(chl_raw)
        kd_raw = b.get('kd')
        if kd_raw is not None and not (isinstance(kd_raw, float) and math.isnan(kd_raw)):
            kd_cop = float(kd_raw)
        zeu_raw = b.get('zeu')
        if zeu_raw is not None and not (isinstance(zeu_raw, float) and math.isnan(zeu_raw)):
            zeu = float(zeu_raw)

    # --- Visibility ---
    vis, kd_base, kd_sed, kd_run = calc_visibility(
        kd_cop, H_ww, T_ww, salinity, depth=depth, chl=chl
    )

    # --- Scores ---
    (dive_score,
     vis_factor, temp_factor, curr_factor,
     wave_factor, sal_factor) = calc_dive_score(
        vis, temp, curr_speed, effective_wave, salinity
    )

    surf_score = calc_surf_score(h, pr, swell_h, swell_per)

    return {
        # raw telemetry
        "h": h, "pr": pr, "d": d,
        "H_ww": H_ww, "T_ww": T_ww,
        "swell_h": swell_h, "swell_per": swell_per,
        "temp": temp, "uo": uo, "vo": vo,
        "curr_speed": curr_speed, "curr_dir": curr_dir,
        "salinity": salinity,
        "chl": chl,
        "kd_cop": kd_cop, "zeu": zeu,
        "effective_wave": effective_wave,
        # derived
        "vis": vis,
        "kd_base": kd_base, "kd_sed": kd_sed, "kd_run": kd_run,
        # scores
        "dive_score": dive_score,
        "vis_factor": vis_factor, "temp_factor": temp_factor,
        "curr_factor": curr_factor, "wave_factor": wave_factor,
        "sal_factor": sal_factor,
        "surf_score": surf_score,
    }


# ---------------------------------------------------------------------------
# Trend calculation
# ---------------------------------------------------------------------------

def calc_trend(score_now, score_3h, score_6h):
    """Compare dive scores to produce improving / worsening / stable."""
    avg_future = (score_3h + score_6h) / 2.0
    delta = avg_future - score_now
    if delta > 5:
        return "improving"
    if delta < -5:
        return "worsening"
    return "stable"


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def get_data():
    target_date, hour = parse_args()

    target = datetime(
        target_date.year, target_date.month, target_date.day,
        hour, 0, 0, tzinfo=timezone.utc
    )

    # -----------------------------------------------------------------------
    # Fetch wave data: 5 days of history + 6h future for trend
    # -----------------------------------------------------------------------
    wave_df = fetch_nearest(
        DS_WAV,
        ["VHM0", "VTPK", "VMDR", "VHM0_WW", "VTM01_WW", "VHM0_SW1", "VTM01_SW1"],
        target - timedelta(days=5),
        target + timedelta(hours=7)
    )

    # -----------------------------------------------------------------------
    # Fetch physics data: ±6h window (covers target, +3h, +6h)
    # -----------------------------------------------------------------------
    phy_df = fetch_nearest(
        DS_PHY,
        ["uo", "vo", "thetao", "so"],
        target - timedelta(hours=3),
        target + timedelta(hours=7)
    )

    # -----------------------------------------------------------------------
    # Fetch BGC data: daily, small window around target
    # -----------------------------------------------------------------------
    bgc_df = fetch_nearest(
        DS_BGC,
        ["chl", "kd", "zeu"],
        target - timedelta(hours=24),
        target + timedelta(hours=24)
    )

    # -----------------------------------------------------------------------
    # Fetch tide data: 13-hour window centred on target (15-min resolution)
    # -----------------------------------------------------------------------
    tide_df = fetch_nearest(
        DS_TIDE,
        ["zos"],
        target - timedelta(hours=6, minutes=30),
        target + timedelta(hours=6, minutes=30)
    )

    # -----------------------------------------------------------------------
    # Extract point data at target, target+3h, target+6h
    # -----------------------------------------------------------------------
    pt_now = extract_point(wave_df, phy_df, bgc_df, target)
    pt_3h  = extract_point(wave_df, phy_df, bgc_df, target + timedelta(hours=3))
    pt_6h  = extract_point(wave_df, phy_df, bgc_df, target + timedelta(hours=6))

    trend = calc_trend(pt_now['dive_score'], pt_3h['dive_score'], pt_6h['dive_score'])

    # -----------------------------------------------------------------------
    # Tide analysis
    # -----------------------------------------------------------------------
    tide = analyse_tide(tide_df, target)

    # -----------------------------------------------------------------------
    # Safety alerts
    # -----------------------------------------------------------------------
    sal = pt_now['salinity']
    alerts = build_safety_alerts(
        pt_now['effective_wave'],
        pt_now['curr_speed'],
        pt_now['temp'],
        pt_now['vis'],
        sal,
        pt_now['H_ww']
    )

    # -----------------------------------------------------------------------
    # Assemble output
    # -----------------------------------------------------------------------
    telemetry = {
        "wave_height_m":        round(pt_now['h'], 2),
        "wave_period_s":        round(pt_now['pr'], 2),
        "wave_dir_deg":         round(pt_now['d'], 1),
        "wind_wave_height_m":   round(pt_now['H_ww'], 2),
        "wind_wave_period_s":   round(pt_now['T_ww'], 2),
        "swell_height_m":       round(pt_now['swell_h'], 2),
        "swell_period_s":       round(pt_now['swell_per'], 2),
        "water_temp_c":         round(pt_now['temp'], 1),
        "salinity_psu":         round(sal, 2) if sal is not None else None,
        "chlorophyll_mg_m3":    round(pt_now['chl'], 3),
        "current_speed_mps":    round(pt_now['curr_speed'], 3),
        "current_dir_deg":      round(pt_now['curr_dir'], 1),
        "effective_wave_m":     round(pt_now['effective_wave'], 2),
        "kd_copernicus":        round(pt_now['kd_cop'], 4) if pt_now['kd_cop'] is not None else None,
        "euphotic_depth_m":     round(pt_now['zeu'], 1) if pt_now['zeu'] is not None else None,
    }

    # Backward-compat alias kept from original output
    telemetry["chlorophyll"] = telemetry["chlorophyll_mg_m3"]

    surf = {
        "score":  pt_now['surf_score'],
        "rating": rating(pt_now['surf_score'])
    }

    dive = {
        "score":        pt_now['dive_score'],
        "visibility_m": pt_now['vis'],
        "rating":       rating(pt_now['dive_score']),
        "factors": {
            "visibility":  pt_now['vis_factor'],
            "temperature": pt_now['temp_factor'],
            "currents":    pt_now['curr_factor'],
            "waves":       pt_now['wave_factor'],
            "salinity":    pt_now['sal_factor'],
        }
    }

    res = {
        "location": {"lat": LAT, "lon": LON},
        "target":   target.strftime("%Y-%m-%d %H:%M UTC"),
        "telemetry": telemetry,
        "tide":  tide,
        "surf":  surf,
        "dive":  dive,
        "trend": trend,
        "safety_alerts": alerts,
    }

    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    get_data()
