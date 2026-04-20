import copernicusmarine
import pandas as pd
import numpy as np
import json
import sys
import math
import ssl
import urllib.request
import urllib.parse

try:
    import certifi
    _SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    _SSL_CTX = ssl.create_default_context()
from datetime import datetime, timedelta, timezone
from pathlib import Path

# Default coordinates of the reference point (A Coruña coast, Galicia).
# Overridden when a `site=<id>` CLI arg matches an entry in dive_sites.json.
LAT, LON = 43.382167, -8.389000

# Compass bearing (degrees, 0=N, clockwise) for the exposure strings in dive_sites.json.
COMPASS_DEG = {
    "N":   0, "NNE":  22.5, "NE":  45, "ENE":  67.5,
    "E":  90, "ESE": 112.5, "SE": 135, "SSE": 157.5,
    "S": 180, "SSW": 202.5, "SW": 225, "WSW": 247.5,
    "W": 270, "WNW": 292.5, "NW": 315, "NNW": 337.5,
}

# Residual-energy multipliers for each site's shelter class. Head-on swell still
# reaches the site at shelter_max; fully oblique swell is attenuated to shelter_min
# (diffraction around headlands means some energy always arrives).
SHELTER_CLASS = {
    "exposed":        {"min": 0.15, "max": 1.00},
    "semi-sheltered": {"min": 0.05, "max": 0.60},
    "sheltered":      {"min": 0.02, "max": 0.30},
}

# Copernicus Marine IBI datasets
DS_WAV   = "cmems_mod_ibi_wav_anfc_0.027deg_PT1H-i"
DS_PHY   = "cmems_mod_ibi_phy_anfc_0.027deg-3D_PT1H-m"
DS_BGC   = "cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1D-m"
DS_OPT   = "cmems_mod_ibi_bgc-optics_anfc_0.027deg_P1D-m"
DS_TIDE  = "cmems_mod_ibi_phy_anfc_0.027deg-2D_PT15M-i"
# Ocean-colour observations — preferred over model optics when available.
# L4 gap-free multisensor 4 km daily product, global coverage.
DS_SAT_KD = "cmems_obs-oc_glo_bgc-transp_nrt_l4-gapfree-multi-4km_P1D"

# Visibility / physics constants
SECCHI_COEFF      = 1.7      # Poole-Atkins coefficient
SEDIMENT_HALFLIFE_H = 24     # hours for exponential wave-history weighting
DEFAULT_HOUR      = 9
DEFAULT_DEPTH     = 12.0     # metres — default dive depth for orbital velocity calc
GRAVITY           = 9.81     # m/s²

# Sediment-resuspension model constants
UB_ALPHA          = 0.3      # kd per (m/s) above threshold
UB_THRESHOLD      = 0.15     # m/s — orbital velocity threshold for resuspension

# Runoff / atmospheric constants
RAIN_LOOKBACK_H   = 72       # hours of precipitation accumulated ahead of target
RAIN_KD_COEFF     = 0.008    # kd per mm·exp(-age/tau) of recent rain
RAIN_TAU_H        = 36       # exponential decay half-life of runoff effect

# Wind-driven surface mixing: extra kd in the top layer, attenuated with depth
WIND_KD_COEFF     = 0.02     # kd per (m/s) over threshold
WIND_THRESHOLD_MS = 3.0      # m/s below which wind effect is negligible
WIND_DEPTH_SCALE  = 5.0      # metres — e-folding depth of surface mixing


def parse_args():
    """Parse arguments: [YYYY-MM-DD | tomorrow] [HH] [site=<id>]
    Examples:
        python3 oceanographic_engine.py                          -> today 9:00, default point
        python3 oceanographic_engine.py 2026-05-01 14 site=el_grelle
    """
    today = datetime.now(timezone.utc).date()
    target_date = today
    hour = DEFAULT_HOUR
    site_id = None

    for arg in sys.argv[1:]:
        if arg == "tomorrow":
            target_date = today + timedelta(days=1)
            continue
        if arg.startswith("site="):
            site_id = arg.split("=", 1)[1].strip() or None
            continue
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

    return target_date, hour, site_id


def load_site(site_id):
    """Load a site definition from dive_sites.json. Returns None if not found."""
    if not site_id:
        return None
    sites_file = Path(__file__).parent / "dive_sites.json"
    try:
        with open(sites_file) as f:
            sites = json.load(f).get("sites", [])
    except FileNotFoundError:
        return None
    for s in sites:
        if s.get("id") == site_id:
            return s
    return None


def _angular_diff(a, b):
    """Smallest angular difference between two bearings in degrees [0, 180]."""
    d = abs((a - b) % 360)
    return d if d <= 180 else 360 - d


def site_shelter_factor(wave_dir, site):
    """Fraction of open-water wave height that actually reaches the site.

    Combines two effects:
      - Geometric: cos²(diff/2), which stays near 1 for small off-angles and
        drops sharply past 90°, matching how a headland shelters a site.
      - Shelter class: residual energy bounds from SHELTER_CLASS, which encode
        that even perfectly-sheltered sites still see some swell via diffraction.
    """
    if not site or wave_dir is None:
        return 1.0
    exp = site.get("exposure")
    if exp not in COMPASS_DEG:
        return 1.0
    diff = _angular_diff(wave_dir, COMPASS_DEG[exp])   # 0..180
    alignment = math.cos(math.radians(diff / 2.0)) ** 2  # 1 head-on → 0 opposite
    bounds = SHELTER_CLASS.get(site.get("shelter", "exposed"),
                               SHELTER_CLASS["exposed"])
    return bounds["min"] + (bounds["max"] - bounds["min"]) * alignment


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
    except Exception as e:
        print(f"[fetch_nearest] {ds_id} vars={vars}: {type(e).__name__}: {e}",
              file=sys.stderr)
        return pd.DataFrame()


def fetch_satellite_kd(target, lat=None, lon=None, R=0.1):
    """Return (kd_effective, source_date, zsd_m) observed by satellite for target.

    Uses the Copernicus L4 gap-free multisensor transparency product (~4 km).
    Prefers the ZSD (Secchi disk depth) variable, which for Case-2 coastal
    waters is more faithful to real in-water optics than the spectral KD490.
    An effective kd is derived via kd = SECCHI_COEFF / ZSD so the rest of the
    pipeline stays consistent with the Poole-Atkins formula.

    The product lags by ~2 days behind real-time, so for future or too-recent
    targets this returns (None, None, None). Matches only the exact target
    date to avoid silently injecting stale data.
    """
    ref_lat = LAT if lat is None else lat
    ref_lon = LON if lon is None else lon
    try:
        df = copernicusmarine.read_dataframe(
            dataset_id=DS_SAT_KD, variables=["KD490", "ZSD"],
            minimum_longitude=ref_lon - R, maximum_longitude=ref_lon + R,
            minimum_latitude=ref_lat - R,  maximum_latitude=ref_lat + R,
            start_datetime=(target - timedelta(days=1)).strftime("%Y-%m-%dT%H:%M:%S"),
            end_datetime  =(target + timedelta(days=1)).strftime("%Y-%m-%dT%H:%M:%S"),
        ).reset_index()
    except Exception as e:
        print(f"[fetch_satellite_kd] {DS_SAT_KD}: {type(e).__name__}: {e}",
              file=sys.stderr)
        return None, None, None
    if df.empty:
        return None, None, None
    target_day = target.astimezone(timezone.utc).date()
    df["obs_day"] = pd.to_datetime(df["time"]).dt.date
    same_day = df[df["obs_day"] == target_day].dropna(subset=["ZSD"])
    if same_day.empty:
        return None, None, None
    same_day = same_day.assign(
        dist=(same_day["latitude"] - ref_lat) ** 2 + (same_day["longitude"] - ref_lon) ** 2
    )
    row = same_day.sort_values("dist").iloc[0]
    zsd = float(row["ZSD"])
    if zsd <= 0:
        return None, None, None
    kd_eff = SECCHI_COEFF / zsd
    return kd_eff, target_day.isoformat(), zsd


def fetch_open_meteo(target, lat=None, lon=None):
    """Return dict with 72h accumulated precipitation (mm, exp-decay weighted)
    and hourly wind speed / direction at target hour.

    Uses Open-Meteo free API (no auth). Returns None on failure.
    """
    ref_lat = LAT if lat is None else lat
    ref_lon = LON if lon is None else lon
    target_utc = target.astimezone(timezone.utc).replace(minute=0, second=0, microsecond=0)
    start = (target_utc - timedelta(hours=RAIN_LOOKBACK_H))
    end   = target_utc
    params = {
        "latitude":  f"{ref_lat:.4f}",
        "longitude": f"{ref_lon:.4f}",
        "hourly":    "precipitation,wind_speed_10m,wind_direction_10m",
        "timezone":  "UTC",
        "start_date": start.strftime("%Y-%m-%d"),
        "end_date":   end.strftime("%Y-%m-%d"),
        "wind_speed_unit": "ms",
    }
    url = "https://api.open-meteo.com/v1/forecast?" + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=15, context=_SSL_CTX) as r:
            data = json.loads(r.read())
    except Exception as e:
        print(f"[fetch_open_meteo] {type(e).__name__}: {e}", file=sys.stderr)
        return None

    h = data.get("hourly", {})
    times = h.get("time", [])
    if not times:
        return None
    precip = h.get("precipitation") or []
    wspd   = h.get("wind_speed_10m") or []
    wdir   = h.get("wind_direction_10m") or []

    # Exponentially-decayed rain accumulation relative to target hour
    accum = 0.0
    for t_str, p in zip(times, precip):
        if p is None:
            continue
        t = datetime.fromisoformat(t_str).replace(tzinfo=timezone.utc)
        age_h = (target_utc - t).total_seconds() / 3600.0
        if age_h < 0 or age_h > RAIN_LOOKBACK_H:
            continue
        accum += p * math.exp(-age_h / RAIN_TAU_H)

    # Wind at target hour (fall back to any available value in ±2h window)
    target_str = target_utc.strftime("%Y-%m-%dT%H:00")
    wind_speed = wind_dir = None
    try:
        idx = times.index(target_str)
        if idx < len(wspd) and wspd[idx] is not None:
            wind_speed = float(wspd[idx])
        if idx < len(wdir) and wdir[idx] is not None:
            wind_dir = float(wdir[idx])
    except ValueError:
        pass

    return {
        "rain_accum_mm":   float(accum),
        "wind_speed_ms":   wind_speed,
        "wind_dir_deg":    wind_dir,
        "rain_window_h":   RAIN_LOOKBACK_H,
    }


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

class MissingBGCDataError(RuntimeError):
    """Raised when the Copernicus BGC/optics data needed for visibility is absent."""


def calc_visibility(kd_copernicus, H_ww, T_ww, salinity,
                    depth=DEFAULT_DEPTH,
                    H_sw=0.0, T_sw=0.0, H_eff=None,
                    rain_accum_mm=0.0, wind_speed_ms=None):
    """Multi-factor underwater visibility model for Case 2 coastal waters.

    Requires a Copernicus kd baseline (satellite or model). Adds:
      - Sediment resuspension from bottom orbital velocity of both wind-wave
        and swell components, scaled by H_eff/H_inst to carry yesterday's
        storm into today's kd.
      - Freshwater-runoff penalty driven by salinity anomaly.
      - Recent-rain runoff penalty (Open-Meteo precipitation accumulation).
      - Wind-driven surface mixing, attenuated with depth.

    Raises `MissingBGCDataError` if kd is unavailable; no fallback is applied
    because that would hide missing data as a plausible number.
    """
    if kd_copernicus is None or not (kd_copernicus > 0):
        raise MissingBGCDataError(
            "kd baseline (satellite or Copernicus optics) unavailable — cannot compute visibility"
        )
    kd_base = float(kd_copernicus)

    # --- Sediment resuspension: use the stronger of wind-wave and swell ---
    U_b_ww = calc_bottom_orbital_velocity(H_ww, T_ww, depth)
    U_b_sw = calc_bottom_orbital_velocity(H_sw, T_sw, depth)
    U_b    = max(U_b_ww, U_b_sw)

    # --- Carry-over of previously resuspended sediment (24h half-life) ---
    H_inst = max(H_ww, H_sw)
    if H_eff is not None and H_inst > 0 and H_eff > H_inst:
        U_b *= H_eff / H_inst

    kd_sediment = UB_ALPHA * max(0.0, U_b - UB_THRESHOLD)

    # --- Salinity-driven runoff penalty ---
    if salinity is not None and salinity < 34.5:
        kd_runoff_sal = min(1.0, 0.15 * (35.5 - salinity))
    else:
        kd_runoff_sal = 0.0

    # --- Recent-rain runoff penalty (exp-decay weighted mm over 72h) ---
    kd_runoff_rain = RAIN_KD_COEFF * max(0.0, float(rain_accum_mm or 0.0))

    # --- Wind-driven surface mixing, attenuated with depth ---
    if wind_speed_ms is not None and wind_speed_ms > WIND_THRESHOLD_MS:
        kd_wind = (WIND_KD_COEFF * (wind_speed_ms - WIND_THRESHOLD_MS)
                   * math.exp(-max(depth, 0.0) / WIND_DEPTH_SCALE))
    else:
        kd_wind = 0.0

    kd_runoff = kd_runoff_sal + kd_runoff_rain
    kd_total  = kd_base + kd_sediment + kd_runoff + kd_wind
    vis = SECCHI_COEFF / max(kd_total, 0.01)
    return (round(max(0.5, min(vis, 15.0)), 1),
            kd_base, kd_sediment, kd_runoff, kd_wind)


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


def extract_point(wave_df, phy_df, bgc_df, target, depth=DEFAULT_DEPTH, opt_df=None,
                  site=None, sat_kd=None, meteo=None):
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
        zeu_raw = b.get('zeu')
        if zeu_raw is not None and not (isinstance(zeu_raw, float) and math.isnan(zeu_raw)):
            zeu = float(zeu_raw)

    if opt_df is not None and not opt_df.empty:
        o = opt_df.copy()
        o["time_dist"] = abs(
            (o["time"] - target.replace(tzinfo=None)).dt.total_seconds()
        )
        if 'depth' in o.columns:
            o['depth_dist'] = (o['depth'] - depth).abs()
            sort_cols = ['time_dist', 'depth_dist', 'dist']
        else:
            sort_cols = ['time_dist', 'dist']
        kd_raw = o.sort_values(sort_cols).iloc[0].get('kd')
        if kd_raw is not None and not (isinstance(kd_raw, float) and math.isnan(kd_raw)):
            kd_cop = float(kd_raw)

    # --- Site sheltering: attenuate swell reaching this specific site ---
    shelter = site_shelter_factor(d, site)  # 1.0 when no site is specified
    H_ww_site  = H_ww  * shelter
    H_sw_site  = swell_h * shelter
    H_eff_site = effective_wave * shelter

    # --- kd baseline: prefer satellite observation over model when available ---
    kd_baseline = sat_kd if sat_kd is not None else kd_cop

    # --- Atmospheric inputs ---
    rain_mm    = (meteo or {}).get("rain_accum_mm") or 0.0
    wind_speed = (meteo or {}).get("wind_speed_ms")

    # --- Visibility ---
    vis, kd_base, kd_sed, kd_run, kd_wnd = calc_visibility(
        kd_baseline, H_ww_site, T_ww, salinity, depth=depth,
        H_sw=H_sw_site, T_sw=swell_per, H_eff=H_eff_site,
        rain_accum_mm=rain_mm, wind_speed_ms=wind_speed,
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
        "kd_base": kd_base, "kd_sed": kd_sed, "kd_run": kd_run, "kd_wnd": kd_wnd,
        "kd_source": "satellite" if sat_kd is not None else ("model" if kd_cop is not None else None),
        # scores
        "dive_score": dive_score,
        "vis_factor": vis_factor, "temp_factor": temp_factor,
        "curr_factor": curr_factor, "wave_factor": wave_factor,
        "sal_factor": sal_factor,
        "surf_score": surf_score,
        # site sheltering
        "site_shelter_factor": shelter,
        # atmospheric inputs (may be None when unavailable)
        "rain_accum_mm": rain_mm,
        "wind_speed_ms": wind_speed,
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
    global LAT, LON

    target_date, hour, site_id = parse_args()
    site = load_site(site_id)
    if site_id and site is None:
        print(f"[oceanographic_engine] unknown site id: {site_id!r}", file=sys.stderr)
        sys.exit(2)
    if site:
        LAT, LON = site["lat"], site["lon"]
        site_depth = float(site.get("depth_m", DEFAULT_DEPTH))
    else:
        site_depth = DEFAULT_DEPTH

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
    # Fetch BGC data: daily, small window around target.
    # `kd` lives in the separate optics dataset; chl and zeu are in the BGC one.
    # -----------------------------------------------------------------------
    bgc_df = fetch_nearest(
        DS_BGC,
        ["chl", "zeu"],
        target - timedelta(hours=24),
        target + timedelta(hours=24)
    )
    opt_df = fetch_nearest(
        DS_OPT,
        ["kd"],
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
    # Observational / atmospheric augmentations (best-effort, never raise)
    # -----------------------------------------------------------------------
    sat_kd, sat_kd_date, sat_zsd = fetch_satellite_kd(target)
    meteo = fetch_open_meteo(target)

    # -----------------------------------------------------------------------
    # Extract point data at target, target+3h, target+6h
    # -----------------------------------------------------------------------
    try:
        pt_now = extract_point(wave_df, phy_df, bgc_df, target, opt_df=opt_df,
                               depth=site_depth, site=site,
                               sat_kd=sat_kd, meteo=meteo)
        pt_3h  = extract_point(wave_df, phy_df, bgc_df, target + timedelta(hours=3),
                               opt_df=opt_df, depth=site_depth, site=site,
                               sat_kd=sat_kd, meteo=meteo)
        pt_6h  = extract_point(wave_df, phy_df, bgc_df, target + timedelta(hours=6),
                               opt_df=opt_df, depth=site_depth, site=site,
                               sat_kd=sat_kd, meteo=meteo)
    except MissingBGCDataError as e:
        err = {
            "error": "missing_bgc_data",
            "message": str(e),
            "target":  target.strftime("%Y-%m-%d %H:%M UTC"),
            "location": {"lat": LAT, "lon": LON},
            "hint": (
                "Copernicus IBI BGC/optics datasets are updated weekly (Thursdays). "
                "Requested date is likely beyond the current forecast horizon."
            ),
        }
        print(json.dumps(err, indent=2))
        print(f"[oceanographic_engine] {err['message']}", file=sys.stderr)
        sys.exit(2)

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
        "kd_satellite":         round(sat_kd, 4) if sat_kd is not None else None,
        "secchi_satellite_m":   round(sat_zsd, 1) if sat_zsd is not None else None,
        "kd_source":            pt_now['kd_source'],
        "sat_observation_date": sat_kd_date,
        "euphotic_depth_m":     round(pt_now['zeu'], 1) if pt_now['zeu'] is not None else None,
        "rain_prev_72h_mm":     round(pt_now['rain_accum_mm'], 2),
        "wind_speed_ms":        round(pt_now['wind_speed_ms'], 2) if pt_now['wind_speed_ms'] is not None else None,
        "kd_components": {
            "baseline": round(pt_now['kd_base'], 4),
            "sediment": round(pt_now['kd_sed'],  4),
            "runoff":   round(pt_now['kd_run'],  4),
            "wind":     round(pt_now['kd_wnd'],  4),
        },
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
    if site:
        res["site"] = {
            "id":        site["id"],
            "name":      site["name"],
            "exposure":  site.get("exposure"),
            "shelter":   site.get("shelter"),
            "depth_m":   site_depth,
            "shelter_factor": round(pt_now["site_shelter_factor"], 3),
        }

    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    get_data()
