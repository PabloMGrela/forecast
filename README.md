# diving-forecast

Oceanographic forecast engine for surf and dive conditions on the Galician coast (A Coruña area). Uses real-time data from [Copernicus Marine Service](https://marine.copernicus.eu/) IBI (Iberian-Biscay-Irish) regional models with a sophisticated visibility model and safety alerting.

## Architecture

| Component | Purpose |
|-----------|---------|
| `oceanographic_engine.py` | Main engine — fetches multi-model data from Copernicus, computes forecasts with visibility and safety analysis |
| `ocean_engine.py` | Lightweight per-site scorer — reads stdin JSON, outputs per-dive-site recommendations |
| `dive_sites.json` | Database of dive and surf sites on the Galician coast with coordinates and depths |
| `surf.sh` | Shell wrapper for `oceanographic_engine.py` |
| `weather.sh` | Atmospheric forecast via wttr.in (no API key needed) |

## Data sources

- **Waves**: `cmems_mod_ibi_wav_anfc_0.027deg_PT1H-i` — Hm0, peak period, direction, wind waves (VHM0_WW), primary swell (VHM0_SW1)
- **Physics**: `cmems_mod_ibi_phy_anfc_0.027deg-3D_PT1H-m` — currents (U/V), temperature, salinity
- **Tides**: `cmems_mod_ibi_phy_anfc_0.027deg-2D_PT15M-i` — sea surface height at 15-min resolution for tidal state
- **Biogeochemistry**: `cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1D-m` — chlorophyll, Kd (light attenuation), euphotic depth

## Visibility model

Multi-factor attenuation for Case 2 coastal waters:

- **Baseline**: Copernicus `kd` from PISCES biogeochemical model
- **Sediment resuspension**: Wind waves increase bottom orbital velocity (depends on height, period, site depth)
- **River runoff**: Salinity anomaly proxy for Galician rías freshwater input
- **Secchi depth**: Poole-Atkins formula: `visibility = 1.7 / kd_total`

## Dive score

Weighted composite of five factors:
- Visibility (35%)
- Currents (25%)
- Waves (20%)
- Temperature (10%)
- Salinity (10%)

Score ≥70 = good, 50–69 = fair, <50 = poor.

## Features

- **Arbitrary target dates**: `python3 scripts/oceanographic_engine.py 2026-05-01 14`
- **Tidal state**: Rising/falling/high/low with time to next event
- **Safety alerts**: Heavy swell, strong currents, cold water, poor visibility, river discharge
- **6-hour trend**: Improving/worsening/stable for dive conditions
- **Per-site analysis**: Pipe to `ocean_engine.py` for per-dive-site recommendations

## Setup

```bash
pip install -r requirements.txt
copernicusmarine login
```

## Usage

```bash
# Today at 9:00 UTC
python3 scripts/oceanographic_engine.py

# Specific date and hour
python3 scripts/oceanographic_engine.py 2026-05-01 14

# Per-site recommendations
python3 scripts/oceanographic_engine.py 2026-05-01 14 | python3 scripts/ocean_engine.py

# Atmospheric forecast
bash scripts/weather.sh
```

## Output example

```json
{
  "target": "2026-04-14 09:00 UTC",
  "telemetry": {
    "wave_height_m": 1.4,
    "swell_primary_m": 1.1,
    "wind_waves_m": 0.5,
    "peak_period_s": 8.2,
    "wave_direction_deg": 215,
    "water_temp_c": 14.3,
    "current_speed_mps": 0.18,
    "current_direction_deg": 95,
    "salinity_psu": 34.8,
    "chlorophyll_mgm3": 0.65,
    "kd": 0.28
  },
  "tide": {
    "state": "rising",
    "height_m": 1.2,
    "next_event": "high_tide",
    "hours_to_next": 2.5
  },
  "surf": {
    "score": 75,
    "rating": "good"
  },
  "dive": {
    "score": 68,
    "rating": "good",
    "visibility_m": 6.1,
    "factors": {
      "visibility": 35,
      "currents": 22,
      "waves": 18,
      "temperature": 8,
      "salinity": 9
    }
  },
  "trend": {
    "6h_direction": "worsening",
    "expected_change": -5
  },
  "safety_alerts": [
    {
      "type": "cold_water",
      "level": "warning",
      "message": "14.3°C — recommend wetsuit"
    }
  ]
}
```
