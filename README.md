# diving-forecast

Oceanographic forecast engine for surf and dive conditions on the Galician coast (A Coruña area). Uses real-time data from [Copernicus Marine Service](https://marine.copernicus.eu/) IBI (Iberian-Biscay-Irish) regional models.

## Scripts

| Script | Description |
|--------|-------------|
| `oceanographic_engine.py` | Main engine — fetches wave, current, temperature and chlorophyll data from Copernicus and computes surf/dive scores |
| `ocean_engine.py` | Lightweight analyzer — takes JSON input and scores per-beach conditions for beginner/returning surfers |
| `surf_copernicus.py` | Snapshot fetcher — pulls current conditions and writes to local cache |
| `surf.sh` | Shell wrapper for `oceanographic_engine.py` |
| `weather.sh` | Atmospheric forecast via wttr.in (no API key needed) |

## Models used

- **Waves**: `cmems_mod_ibi_wav_anfc_0.027deg_PT1H-i` — significant wave height, peak period, mean direction
- **Physics**: `cmems_mod_ibi_phy_anfc_0.027deg-3D_PT1H-m` — currents (U/V), water temperature
- **Biogeochemistry**: `cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1D-m` — chlorophyll-a

## Dive visibility model

Underwater visibility is estimated using a multi-factor attenuation model for Case 2 coastal waters:

```
Kd_total = Kd_baseline + Kd_chl(chl) + Kd_waves(effective_wave) + Kd_currents(speed)
visibility_m = 1.4 / Kd_total
```

`effective_wave` is an exponentially weighted average of wave height over the past 5 days (half-life = 24h), modelling sediment settling after storm events.

## Setup

```bash
pip install -r requirements.txt
```

Copernicus Marine credentials are needed. Configure them once with:

```bash
copernicusmarine login
# or set env vars:
# COPERNICUSMARINE_SERVICE_USERNAME / COPERNICUSMARINE_SERVICE_PASSWORD
```

## Usage

```bash
# Today at 9:00 UTC
python3 scripts/oceanographic_engine.py

# Tomorrow at 7:00 UTC
python3 scripts/oceanographic_engine.py tomorrow 7

# Atmospheric forecast
bash scripts/weather.sh
```

Output example:
```json
{
  "target": "2026-04-14 09:00 UTC",
  "telemetry": {
    "wave_height_m": 1.2,
    "water_temp_c": 14.3,
    "chlorophyll": 0.8,
    "current_speed_mps": 0.12,
    "effective_wave_m": 1.05
  },
  "surf": { "score": 72, "rating": "good" },
  "dive": { "score": 68, "visibility_m": 4.2, "rating": "good" }
}
```
