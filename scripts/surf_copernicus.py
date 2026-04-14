import copernicusmarine
import pandas as pd
import numpy as np
import json
import os
import xarray as xr
from datetime import datetime, timedelta, timezone

DS_WAV = "cmems_mod_ibi_wav_anfc_0.027deg_PT1H-i"
DS_PHY = "cmems_mod_ibi_phy_anfc_0.027deg-3D_PT1H-m"
DS_BGC = "cmems_mod_ibi_bgc_anfc_0.027deg-3D_P1D-m"

LAT, LON = 43.382167, -8.389000

CACHE_FILE = os.path.join(os.path.dirname(__file__), "cache", "copernicus.json")


def fetch_all():
    now = datetime.now(timezone.utc)
    R = 0.1  # ~10km radius

    def get_var(ds_id, vars):
        try:
            ds = copernicusmarine.open_dataset(
                dataset_id=ds_id,
                minimum_longitude=LON-R, maximum_longitude=LON+R,
                minimum_latitude=LAT-R, maximum_latitude=LAT+R,
                start_datetime=(now - timedelta(hours=6)).strftime("%Y-%m-%dT%H:%M:%S"),
                end_datetime=(now + timedelta(hours=6)).strftime("%Y-%m-%dT%H:%M:%S")
            )
            point_ds = ds.sel(time=now, method='nearest').load()
            for v in vars:
                if v in point_ds:
                    data = point_ds[v].values
                    if len(data.shape) > 0:
                        valid = data[~np.isnan(data)]
                        if len(valid) > 0:
                            return float(valid[0])
            return 0
        except:
            return 0

    h = get_var(DS_WAV, ["VHM0"])
    p = get_var(DS_WAV, ["VTPK"])
    d = get_var(DS_WAV, ["VMDR"])
    t = get_var(DS_PHY, ["thetao"])
    uo = get_var(DS_PHY, ["uo"])
    vo = get_var(DS_PHY, ["vo"])
    chl = get_var(DS_BGC, ["chl"])

    res = {
        "wave_height": h, "wave_period": p, "wave_dir": d,
        "temp": t, "chl": chl, "current_speed": np.sqrt(uo**2 + vo**2)
    }

    os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)
    with open(CACHE_FILE, 'w') as f:
        json.dump([res], f)
    print("Ocean telemetry captured.")


if __name__ == "__main__":
    fetch_all()
