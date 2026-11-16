import truststore

# Use Windows' native certificate trust store
truststore.inject_into_ssl()

import requests
import pandas as pd
import math
from datetime import timedelta
from itertools import combinations

"""*******************************************************************************************************************************************************
 Date__(UT)__HR:MN, , , R.A._(ICRF), DEC_(ICRF), R.A._(a-app), DEC_(a-app), Azi_(a-app), Elev_(a-app),    APmag,  S-brt,             delta,     deldot,
*******************************************************************************************************************************************************"""

cols = ["planet", "date", "flag1", "flag2",
        "ra_icrf", "dec_icrf", "ra_app", "dec_app",
        "azi", "elev", "apmag", "sbrt", "delta", "deldot", "extra"]
URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
JST_OFFSET = timedelta(hours=9)
TIME_FMT = "%Y-%m-%d %H:%M"
PLANETS = {"10":"Sun","199": "Mercury", "299": "Venus", "499": "Mars", "599": "Jupiter"}

BASE_PARAMS = {
    "format": "json",
    "OBJ_DATA": "YES",
    "MAKE_EPHEM": "YES",
    "EPHEM_TYPE": "OBSERVER",
    "CENTER": "coord@399",
    "COORD_TYPE": "GEODETIC",
    "SITE_COORD": "130.87521,33.88342,0",
    "QUANTITIES": "1,2,4,9,20",
    "CSV_FORMAT": "YES",
    "ANG_FORMAT": "DEG",
}
FMT = "%Y-%b-%d %H:%M"

def local_to_utc_gz1(local_dt):
    #JST datetime -> UTC datetime (what Horizons expects)
    return local_dt - JST_OFFSET

def build_params(body, start_utc, stop_utc, step, **overrides):
    #Extra keywords override or add keys.
    params = dict(BASE_PARAMS)
    params.update({
        "COMMAND": body,
        "START_TIME": start_utc.strftime(TIME_FMT),
        "STOP_TIME": stop_utc.strftime(TIME_FMT),
        "STEP_SIZE": step,
    })
    params.update(overrides)
    return params
def fetch_rows(params):
    """One Horizons request -> list of raw CSV lines between $$SOE and $$EOE."""
    response = requests.get(URL, params=params, timeout=30)
    response.raise_for_status()
    result = response.json()["result"]
    if "$$SOE" not in result:
        raise ValueError(f"No ephemeris in response:\n{result[:500]}")
    start = result.index("$$SOE") + len("$$SOE")
    end = result.index("$$EOE")
    return result[start:end].strip().splitlines()

def JPL_request(start_utc, stop_utc, step,outfile, bodies=PLANETS,**overrides):
    lines = []
    for body in bodies:
        try:
            rows = fetch_rows(build_params(body, start_utc, stop_utc, step, **overrides))
        except (requests.RequestException, ValueError, KeyError) as e:
            print(f"Skipping {body}: {e}")
            continue
        lines += [f"{body},{row}" for row in rows]
        print(f"{body}: {len(rows)} rows")

    csv_data = "\n".join(lines) + "\n"
    with open(outfile, "w", encoding="utf-8") as f:
        f.write(csv_data)
    return csv_data





#Examples
#Wide Scan
#JPL_request(datetime(2026, 11, 8, 0, 0), datetime(2026, 11, 18, 23, 50))
#JPL_request(local_to_utc_gz1(datetime(2026, 11, 16, 5, 00)), local_to_utc_gz1(datetime(2026, 11, 16, 6, 30)), step="1 min",outfile="event_ephemeris_candidate_gz1.csv",SKIP_DAYLT="YES")
#angular_separation_analysis("event_ephemeris_candidate_gz1.csv","angular_separation_candidate_gz1.csv")   
#ts = per_timestamp("angular_separation_candidate_gz1.csv")
#nights = best_nights(ts,min_elev=5, max_sun=-12)
#print(nights)
#night = nights.index[0]            # a datetime.date, the best night
#d, b = best_time_window(ts, night, metric="S", tol=100,min_elev=5, max_sun=-12,max_gap_h=1.5)                    # worst-pair spread
#print(d)
#print(b)
#PEAK = d.loc[d["S"].idxmin(), "dt"]
#print("this is peak:",PEAK)
#print(len(d))
#print(d[["dt", "S", "S_mean", "S_min", "min_elev", "sun_elev"]].to_string())