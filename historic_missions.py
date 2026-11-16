"""
historic_missions.py  -  real NASA Jupiter missions, analysed with the homework modules.

    python historic_missions.py          (downloads once into history_cache/, then reuses it)

Pipeline (all your code):
  Horizon.JPL_request / fetch_rows  -> real heliocentric state vectors of each spacecraft + planets
  Orbital_analyzer.indep_backup     -> osculating a, e, i, specific energy of the spacecraft every 2 days
  Orbital_analyzer.kepler_state     -> (not needed here: real JPL planet vectors are used directly)
  transfer.lambert                  -> "what would a DIRECT shot have cost on the same launch day?"
Writes history.js for index.html.
"""
import os, sys, types, json, math
from datetime import datetime, timedelta
import numpy as np, pandas as pd

import Horizon                                     # your JPL helper (also injects truststore)
from Horizon import JPL_request, fetch_rows, build_params

sys.modules.setdefault("Jupiter", types.SimpleNamespace(
    NAMES={399: "Earth", 599: "Jupiter"}, COLORS={399: "tab:green", 599: "tab:orange"}))
import matplotlib
from Orbital_analyzer import indep_backup
matplotlib.use("Agg", force=True)
import transfer as T                               # your Lambert solver + constants

CACHE = "history_cache"; os.makedirs(CACHE, exist_ok=True)
VEC = dict(EPHEM_TYPE="VECTORS", CENTER="500@10", REF_PLANE="ECLIPTIC", REF_SYSTEM="ICRF", VEC_TABLE="2",
           VEC_CORR="NONE", OUT_UNITS="AU-D", CSV_FORMAT="YES", OBJ_DATA="NO",
           QUANTITIES=None, SITE_COORD=None, ANG_FORMAT=None, CAL_FORMAT=None,
           ELEV_CUT=None, SKIP_DAYLT=None, EXTRA_PREC=None, COORD_TYPE=None)
AU_KM, DAY_S = 1.495978707e8, 86400.0
E_KMS2 = (AU_KM / DAY_S) ** 2                      # AU^2/day^2 -> km^2/s^2
MU_KM = {"299": 324859.0, "399": 398600.4, "599": 1.26687e8, "699": 3.7931e7}
RADIUS = {"299": 6052, "399": 6378, "599": 71492, "699": 60268}
PNAME = {"299": "Venus", "399": "Earth", "599": "Jupiter", "699": "Saturn"}

#          id     name           launch (UTC, just after JPL coverage starts)  stop         colour
MISSIONS = [
    ("-23", "Pioneer 10",   "1972-03-03 12:00", "1974-03-01", "#9fd3ff",
     "First spacecraft to Jupiter. Brute force: a direct shot on a huge launch, no gravity assists."),
    ("-24", "Pioneer 11",   "1973-04-06 12:00", "1976-01-01", "#7fb0ff",
     "Followed Pioneer 10, then used Jupiter's gravity to bend its path toward Saturn."),
    ("-32", "Voyager 2",    "1977-08-20 18:00", "1981-10-01", "#c6a0ff",
     "Grand Tour alignment: Jupiter flung it on to Saturn (and later Uranus and Neptune)."),
    ("-31", "Voyager 1",    "1977-09-05 18:00", "1981-01-01", "#ff9ad5",
     "Launched after Voyager 2 but on a faster path, so it reached Jupiter first. Jupiter pushed it onto an escape trajectory."),
    ("-77", "Galileo",      "1989-10-19 06:00", "1995-12-15", "#ffd36b",
     "Too heavy for a direct shot on its rocket, so it flew Venus, Earth, Earth (VEEGA) to reach Jupiter."),
    ("-55", "Ulysses",      "1990-10-07 00:00", "1994-01-01", "#7dffcf",
     "Used Jupiter to leave the plane of the planets entirely, to fly over the Sun's poles."),
    ("-82", "Cassini",      "1997-10-15 12:00", "2001-06-30", "#ffb08a",
     "Venus, Venus, Earth, then a Jupiter flyby on the way to Saturn."),
    ("-98", "New Horizons", "2006-01-20 00:00", "2007-12-31", "#ff6b6b",
     "Fastest launch ever. Reached Jupiter in 13 months, then used it to speed up toward Pluto."),
    ("-61", "Juno",         "2011-08-05 18:00", "2016-07-10", "#5cff9a",
     "Launched too weakly ON PURPOSE, fired its engine in deep space to come back to Earth, then used Earth to reach Jupiter."),
]
PL_START, PL_STOP, PL_STEP = datetime(1972, 3, 1), datetime(2016, 8, 1), 4      # planets every 4 d
SC_STEP = 2                                                                         # spacecraft every 2 d
FMT = "%Y-%m-%d %H:%M"


def jd2dt(jd): return datetime(2000, 1, 1, 12) + timedelta(days=float(jd) - 2451545.0)


def vectors(body, start, stop, step, tag, center="500@10"):
    """Your JPL_request, cached. Returns DataFrame Body,JDTDB,Date,x,y,z,vx,vy,vz."""
    path = os.path.join(CACHE, f"{tag}.csv")
    if not os.path.exists(path):
        JPL_request(start, stop, step, path, bodies={body: tag}, **{**VEC, "CENTER": center})
    rows = []
    for line in open(path, encoding="utf-8"):
        p = [s.strip() for s in line.split(",")]
        if len(p) < 9 or not p[1]:
            continue
        rows.append([p[0], float(p[1]), p[2]] + [float(v) for v in p[3:9]])
    if not rows:
        raise RuntimeError(f"no data for {tag}; delete {path} and retry")
    return pd.DataFrame(rows, columns=["Body", "JDTDB", "Date", "x", "y", "z", "vx", "vy", "vz"])


# ------------------------------------------------------------------ planets (one request each)
planets = {b: vectors(b, PL_START, PL_STOP, f"{PL_STEP} d", f"planet_{b}") for b in PNAME}
PJD0 = planets["399"].JDTDB.iloc[0]


def planet_state(b, jd):
    """Linear interpolation of the real JPL planet vectors (fine at 4-day spacing for a heliocentric view)."""
    df = planets[b]; f = (jd - PJD0) / PL_STEP; i = int(np.clip(np.floor(f), 0, len(df) - 2)); q = f - i
    a, c = df.iloc[i], df.iloc[i + 1]
    r = np.array([a.x + (c.x - a.x) * q, a.y + (c.y - a.y) * q, a.z + (c.z - a.z) * q])
    v = np.array([a.vx + (c.vx - a.vx) * q, a.vy + (c.vy - a.vy) * q, a.vz + (c.vz - a.vz) * q])
    return r, v


def closest_approach(sid, body, jd_guess):
    """Fine search with Horizons centred on the planet: exact flyby date and altitude."""
    t = jd2dt(jd_guess)
    df = vectors(sid, t - timedelta(days=3), t + timedelta(days=3), "10 m", f"ca_{sid}_{body}_{int(jd_guess)}",
                 center=f"500@{body}")
    r = np.sqrt(df.x**2 + df.y**2 + df.z**2) * AU_KM
    k = int(r.idxmin())
    # periapsis of the planet-centred osculating conic through the closest sample: rp = a(1-e), exact for 2-body
    R = df.loc[k, ["x", "y", "z"]].values.astype(float) * AU_KM
    V = df.loc[k, ["vx", "vy", "vz"]].values.astype(float) * AU_KM / DAY_S
    mu = MU_KM[body]; h = np.cross(R, V); ev = np.cross(V, h) / mu - R / np.linalg.norm(R)
    e = np.linalg.norm(ev); p = h @ h / mu; rp = p / (1 + e)
    return float(df.JDTDB[k]), float(min(rp, r[k])) - RADIUS[body]


def launch_vinf(sid, jd_launch):
    """Hyperbolic excess speed from a geocentric JPL state 5 days after launch: v_inf^2 = v^2 - 2 mu / r."""
    t = jd2dt(jd_launch + 5)
    df = vectors(sid, t, t + timedelta(hours=1), "1 h", f"geo_{sid}", center="500@399")
    r = math.sqrt(df.x[0]**2 + df.y[0]**2 + df.z[0]**2) * AU_KM
    v = math.sqrt(df.vx[0]**2 + df.vy[0]**2 + df.vz[0]**2) * AU_KM / DAY_S
    return math.sqrt(max(v * v - 2 * MU_KM["399"] / r, 0))


def direct_lambert(jd_launch):
    """transfer.lambert from Earth on the real launch day to Jupiter on every arrival day, keep the cheapest."""
    rE, vE = planet_state("399", jd_launch); best = (np.inf, None); curve = []
    for dd in range(300, 2400, 10):
        rJ, _ = planet_state("599", jd_launch + dd)
        dth = np.arctan2(np.cross(rE, rJ)[2], rE @ rJ) % (2 * np.pi)
        if abs(dth - np.pi) < np.radians(1): continue
        sol = T.lambert(rE, rJ, dd / T.YR)
        if sol is None: continue
        c3 = (np.linalg.norm(sol[0] - vE * T.YR) * T.KMS) ** 2
        curve.append([dd, round(float(c3), 1)])
        if c3 < best[0]: best = (c3, dd)
    return best[0], best[1], curve


out = []
for sid, name, start, stop, col, blurb in MISSIONS:
    print(f"\n=== {name} ({sid})")
    t0, t1 = datetime.strptime(start, FMT), datetime.strptime(stop, FMT[:8])
    sc = vectors(sid, t0, t1, f"{SC_STEP} d", f"sc_{sid}")

    # ---- your element derivation on the real spacecraft vectors
    base = os.path.join(CACHE, f"sc_{sid}_vec.csv"); elo = os.path.join(CACHE, f"sc_{sid}_elements.csv")
    sc.assign(Body=sid).to_csv(base, index=False)
    indep_backup(base, elo)
    el = pd.read_csv(elo)

    # ---- flybys: local minima of the distance to each planet
    events = []
    for b, lim in (("299", .08), ("399", .08), ("599", .5), ("699", .5)):
        P = np.array([planet_state(b, j)[0] for j in sc.JDTDB])
        d = np.linalg.norm(sc[["x", "y", "z"]].values - P, axis=1)
        for k in range(1, len(d) - 1):
            if d[k] < lim and d[k] <= d[k - 1] and d[k] <= d[k + 1] and sc.JDTDB[k] - sc.JDTDB[0] > 20:
                jd_ca, alt = closest_approach(sid, b, sc.JDTDB[k])
                w = 15 if b in ("299", "399") else 30                          # samples (2 d) away from the planet
                i0 = max(0, k - w); i1 = min(len(sc) - 1, k + w)
                g = lambda col, i: float(el[col][i])
                ev_ = dict(body=PNAME[b], jd=jd_ca, date=jd2dt(jd_ca).strftime("%Y-%m-%d"), alt_km=alt,
                           a0=g("Semimajor Axis", i0), a1=g("Semimajor Axis", i1), e0=g("Eccentricity", i0),
                           e1=g("Eccentricity", i1), i0=g("Inclination", i0), i1=g("Inclination", i1),
                           en0=g("Specific Energy", i0) * E_KMS2, en1=g("Specific Energy", i1) * E_KMS2)
                events.append(ev_)
                print(f"  {PNAME[b]:8s} {jd2dt(jd_ca):%Y-%m-%d}  altitude {alt:12,.0f} km   a {ev_['a0']:7.2f} -> {ev_['a1']:7.2f} AU"
                      f"   e {ev_['e0']:.3f} -> {ev_['e1']:.3f}   i {ev_['i0']:5.1f} -> {ev_['i1']:5.1f}   "
                      f"energy {ev_['en0']:7.1f} -> {ev_['en1']:7.1f} km2/s2")
    # ---- engine burns: energy jumps with no planet nearby (deep-space manoeuvres)
    E = el["Specific Energy"].values * E_KMS2
    near = np.min([np.linalg.norm(sc[["x", "y", "z"]].values - np.array([planet_state(b, j)[0] for j in sc.JDTDB]), axis=1)
                   for b in PNAME], axis=0)
    for k in range(10, len(E) - 1):
        if abs(E[k + 1] - E[k]) > 3 and near[k] > 0.3 and near[k + 1] > 0.3:
            jd_b = float(sc.JDTDB[k] + SC_STEP / 2)
            events.append(dict(body="Engine burn", jd=jd_b, date=jd2dt(jd_b).strftime("%Y-%m-%d"), alt_km=None,
                               a0=float(el["Semimajor Axis"][k]), a1=float(el["Semimajor Axis"][k + 1]),
                               e0=float(el["Eccentricity"][k]), e1=float(el["Eccentricity"][k + 1]),
                               i0=float(el["Inclination"][k]), i1=float(el["Inclination"][k + 1]),
                               en0=float(E[k]), en1=float(E[k + 1])))
            print(f"  engine burn ~{jd2dt(jd_b):%Y-%m-%d}: energy {E[k]:.1f} -> {E[k+1]:.1f} km2/s2")
    events.sort(key=lambda e: e["jd"])
    jup = next(e for e in events if e["body"] == "Jupiter")

    vinf = launch_vinf(sid, sc.JDTDB[0])
    dc3, dtof, dcurve = direct_lambert(sc.JDTDB[0])
    tof = jup["jd"] - sc.JDTDB[0]
    print(f"  launch C3 (JPL) {vinf**2:6.1f} km2/s2   time to Jupiter {tof:5.0f} d   "
          f"direct Lambert that day: C3 {dc3:6.1f} (TOF {dtof} d)")

    out.append(dict(
        id=sid, name=name, color=col, blurb=blurb, jd0=float(sc.JDTDB[0]), step=SC_STEP,
        launch=jd2dt(sc.JDTDB[0]).strftime("%Y-%m-%d"),
        sc=[round(float(v), 4) for v in sc[["x", "y", "z"]].values.ravel()],
        spd=[round(float(v) * AU_KM / DAY_S, 2) for v in np.sqrt(sc.vx**2 + sc.vy**2 + sc.vz**2)],
        en=[round(float(v) * E_KMS2, 1) for v in el["Specific Energy"]],
        a=[None if not np.isfinite(v) else round(float(v), 3) for v in el["Semimajor Axis"]],
        e=[round(float(v), 4) for v in el["Eccentricity"]],
        inc=[round(float(v), 2) for v in el["Inclination"]],
        events=events, direct_curve=dcurve, c3=vinf**2, vinf=vinf, tof_jup=tof, direct_c3=dc3, direct_tof=dtof,
    ))

# planets for the browser
pl = {b: [round(float(v), 4) for v in planets[b][["x", "y", "z"]].values.ravel()] for b in PNAME}
with open("history.js", "w", encoding="utf-8") as f:
    f.write("// generated by historic_missions.py (Horizon.py + Orbital_analyzer.indep_backup + transfer.lambert)\n")
    f.write("var HIST = " + json.dumps(dict(generated=datetime.now().strftime("%Y-%m-%d %H:%M"),
                                            pjd0=float(PJD0), pstep=PL_STEP, planets=pl, missions=out),
                                       separators=(",", ":")) + ";\n")
print("\nwrote history.js")
