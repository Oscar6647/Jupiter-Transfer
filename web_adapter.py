import sys, types, io, json, contextlib
from datetime import datetime, timedelta
import numpy as np

# --- Orbital_analyzer imports Horizon (the JPL request helper). Only kepler_state/K are used
#     here, so if Horizon.py isn't next to this file a stand-in is enough.
try:
    import Horizon  
except ImportError:
    sys.modules["Horizon"] = types.SimpleNamespace(JPL_request=None)

# --- Transfer.py does `from Jupiter import NAMES, COLORS`. Importing the real Jupiter.py runs its
#     plots at import time, so provide just those two dicts (same values as Jupiter.py).
sys.modules["Jupiter"] = types.SimpleNamespace(
    NAMES={399: "Earth", 599: "Jupiter"},
    COLORS={399: "tab:green", 599: "tab:orange"},
)

import matplotlib
matplotlib.use("Agg", force=True)
import transfer as T                      # <- your module: state, lambert, transfer, MU, EL, JD0 ...
from Orbital_analyzer import kepler_state, K

YR, KMS, MU, EL, JD0 = T.YR, T.KMS, T.MU, T.EL, T.JD0


def captured(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = fn(*a, **k)
    return out, buf.getvalue()


def jd2iso(jd):
    return (datetime(2000, 1, 1, 12) + timedelta(days=float(jd) - 2451545.0)).strftime("%Y-%m-%d")


def lst(v, nd=6):
    return [round(float(x), nd) for x in np.ravel(v)]


# ---------------------------------------------------------------- 1. Hohmann (textbook + real a)
def hohmann(R1, R2):
    a_h = (R1 + R2) / 2
    e_h = (R2 - R1) / (R2 + R1)
    tof_h = np.pi * np.sqrt(a_h**3 / MU)
    nE, nJ = np.sqrt(MU / R1**3), np.sqrt(MU / R2**3)
    phi_h = np.degrees(np.pi - nJ * tof_h)
    syn = 2 * np.pi / (nE - nJ)
    vE = np.sqrt(MU / R1); vp = np.sqrt(MU * (2 / R1 - 1 / a_h))
    vJ = np.sqrt(MU / R2); va = np.sqrt(MU * (2 / R2 - 1 / a_h))
    return dict(a=a_h, e=e_h, tof=tof_h, phi=phi_h, syn=syn,
                dv1=(vp - vE) * KMS, dv2=(vJ - va) * KMS, C3=((vp - vE) * KMS) ** 2)


hoh_text = hohmann(1.0, 5.2)
hoh_real = hohmann(EL[399].A, EL[599].A)

# ---------------------------------------------------------------- 2. transfer() runs, printed output kept
runs = {}
for assumed in (True, False):
    for tof in (1.0, 1.6, 2.0, 2.73):
        try:
            res, text = captured(T.transfer, assumed, tof)
        except Exception as ex:                      # Lambert can fail for odd geometries
            res, text = None, f"failed: {ex}"
        key = f"{'assumed' if assumed else 'real'}_{tof}"
        if res:
            runs[key] = dict(a=res["a"], e=res["e"], dtheta=np.degrees(res["dtheta"]),
                             r1=lst(res["r1"]), r2=lst(res["r2"]), tof=tof, log=text.strip())
        else:
            runs[key] = dict(log=text)

# ---------------------------------------------------------------- 3. real phase angle "now"
rE0, vE0 = T.state(399, JD0)
rJ0, vJ0 = T.state(599, JD0)
phase_now = np.degrees(np.arctan2(rJ0[1], rJ0[0]) - np.arctan2(rE0[1], rE0[0])) % 360
nE, nJ = np.sqrt(MU / EL[399].A**3), np.sqrt(MU / EL[599].A**3)
t_wait = np.radians((phase_now - hoh_real["phi"]) % 360) / (nE - nJ)       # years

# ---------------------------------------------------------------- 4. porkchop (same grid as launch_window)
launch_days = np.arange(0, 700, 7)
arrive_days = np.arange(365, 365 * 5, 10)
C3 = np.full((len(arrive_days), len(launch_days)), np.nan)
VIN = C3.copy()
for i, dl in enumerate(launch_days):
    rE, vE = T.state(399, JD0 + dl)
    for j, da in enumerate(arrive_days):
        tof = (da - dl) / YR
        if tof < 0.5:
            continue
        rJ, vJ = T.state(599, JD0 + da)
        dth = np.arctan2(np.cross(rE, rJ)[2], rE @ rJ) % (2 * np.pi)
        if abs(dth - np.pi) < np.radians(1.0):
            continue
        sol = T.lambert(rE, rJ, tof)
        if sol is None:
            continue
        v1, v2 = sol
        C3[j, i] = (np.linalg.norm(v1 - vE) * KMS) ** 2
        VIN[j, i] = np.linalg.norm(v2 - vJ) * KMS

j, i = np.unravel_index(np.nanargmin(C3), C3.shape)
dl, da = launch_days[i], arrive_days[j]
rE, vE = T.state(399, JD0 + dl); rJ, vJ = T.state(599, JD0 + da)
v1, v2 = T.lambert(rE, rJ, (da - dl) / YR)
best = dict(i=int(i), j=int(j), jd_l=float(JD0 + dl), jd_a=float(JD0 + da),
            launch=jd2iso(JD0 + dl), arrive=jd2iso(JD0 + da), tof_d=int(da - dl),
            C3=float(C3[j, i]), vinf=float(VIN[j, i]),
            r1=lst(rE, 9), v1=lst(v1, 9), r2=lst(rJ, 9), v2=lst(v2, 9),
            vE=lst(vE, 9), vJ=lst(vJ, 9))
a_b = 1 / (2 / np.linalg.norm(rE) - v1 @ v1 / MU)
e_b = np.linalg.norm(np.cross(v1, np.cross(rE, v1)) / MU - rE / np.linalg.norm(rE))
best.update(a=float(a_b), e=float(e_b))

# ---------------------------------------------------------------- 5. ephemeris check table for the JS port
check = []
for d in (0, 200, 500, 1000, 1500):
    p3, _ = T.state(399, JD0 + d); p5, _ = T.state(599, JD0 + d)
    check.append(dict(d=d, earth=lst(p3, 9), jup=lst(p5, 9)))

el_out = {str(b): dict(A=float(e.A), EC=float(e.EC), IN=float(e.IN), OM=float(e.OM),
                       W=float(e.W), MA=float(e.MA), JDTDB=float(e.JDTDB)) for b, e in EL.items()}

data = dict(
    generated=datetime.now().strftime("%Y-%m-%d %H:%M"),
    K=K, MU=MU, YR=YR, KMS=KMS, JD0=JD0, epoch=jd2iso(JD0), EL=el_out,
    assumed=dict(R=T.ASSUMED_R, phi=T.ASSUMED_PHI),
    hohmann=dict(textbook=hoh_text, real=hoh_real),
    phase_now=phase_now, t_wait_yr=t_wait, window=jd2iso(JD0 + t_wait * YR), jd_window=JD0 + t_wait * YR,
    runs=runs,
    porkchop=dict(launch_days=launch_days.tolist(), arrive_days=arrive_days.tolist(),
                  C3=[[None if np.isnan(x) else round(float(x), 2) for x in row] for row in C3],
                  VIN=[[None if np.isnan(x) else round(float(x), 3) for x in row] for row in VIN],
                  valid=int(np.count_nonzero(~np.isnan(C3))), total=int(C3.size)),
    best=best, check=check,
)


def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, (np.integer,)): return int(o)
    return o


with open("data.js", "w", encoding="utf-8") as f:
    f.write("// generated by export_web_data.py from Transfer.py + Orbital_analyzer.py. Do not edit.\n")
    f.write("var DATA = " + json.dumps(clean(data), separators=(",", ":")) + ";\n")

print(f"epoch {jd2iso(JD0)}  phase now {phase_now:.1f} deg  Hohmann phi {hoh_real['phi']:.1f} deg")
print(f"Hohmann textbook: a={hoh_text['a']:.3f} e={hoh_text['e']:.3f} TOF={hoh_text['tof']:.2f} yr "
      f"phi={hoh_text['phi']:.1f} syn={hoh_text['syn']:.3f} dv1={hoh_text['dv1']:.2f} dv2={hoh_text['dv2']:.2f}")
print(f"Window opens {data['window']} ({t_wait*365.25:.0f} d)")
print(f"Porkchop valid {data['porkchop']['valid']}/{data['porkchop']['total']}")
print(f"Best: launch {best['launch']} arrive {best['arrive']} TOF {best['tof_d']} d "
      f"C3={best['C3']:.1f} vinf={best['vinf']:.2f}  a={a_b:.3f} e={e_b:.3f}")
for k, v in runs.items():
    print(k, {kk: round(vv, 4) for kk, vv in v.items() if isinstance(vv, float)})
print("wrote mission_data.js")
