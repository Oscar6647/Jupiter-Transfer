import numpy as np, matplotlib, pandas as pd
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from scipy.optimize import brentq
from Orbital_analyzer import kepler_state, K
from Jupiter import NAMES, COLORS


launch_days = np.arange(-100, 700, 7)
arrive_days = np.arange(300, 365*5, 10)  
YR = 365.25
MU = K**2 * YR**2                    # AU^3/yr^2 (differs from 4*pi^2 by about 2e-5)
ASSUMED_R, ASSUMED_PHI = {399: 1.0, 599: 5.2}, 120.0

EL  = {int(b): g.iloc[0] for b, g in pd.read_csv("planetary_orbital_elements.csv").groupby("Body")}
JD0 = EL[399].JDTDB                  # departure epoch ("now")

KMS = 4.74047
C3  = np.full((len(arrive_days), len(launch_days)), np.nan)
VIN = C3.copy()

ASSUMED = True      # True -> textbook values, False -> real JPL values
TOF     = 1.6      # years

# ---------- profile 

def state(body,jd): 
    p, v = kepler_state(EL[body], jd)
    return p[0], v[0] * YR

def radii(assumed):
    return ASSUMED_R if assumed else {b: EL[b].A for b in (399, 599)}

def setup(assumed, tof):
    rE, _ = state(399, JD0)
    thE = np.arctan2(rE[1], rE[0])

    if assumed:
        R1, R2 = ASSUMED_R[399], ASSUMED_R[599]
        nJ = np.sqrt(MU / R2**3)
        thJ = thE + np.radians(ASSUMED_PHI) + nJ*tof
        r1 = R1*np.array([
            np.cos(thE), np.sin(thE), 0.0
        ])
        r2 = R2*np.array([
            np.cos(thJ), np.sin(thJ), 0.0
        ])
        # Circular Earth velocity
        vc = np.sqrt(MU/R1)
        vE = vc*np.array([
            -np.sin(thE), np.cos(thE), 0.0
        ])
    else:
        r1 = rE
        r2, _ = state(599, JD0 + tof*YR)
        _, vE = state(399, JD0)
    return r1, r2, vE

# ---------- Lambert (universal variables, Bate-Mueller-White)
def stumpff(z):
    if z>1e-8: s=np.sqrt(z); return (s-np.sin(s))/s**3,(1-np.cos(s))/z
    if z<-1e-8: s=np.sqrt(-z); return (np.sinh(s)-s)/s**3,(np.cosh(s)-1)/(-z)
    return 1/6,0.5

def lambert(r1,r2,dt):
    n1,n2=np.linalg.norm(r1),np.linalg.norm(r2)
    c=np.dot(r1,r2)/(n1*n2); dnu=np.arccos(np.clip(c,-1,1))
    if np.cross(r1,r2)[2]<0: dnu=2*np.pi-dnu
    A=np.sin(dnu)*np.sqrt(n1*n2/(1-np.cos(dnu)))
    def tof(z):
        S,C=stumpff(z); y=n1+n2+A*(z*S-1)/np.sqrt(C)
        if y<0: return -1e9
        return ((y/C)**1.5*S + A*np.sqrt(y))/np.sqrt(MU)-dt
    try: z=brentq(tof,-30,4*np.pi**2-1e-3,xtol=1e-10)
    except ValueError: return None
    S,C=stumpff(z); y=n1+n2+A*(z*S-1)/np.sqrt(C)
    f=1-y/n1; g=A*np.sqrt(y/MU); gd=1-y/n2
    return (r2-f*r1)/g, (gd*r2-r1)/g


def transfer(assumed, tof=1.6):
    R = radii(assumed); R1, R2 = R[399], R[599]

    # ---- Hohmann
    a_h = (R1 + R2)/2
    e_h = (R2 - R1)/(R2 + R1)
    tof_h = np.pi*np.sqrt(a_h**3/MU)
    nE, nJ = np.sqrt(MU/R1**3), np.sqrt(MU/R2**3)
    phi_h = np.degrees(np.pi - nJ*tof_h)
    syn = 2*np.pi/(nE - nJ)
    print(f"Hohmann: a={a_h:.3f} e={e_h:.3f} TOF={tof_h:.2f} yr phi={phi_h:.1f} deg synodic={syn:.3f} yr")

    # ---- chosen-TOF transfer
    r1, r2, vE = setup(assumed, tof)
    n1, n2 = np.linalg.norm(r1), np.linalg.norm(r2)
    c = np.linalg.norm(r2 - r1); s = (n1 + n2 + c)/2
    dth = np.arctan2(np.cross(r1, r2)[2], r1 @ r2) % (2*np.pi)    # prograde transfer angle

    a_min = s/2
    beta = 2*np.arcsin(np.sqrt((s - c)/s))
    sgn = -1 if dth < np.pi else 1                                # sign flips past 180 deg
    tof_min = np.sqrt(a_min**3/MU)*(np.pi + sgn*(beta - np.sin(beta)))

    v1, v2 = lambert(r1, r2, tof)
    a = 1/(2/n1 - (v1 @ v1)/MU)
    ev = np.cross(v1, np.cross(r1, v1))/MU - r1/n1
    e = np.linalg.norm(ev)
    nu_dep = (np.arctan2(r1[1], r1[0]) - np.arctan2(ev[1], ev[0])) % (2*np.pi)
    dv = np.linalg.norm(v1 - vE)                                  # AU/yr, times 4.74047 for km/s
    print("\nMINIMUM-ENERGY TRANSFER")
    print(f"a_min = {a_min:.4f} AU")
    print(f"TOF_min = {tof_min:.4f} years")

    print("\nSPECIFIED-TIME TRANSFER")
    print(f"TOF = {tof:.4f} years")
    print(f"a = {a:.4f} AU")
    print(f"e = {e:.4f}")
    print(f"Departure delta-v = {dv*4.74047:.3f} km/s")
    return dict(a=a, e=e, nu_dep=nu_dep, dtheta=dth, r1=r1, r2=r2, tof=tof)

# ---------- porkchop: launch vs arrival

def launch_window():
    KMS = 4.74047
    launch_days = np.arange(0, 700, 7)               # days after JD0
    arrive_days = np.arange(365, 365*5, 10)

    t0 = datetime(2000, 1, 1, 12) + timedelta(days=JD0 - 2451545.0)
    launch = [t0 + timedelta(days=float(d)) for d in launch_days]
    arrive = [t0 + timedelta(days=float(d)) for d in arrive_days]

    nE, nJ = np.sqrt(MU/1.0**3), np.sqrt(MU/5.2**3)
    phi_now = np.radians(82.4)
    phi_h   = np.radians(97.1)
    t_wait  = ((phi_now - phi_h) % (2*np.pi)) / (nE - nJ)
    print(f"Hohmann window opens in {t_wait:.2f} yr -> {t_wait*365.25:.0f} days after JD0")

    C3  = np.full((len(arrive_days), len(launch_days)), np.nan)
    VIN = C3.copy()

    for i, dl in enumerate(launch_days):
        rE, vE = state(399, JD0 + dl)                        # Earth at launch
        for j, da in enumerate(arrive_days):
            tof = (da - dl) / YR
            if tof < 0.5: continue
            rJ, vJ = state(599, JD0 + da)                    # Jupiter at arrival
            dth = np.arctan2(np.cross(rE, rJ)[2], rE @ rJ) % (2*np.pi)
            if abs(dth - np.pi) < np.radians(1.0): continue
            sol = lambert(rE, rJ, tof)
            if sol is None: continue
            v1, v2 = sol
            C3[j, i]  = (np.linalg.norm(v1 - vE) * KMS)**2
            VIN[j, i] = np.linalg.norm(v2 - vJ) * KMS

    print("valid points:", np.count_nonzero(~np.isnan(C3)), "of", C3.size)   # sanity check

    k = np.nanargmin(C3); j, i = np.unravel_index(k, C3.shape)
    print(f"Best: launch {launch[i].date()} arrive {arrive[j].date()} "
          f"TOF {(arrive[j]-launch[i]).days} d  C3={C3[j,i]:.1f} km2/s2  v_inf_arr={VIN[j,i]:.1f} km/s")

    Cp = np.where(C3 < 300, C3, np.nan)
    fig, ax = plt.subplots(figsize=(9, 6))
    cs = ax.contourf(launch, arrive, Cp, levels=np.arange(60, 301, 10), cmap="viridis_r")
    t_h = t_wait * 365.25                       # days after JD0
    ax.plot(t0 + timedelta(days=t_h), t0 + timedelta(days=t_h + 2.73*YR), "wo", mec="k", label="Hohmann (TOF 2.73 yr)")
    ax.plot(launch[0], t0 + timedelta(days=1.6*YR), "m^", label="depart now, TOF 1.6 yr")
    ax.contour(launch, arrive, Cp, levels=[80, 100, 150, 200], colors="w", linewidths=.6)
    ax.plot(launch[i], arrive[j], "r*", ms=14, label=f"min C3 = {C3[j,i]:.0f} km²/s²")
    plt.colorbar(cs, label="launch C3 (km²/s²)")
    ax.set_xlabel("launch date"); ax.set_ylabel("Jupiter arrival date")
    ax.set_title("Porkchop plot: Earth → Jupiter")
    ax.legend(); fig.autofmt_xdate()
    plt.tight_layout(); plt.savefig("jupiter_porkchop.png", dpi=150); plt.close()

def plot_transfer(res, assumed, hohmann=True):
    r1, r2 = res["r1"][:2], res["r2"][:2]
    a, e, nu_dep, dth, tof = res["a"], res["e"], res["nu_dep"], res["dtheta"], res["tof"]
    th_dep = np.arctan2(r1[1], r1[0])

    # ---- transfer conic: focus at the Sun, rotated so the departure point lies on it
    def conic(a, e, omega, nu):
        r = a*(1 - e**2) / (1 + e*np.cos(nu))
        return np.c_[r*np.cos(nu + omega), r*np.sin(nu + omega)]

    omega = th_dep - nu_dep
    nu_dep = (nu_dep + np.pi) % (2*np.pi) - np.pi          # map to (-pi, pi]; matters for hyperbolas

    if e < 1:
        nu_full = np.linspace(0, 2*np.pi, 600)
    else:
        nu_inf = np.arccos(-1/e)                           # asymptote angle
        nu_full = np.linspace(-0.97*nu_inf, 0.97*nu_inf, 600)

    omega = th_dep - nu_dep
    full = conic(a, e, omega, nu_full)
    arc  = conic(a, e, omega, np.linspace(nu_dep, nu_dep + dth, 300))

    # check: the arc must end on Jupiter's arrival position
    miss = np.linalg.norm(conic(a, e, omega, np.array([nu_dep + dth]))[0] - r2)
    print(f"arrival miss distance: {miss:.2e} AU")

    # ---- planet orbits and Jupiter at departure
    fig, ax = plt.subplots(figsize=(8, 8))
    if assumed:
        th = np.linspace(0, 2*np.pi, 400)
        for body, R in ASSUMED_R.items():
            ax.plot(R*np.cos(th), R*np.sin(th), color=COLORS[body], lw=0.8,
                    label=f"{NAMES[body]} orbit (R = {R} AU)")
        thJ0 = th_dep + np.radians(ASSUMED_PHI)
        pJ0 = ASSUMED_R[599] * np.array([np.cos(thJ0), np.sin(thJ0)])
    else:
        for body in (399, 599):
            period = 2*np.pi * EL[body].A**1.5 / K                      # days
            orb, _ = kepler_state(EL[body], JD0 + np.linspace(0, period, 600))
            ax.plot(orb[:, 0], orb[:, 1], color=COLORS[body], lw=0.8,
                    label=f"{NAMES[body]} orbit (JPL)")
        pJ0 = state(599, JD0)[0][:2]

    # ---- Jupiter's path during the flight (departure -> arrival)
    ang0 = np.arctan2(pJ0[1], pJ0[0]); ang1 = np.arctan2(r2[1], r2[0])
    sweep = (ang1 - ang0) % (2*np.pi)
    angs = ang0 + np.linspace(0, sweep, 100)
    RJ = np.linalg.norm(pJ0)
    ax.plot(RJ*np.cos(angs), RJ*np.sin(angs), color="tab:orange", lw=3, alpha=0.4,
            label="Jupiter during flight")

    # ---- optional Hohmann comparison
    if hohmann:
        R = radii(assumed); R1, R2 = R[399], R[599]
        a_h, e_h = (R1 + R2)/2, (R2 - R1)/(R2 + R1)
        hz = conic(a_h, e_h, th_dep, np.linspace(0, np.pi, 300))
        ax.plot(*hz.T, ":", color="0.4", lw=1.2, label="Hohmann (for comparison)")

    # ---- markers
    ax.plot(*full.T, "--", color="0.6", lw=0.8, label="transfer ellipse")
    ax.plot(*arc.T, color="tab:red", lw=2, label=f"transfer arc (TOF = {tof:.2f} yr)")
    ax.scatter(0, 0, c="gold", s=120, edgecolors="k", zorder=5, label="Sun")
    ax.scatter(*r1, c="tab:green", s=70, edgecolors="k", zorder=5, label="Earth at departure")
    ax.scatter(*pJ0, c="tab:orange", s=40, marker="s", edgecolors="k", zorder=5, label="Jupiter at departure")
    ax.scatter(*r2, c="tab:orange", s=90, edgecolors="k", zorder=5, label="Jupiter at arrival")

    lim = 1.2 * max(5.5, np.abs(arc).max())
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect("equal"); ax.grid(alpha=0.3)
    ax.set_xlabel("X (AU)"); ax.set_ylabel("Y (AU)")
    ax.set_title(f"Earth → Jupiter transfer, TOF = {tof} yr ({'assumed' if assumed else 'real JPL'})")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    return fig


def run(assumed=ASSUMED, tof=TOF, save=True):
    res = transfer(assumed=assumed, tof=tof)          # prints Hohmann, a_min, TOF_min, a, e, dv
    #fig = plot_transfer(res, assumed)
    #if save:
        #fig.savefig(f"transfer_{'assumed' if assumed else 'real'}_{tof}yr.png", dpi=150)
    #return res, fig

if __name__ == "__main__":
    run()
    #plt.show()
