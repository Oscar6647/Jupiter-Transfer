import numpy as np, pandas as pd, matplotlib,os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Orbital_analyzer import JPL_Orbital_Elements_request,JPL_Orbital_Position_request,kepler_state,K

elements_csv="planetary_orbital_elements.csv"
positions_csv= "planetary_orbital_positions.csv"

# ---- 1. Orbital elements JPL
if not os.path.exists(elements_csv):
    JPL_Orbital_Position_request()
    JPL_Orbital_Elements_request()
    df = pd.read_csv(elements_csv)
# ---- 2. plot
NAMES  = {399: "Earth",      599: "Jupiter"}
COLORS = {399: "tab:green",  599: "tab:orange"}

# assignment-constrained values
ASSUMED_R   = {399: 1.0, 599: 5.2}            # AU, circular coplanar orbits
ASSUMED_PHI = 120.0                           # deg, Jupiter ahead of Earth at departure

def plot_orbits(assumed):
    els = pd.read_csv(elements_csv)
    vec = pd.read_csv(positions_csv)

    fig = plt.figure(figsize=(14, 6.5))
    ax  = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122)
    lon = {}

    for body, g in els.groupby("Body"):
        body = int(body)
        el = g.iloc[0]                                    # departure epoch
        v = vec[vec.Body == body]
        lon[body] = np.arctan2(v.y.iloc[0], v.x.iloc[0])
        if not assumed:                                   # real values, drawn per body
            period = 2*np.pi * el.A**1.5 / K
            full, _ = kepler_state(el, el.JDTDB + np.linspace(0, period, 2000))
            ax.plot(*full.T, color=COLORS[body], lw=0.8, label=f"{NAMES[body]} orbit (JPL)")
            ax.scatter(v.x.iloc[0], v.y.iloc[0], v.z.iloc[0], c=COLORS[body],
                       s=50, edgecolors="k", label=f"{NAMES[body]} at departure (JPL)")
            r = np.sqrt(v.x**2 + v.y**2 + v.z**2)
            ax2.plot((v.JDTDB - v.JDTDB.iloc[0]) / 365.25, r,
                     color=COLORS[body], label=NAMES[body])
        
    if assumed:
        th = np.linspace(0, 2*np.pi, 400)
        for body, R in ASSUMED_R.items():
            ax.plot(R*np.cos(th), R*np.sin(th), 0, ls="--", color=COLORS[body], lw=1.0,
                    label=f"{NAMES[body]} circular (R = {R} AU)")
            ax2.axhline(R, ls="--", color=COLORS[body], lw=0.8)

        thE = lon[399]                          # Earth keeps its real departure longitude
        thJ = thE + np.radians(ASSUMED_PHI)     # Jupiter placed phi ahead of it
        for body, ang in [(399, thE), (599, thJ)]:
            R = ASSUMED_R[body]
            ax.scatter(R*np.cos(ang), R*np.sin(ang), 0, marker="X", s=90, c=COLORS[body],
                        edgecolors="k", label=f"{NAMES[body]} at departure (assumed)")
        print(f"Assumed phase: {ASSUMED_PHI:.1f} deg")
    else:    
        real_phase = np.degrees(lon[599] - lon[399]) % 360
        print(f"Real phase: Jupiter ahead of Earth by {real_phase:.1f} deg")

    ax.scatter([0], [0], [0], c="gold", s=80, edgecolors="k", label="Sun")
    ax.set_xlabel("X (au)"); ax.set_ylabel("Y (au)"); ax.set_zlabel("Z (au)")
    ax.set_title("Earth and Jupiter, heliocentric ecliptic (J2000)")
    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6); ax.set_zlim(-0.5, 0.5)
    ax.set_box_aspect((12, 12, 1))  
    ax.legend(fontsize=7, loc="upper left"); ax.view_init(22, -50)

    ax2.set_xlabel("years since departure"); ax2.set_ylabel("distance from Sun (au)")
    ax2.set_title("Heliocentric distance (JPL; dashed = assumed radius)" if assumed else "Heliocentric distance (JPL)"); 
    ax2.grid(alpha=.3); ax2.legend()
    fig.tight_layout()
    return fig

plot_orbits(assumed=True)
plt.savefig("earth_jupiter_textbook.png", dpi=150)
plt.show()

plot_orbits(assumed=False)
plt.savefig("earth_jupiter_real.png", dpi=150)
plt.show()