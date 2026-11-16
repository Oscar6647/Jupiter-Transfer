// physics.js - line-by-line JS port of Orbital_analyzer.kepler_state and Transfer.state / setup / lambert.
// Constants and elements come from mission_data.js (DATA), which export_web_data.py writes from the Python modules.
(function (G) {
  const D = G.DATA, { K, MU, YR, KMS, JD0 } = D, D2R = Math.PI / 180, TAU = 2 * Math.PI;
  const nrm = v => Math.hypot(v[0], v[1], v[2]), dot = (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
  const sub = (p, q) => [p[0] - q[0], p[1] - q[1], p[2] - q[2]], scl = (p, k) => [p[0] * k, p[1] * k, p[2] * k];
  const cross = (p, q) => [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]];

  // Orbital_analyzer.solve_kepler
  function solveKepler(M, e) {
    M = ((M % TAU) + TAU) % TAU; let E = e < 0.8 ? M : Math.PI;
    for (let i = 0; i < 100; i++) { const d = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E)); E -= d; if (Math.abs(d) < 1e-12) break; }
    return E;
  }
  // Orbital_analyzer.kepler_state  (+ Transfer.state: velocity converted to AU/yr)
  function keplerState(el, jd) {
    const a = el.A, e = el.EC, inc = el.IN * D2R, Om = el.OM * D2R, w = el.W * D2R, n = K / a ** 1.5;
    const E = solveKepler(el.MA * D2R + n * (jd - el.JDTDB), e);
    const x = a * (Math.cos(E) - e), y = a * Math.sqrt(1 - e * e) * Math.sin(E), f = n / (1 - e * Math.cos(E));
    const vx = -a * Math.sin(E) * f, vy = a * Math.sqrt(1 - e * e) * Math.cos(E) * f;
    const cO = Math.cos(Om), sO = Math.sin(Om), ci = Math.cos(inc), si = Math.sin(inc), cw = Math.cos(w), sw = Math.sin(w);
    const R = [[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci], [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci], [sw * si, cw * si]];
    return [R.map(r => r[0] * x + r[1] * y), R.map(r => (r[0] * vx + r[1] * vy) * YR)];
  }
  // Transfer.setup(assumed=True): circular coplanar orbits, Jupiter ASSUMED_PHI ahead of Earth at JD0
  const thE0 = (() => { const p = keplerState(D.EL[399], JD0)[0]; return Math.atan2(p[1], p[0]); })();
  const TBR = { 399: D.assumed.R[399], 599: D.assumed.R[599] }, TBT = { 399: thE0, 599: thE0 + D.assumed.phi * D2R };
  function circState(b, jd) {
    const R = TBR[b], n = Math.sqrt(MU / R ** 3), t = TBT[b] + n * (jd - JD0) / YR;
    return [[R * Math.cos(t), R * Math.sin(t), 0], [-R * n * Math.sin(t), R * n * Math.cos(t), 0]];
  }
  const state = (b, jd, textbook) => textbook ? circState(b, jd) : keplerState(D.EL[b], jd);

  // Transfer.stumpff / Transfer.lambert  (brentq replaced by bisection on the same bracket)
  function stumpff(z) {
    if (z > 1e-8) { const s = Math.sqrt(z); return [(s - Math.sin(s)) / s ** 3, (1 - Math.cos(s)) / z]; }
    if (z < -1e-8) { const s = Math.sqrt(-z); return [(Math.sinh(s) - s) / s ** 3, (Math.cosh(s) - 1) / (-z)]; }
    return [1 / 6, 0.5];
  }
  function lambert(r1, r2, dt) {
    const n1 = nrm(r1), n2 = nrm(r2);
    let dnu = Math.acos(Math.max(-1, Math.min(1, dot(r1, r2) / (n1 * n2))));
    if (cross(r1, r2)[2] < 0) dnu = TAU - dnu;                       // prograde transfer
    const A = Math.sin(dnu) * Math.sqrt(n1 * n2 / (1 - Math.cos(dnu)));
    const Y = z => { const [S, C] = stumpff(z); return n1 + n2 + A * (z * S - 1) / Math.sqrt(C); };
    const tof = z => { const [S, C] = stumpff(z), y = Y(z); if (y < 0) return -1e9; return ((y / C) ** 1.5 * S + A * Math.sqrt(y)) / Math.sqrt(MU) - dt; };
    let lo = -30, hi = 4 * Math.PI ** 2 - 1e-3, fl = tof(lo), fh = tof(hi);
    if (!(fl * fh < 0)) return null;
    for (let k = 0; k < 200 && hi - lo > 1e-12; k++) { const m = (lo + hi) / 2, fm = tof(m); if (fl * fm <= 0) hi = m; else { lo = m; fl = fm; } }
    const z = (lo + hi) / 2, y = Y(z), f = 1 - y / n1, g = A * Math.sqrt(y / MU), gd = 1 - y / n2;
    return [scl(sub(r2, scl(r1, f)), 1 / g), scl(sub(scl(r2, gd), r1), 1 / g)];
  }
  // two-body RK4 (AU, yr)
  function propagate(r, v, T, N) {
    const h = T / N, P = [], V = []; let s = [...r, ...v];
    const F = s => { const q = Math.hypot(s[0], s[1], s[2]), k = -MU / q ** 3; return [s[3], s[4], s[5], k * s[0], k * s[1], k * s[2]]; };
    for (let k = 0; ; k++) {
      P.push(s.slice(0, 3)); V.push(s.slice(3)); if (k === N) break;
      const k1 = F(s), k2 = F(s.map((x, i) => x + h / 2 * k1[i])), k3 = F(s.map((x, i) => x + h / 2 * k2[i])), k4 = F(s.map((x, i) => x + h * k3[i]));
      s = s.map((x, i) => x + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
    }
    return { P, V };
  }
  // conic elements of a state, as in Transfer.transfer()
  function elements(r, v) {
    const n = nrm(r), a = 1 / (2 / n - dot(v, v) / MU), ev = sub(scl(cross(v, cross(r, v)), 1 / MU), scl(r, 1 / n));
    return { a, e: nrm(ev), ev };
  }
  // Hohmann numbers, as in Transfer.transfer()
  function hohmann(R1, R2) {
    const a = (R1 + R2) / 2, e = (R2 - R1) / (R2 + R1), tof = Math.PI * Math.sqrt(a ** 3 / MU);
    const nE = Math.sqrt(MU / R1 ** 3), nJ = Math.sqrt(MU / R2 ** 3);
    return { a, e, tof, phi: (Math.PI - nJ * tof) / D2R, syn: TAU / (nE - nJ), nE, nJ, vp: Math.sqrt(MU * (2 / R1 - 1 / a)) };
  }
  G.PHYS = { D2R, TAU, K, MU, YR, KMS, JD0, nrm, dot, sub, scl, cross, solveKepler, keplerState, circState, state, stumpff, lambert, propagate, elements, hohmann };
})(typeof window !== "undefined" ? window : globalThis);
