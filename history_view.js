// history_view.js - "NASA missions" mode. Real JPL trajectories (history.js, made by historic_missions.py).
// Runs after the main script and reuses its scene, camera, helpers and panels.
var HV = (() => {
  const H = HIST, M = H.missions, st = { on: false, cur: -1, base: 0, span: 0, ev: 0, flashT: 0, flashTx: '', evShown: -1 };
  const hexRGBA = h => { const n = parseInt(h.slice(1), 16); return `rgba(${n >> 16},${n >> 8 & 255},${n & 255},1)`; };
  const PN = H.planets['399'].length / 3, PEND = H.pjd0 + (PN - 1) * H.pstep;
  const pl = (b, t) => { const a = H.planets[b], f = (t - H.pjd0) / H.pstep, i = Math.max(0, Math.min(PN - 2, Math.floor(f))), q = Math.max(0, Math.min(1, f - i));
    return [0, 1, 2].map(k => a[3 * i + k] + (a[3 * i + 3 + k] - a[3 * i + k]) * q); };
  const N = m => m.sc.length / 3, mEnd = m => m.jd0 + (N(m) - 1) * m.step;
  const scAt = (m, t) => { const f = (t - m.jd0) / m.step, n = N(m) - 1; if (f < 0 || f > n) return null; const i = Math.min(n - 1, Math.floor(f)), q = f - i, a = m.sc;
    return [0, 1, 2].map(k => a[3 * i + k] + (a[3 * i + 3 + k] - a[3 * i + k]) * q); };
  const idx = (m, t) => Math.max(0, Math.min(N(m) - 1, Math.round((t - m.jd0) / m.step)));
  const CAPTURED = { '-77': 1, '-61': 1 };
  const evName = (m, e) => e.body == 'Engine burn' ? 'Engine burn (deep-space manoeuvre)' : e.body == 'Jupiter' && CAPTURED[m.id] ? 'Jupiter orbit insertion' : e.body + ' flyby';
  const aTxt = v => v == null || v < 0 ? 'hyperbolic' : v.toFixed(2) + ' AU';

  // ---------- 3D objects
  const venS = glow('rgba(255,230,170,1)', .022), satS = glow('rgba(235,205,140,1)', .04);
  const orbV = line(0xffe0a0, .3), orbS = line(0xe0c080, .35);
  const full = M.map(m => { const l = line(parseInt(m.color.slice(1), 16), .5), pts = []; for (let i = 0; i < m.sc.length; i += 3) pts.push([m.sc[i], m.sc[i + 1], m.sc[i + 2]]); l.set(pts); l.pts = pts; return l; });
  const spr = M.map(m => glow(hexRGBA(m.color), .024));
  const trail = line(0xffffff, 1);
  const OWN = [venS, satS, orbV, orbS, trail, ...full, ...spr]; OWN.forEach(o => o.visible = false);
  const MAIN = [orbE, orbJ, arcL, fullL, hohL, phE, phJ, trailJ, mkL, mkA, craft];
  function setOrb(l, b, T, t0) { t0 = Math.max(H.pjd0, Math.min(t0, PEND - T)); const p = []; for (let i = 0; i <= 300; i++) p.push(pl(b, t0 + i / 300 * T)); l.set(p); }

  // ---------- panels
  $('hmb').innerHTML = '<button data-i="-1">All 9</button>' + M.map((m, i) => `<button data-i="${i}" style="border-left:4px solid ${m.color}">${m.name}</button>`).join('');
  $('hmb').querySelectorAll('button').forEach(b => b.onclick = () => { scene = -1; $('cap').style.display = 'none'; enter(+b.dataset.i); });
  function panel() {
    $('hmb').querySelectorAll('button').forEach(b => b.classList.toggle('on', +b.dataset.i == st.cur));
    if (st.cur < 0) {
      $('hmi').innerHTML = `Every line is a real trajectory from JPL Horizons. Each mission's actual launch energy is compared with what a direct Lambert transfer would have needed on the same day.<br><br>` +
        M.map(m => `<span style="color:${m.color}">●</span> ${m.name}: C3 <b>${m.c3.toFixed(0)}</b> vs direct <b>${m.direct_c3.toFixed(0)}</b>`).join('<br>');
      return;
    }
    const m = M[st.cur], ratio = m.direct_c3 / m.c3;
    const verdict = ratio > 1.6 ? `<span class="good">${ratio.toFixed(1)}× less launch energy than any direct shot that day</span>` :
      ratio < 0.85 ? `<span class="warn">paid ${(1 / ratio).toFixed(1)}× more than the cheapest direct shot, to go faster</span>` : `<span class="ok2">matches my Lambert direct transfer</span>`;
    $('hmi').innerHTML = `<div style="color:var(--fg);margin:4px 0 6px">${m.blurb}</div>
<div class="g"><span>Launch</span><b>${m.launch}</b><span>Launch C3 (measured, JPL)</span><b>${m.c3.toFixed(1)} km²/s²</b>
<span>Direct shot that day (Lambert)</span><b>${m.direct_c3.toFixed(1)} km²/s²</b><span>Time to Jupiter</span><b>${(m.tof_jup / 365.25).toFixed(2)} yr</b></div>
<div style="margin:5px 0">${verdict}</div><div id="hev">` +
      m.events.map((e, k) => `<div class="ev" data-k="${k}"><b>${e.date}</b> ${evName(m, e)}${e.alt_km != null ? ' · ' + Math.round(e.alt_km).toLocaleString() + ' km' : ''}<br>
<span>a ${aTxt(e.a0)} → ${aTxt(e.a1)} · e ${e.e0.toFixed(2)} → ${e.e1.toFixed(2)}${Math.abs(e.i1 - e.i0) > 8 ? ' · i ' + e.i0.toFixed(0) + '° → ' + e.i1.toFixed(0) + '°' : ''}</span></div>`).join('') + '</div>';
    $('hev').querySelectorAll('.ev').forEach(d => d.onclick = () => { const e = m.events[+d.dataset.k]; jd = e.jd - 25; sl.value = jd - st.base; });
    st.evShown = -1;
  }

  // ---------- C3 vs time-of-flight chart (lambert curves from historic_missions.py)
  const hc = $('hcv').getContext('2d'), CW = 300, CHh = 250, X = d => 34 + d / 2400 * 256, Y = c => CHh - 26 - Math.min(c, 180) / 180 * (CHh - 44);
  function drawC3() {
    const g = hc; g.fillStyle = '#050714'; g.fillRect(0, 0, CW, CHh); g.font = '10px system-ui'; g.lineWidth = 1;
    for (const c of [50, 100, 150]) { g.strokeStyle = 'rgba(140,160,255,.15)'; g.beginPath(); g.moveTo(34, Y(c)); g.lineTo(CW - 6, Y(c)); g.stroke(); g.fillStyle = '#8d96c4'; g.fillText(c, 8, Y(c) + 3); }
    for (const d of [1, 2, 3, 4, 5, 6]) { const x = X(d * 365.25); g.strokeStyle = 'rgba(140,160,255,.15)'; g.beginPath(); g.moveTo(x, 18); g.lineTo(x, CHh - 26); g.stroke(); g.fillStyle = '#8d96c4'; g.fillText(d + ' yr', x - 8, CHh - 12); }
    g.fillText('launch C3 (km²/s²)', 36, 12); g.fillText('time to Jupiter →', CW - 92, CHh - 2);
    const sel = st.cur >= 0 ? [M[st.cur]] : M, curve = new Map();
    sel.forEach(m => m.direct_curve.forEach(([d, c]) => { if (!(curve.get(d) <= c)) curve.set(d, c); }));
    const pts = [...curve].sort((a, b) => a[0] - b[0]);
    g.fillStyle = 'rgba(255,90,90,.12)'; g.beginPath(); g.moveTo(X(pts[0][0]), Y(0)); pts.forEach(([d, c]) => g.lineTo(X(d), Y(c))); g.lineTo(X(pts[pts.length - 1][0]), Y(0)); g.fill();
    g.strokeStyle = '#ff5a5a'; g.lineWidth = 2; g.beginPath(); pts.forEach(([d, c], k) => k ? g.lineTo(X(d), Y(c)) : g.moveTo(X(d), Y(c))); g.stroke();
    g.fillStyle = 'rgba(255,140,140,.9)'; g.fillText(st.cur >= 0 ? 'direct shot, same launch day (Lambert)' : 'cheapest direct shot (Lambert)', 40, Y(170) + 4);
    g.fillText('below the line: impossible without help', 120, Y(8));
    M.forEach((m, k) => { const on = k == st.cur || st.cur < 0, x = X(m.tof_jup), y = Y(m.c3); g.globalAlpha = on ? 1 : .3;
      g.fillStyle = m.color; g.beginPath(); g.arc(x, y, k == st.cur ? 6 : 4, 0, TAU); g.fill();
      if (k == st.cur) { g.strokeStyle = '#fff'; g.lineWidth = 2; g.beginPath(); g.arc(x, y, 9, 0, TAU); g.stroke(); }
      if (on) { g.fillStyle = m.color; g.fillText(m.name, x + 7, y + (m.name == 'Voyager 1' ? -5 : m.name == 'Voyager 2' ? 9 : 3)); } g.globalAlpha = 1; });
  }
  $('hcv').onclick = ev => { const b = ev.target.getBoundingClientRect(), x = (ev.clientX - b.left) / b.width * CW, y = (ev.clientY - b.top) / b.height * CHh;
    let best = -1, bd = 15; M.forEach((m, k) => { const d = Math.hypot(X(m.tof_jup) - x, Y(m.c3) - y); if (d < bd) { bd = d; best = k; } }); if (best >= 0) enter(best); };

  // ---------- energy / semi-major-axis chart (indep_backup output)
  const chB = document.createElement('canvas'); chB.width = 640; chB.height = 300; let chX = 2300;
  const Emin = -560, Emax = 160, yE = v => 28 + 112 * (Emax - Math.max(Emin, Math.min(Emax, v))) / (Emax - Emin), yA = v => 178 + 112 * (1 - Math.min(10, v == null || v < 0 ? 10 : v) / 10);
  function buildCh() {
    const g = chB.getContext('2d'), all = st.cur < 0, list = all ? M : [M[st.cur]]; g.clearRect(0, 0, 640, 300); g.font = '18px system-ui'; g.lineWidth = 1;
    chX = all ? 2300 : (N(list[0]) - 1) * list[0].step; const X2 = d => d / chX * 640;
    g.fillStyle = '#8d96c4'; g.fillText('Orbital energy around the Sun (km²/s²)', 8, 20); g.strokeStyle = 'rgba(140,160,255,.25)'; g.strokeRect(0, 28, 640, 112);
    g.fillText('Semi-major axis a (AU)  ·  top edge = hyperbolic', 8, 170); g.strokeRect(0, 178, 640, 112);
    g.setLineDash([6, 6]); g.strokeStyle = '#ff6a5a'; g.beginPath(); g.moveTo(0, yE(0)); g.lineTo(640, yE(0)); g.stroke();
    g.strokeStyle = '#fff'; g.beginPath(); g.moveTo(0, yA(3.10)); g.lineTo(640, yA(3.10)); g.stroke(); g.setLineDash([]);
    g.fillStyle = '#ff6a5a'; g.fillText('escape from the Sun ↑', 445, yE(0) - 5); g.fillStyle = '#fff'; g.fillText('my Hohmann a = 3.10', 455, yA(3.10) - 5);
    list.forEach(m => { g.strokeStyle = m.color; g.lineWidth = all ? 2 : 3;
      g.beginPath(); m.en.forEach((v, i) => { const x = X2(i * m.step); i ? g.lineTo(x, yE(v)) : g.moveTo(x, yE(v)); }); g.stroke();
      g.beginPath(); m.a.forEach((v, i) => { const x = X2(i * m.step); i ? g.lineTo(x, yA(v)) : g.moveTo(x, yA(v)); }); g.stroke(); g.lineWidth = 1;
      if (!all) m.events.forEach(e => { const x = X2(e.jd - m.jd0); g.strokeStyle = 'rgba(255,255,255,.35)'; g.beginPath(); g.moveTo(x, 28); g.lineTo(x, 290); g.stroke();
        g.fillStyle = '#fff'; g.font = '16px system-ui'; g.fillText(e.body == 'Engine burn' ? '🔥' : e.body[0], x + 3, 46); g.font = '18px system-ui'; }); });
    g.fillStyle = '#8d96c4'; g.font = '14px system-ui'; g.fillText(all ? 'days since each launch →' : 'mission time →', 500, 298);
  }
  function drawCh(m) {
    const CX = $('chc').getContext('2d'); CX.clearRect(0, 0, 640, 300); CX.drawImage(chB, 0, 0); if (!m) return;
    const d = jd - m.jd0; if (d < 0 || d > chX) return; const x = d / chX * 640, i = idx(m, jd);
    CX.strokeStyle = '#fff'; CX.lineWidth = 2; CX.beginPath(); CX.moveTo(x, 28); CX.lineTo(x, 290); CX.stroke();
    CX.fillStyle = m.color; CX.beginPath(); CX.arc(x, yE(m.en[i]), 7, 0, TAU); CX.fill(); CX.beginPath(); CX.arc(x, yA(m.a[i]), 7, 0, TAU); CX.fill();
  }

  // ---------- enter / exit
  function enter(i, o = {}) {
    if (!st.on) { st.on = true; MAIN.forEach(x => x.visible = false); $('pk').style.display = 'none'; $('ds').style.display = 'none'; $('hc').style.display = 'block'; $('hm').style.display = 'block'; $('bh').classList.add('on'); }
    st.cur = i; const all = i < 0, m = all ? null : M[i];
    st.base = all ? H.pjd0 : m.jd0 - 15; st.span = all ? PEND - H.pjd0 : mEnd(m) - st.base;
    sl.max = st.span; jd = o.t != null ? o.t : st.base; sl.value = jd - st.base;
    const t0 = all ? H.pjd0 : m.jd0; setOrb(orbE, '399', 365.25, t0); setOrb(orbJ, '599', 4332.6, t0); setOrb(orbV, '299', 224.7, t0); setOrb(orbS, '699', 10759, t0);
    [orbE, orbJ, orbV, orbS, venS, satS].forEach(x => x.visible = true);
    full.forEach((l, k) => { l.visible = all || k == i; l.material.opacity = all ? .75 : .28; });
    trail.visible = !all; if (m) trail.material.color.set(m.color);
    st.ev = 0; if (m) while (st.ev < m.events.length && m.events[st.ev].jd < jd) st.ev++;
    distT = o.d || (all ? 26 : m.events.some(e => e.body == 'Saturn') ? 26 : 15); ph = o.ph != null ? o.ph : .75;
    cine = !!o.follow; $('bf').classList.toggle('on', cine);
    $('rate').value = o.rate || (all ? 365 : 30); setPlay(o.play !== false); $('jc').style.display = 'none';
    panel(); buildCh(); drawC3();
  }
  function exit() {
    if (!st.on) return; st.on = false; OWN.forEach(x => x.visible = false); MAIN.forEach(x => x.visible = true);
    sl.max = SPAN; ['pk', 'ds'].forEach(id => $(id).style.display = ''); $('hc').style.display = 'none'; $('hm').style.display = 'none'; $('bh').classList.remove('on');
    $('jc').style.display = 'none'; select(mode, true); setT(JD0);
  }
  $('bh').onclick = () => { scene = -1; $('cap').style.display = 'none'; st.on ? exit() : enter(-1); };

  // ---------- per-frame
  function update(now) {
    const dt = Math.min(.05, (now - last) / 1000); last = now; const all = st.cur < 0, m = all ? null : M[st.cur];
    let k = 1; if (cine && m) for (const e of m.events) { const g = Math.abs(e.jd - jd); if (g < 14) k = Math.min(k, Math.max(.02, g / 14)); }
    if (playing) { jd += dt * (+$('rate').value) * k; if (jd >= st.base + st.span) { jd = st.base + st.span; setPlay(false); } sl.value = jd - st.base; }
    dist += (distT - dist) * .06;
    const pE = pl('399', jd), pJ = pl('599', jd), pV = pl('299', jd), pS = pl('699', jd), cp = m ? scAt(m, jd) : null;
    const W = innerWidth, Hh = innerHeight; if (cv.width !== Math.floor(W * R.getPixelRatio()) || cv.height !== Math.floor(Hh * R.getPixelRatio())) { R.setSize(W, Hh, false); cam.aspect = W / Hh; cam.updateProjectionMatrix(); }
    if (cine && cp) { th += dt * .1; let near = 9; [pV, pE, pJ, pS].forEach(p => near = Math.min(near, nrm(sub(p, cp)))); tgtT.set(cp[0], cp[1], cp[2]); distT = near < .35 ? Math.max(.012, near * 2.5) : 2.4; }
    else if (!cine) tgtT.set(0, 0, 0);
    tgt.lerp(tgtT, .07);
    cam.position.set(tgt.x + dist * Math.sin(ph) * Math.cos(th), tgt.y + dist * Math.sin(ph) * Math.sin(th), tgt.z + dist * Math.cos(ph)); cam.lookAt(tgt);
    const sz = (o, f) => o.scale.setScalar(dist * o.userData.size * (f || 1));
    sz(sun); earthS.position.set(...pE); sz(earthS); jupS.position.set(...pJ); sz(jupS); venS.position.set(...pV); sz(venS); satS.position.set(...pS); sz(satS);
    M.forEach((mm, j) => { const p = (all || j == st.cur) ? scAt(mm, jd) : null; spr[j].visible = !!p; if (p) { spr[j].position.set(...p); sz(spr[j], 1.3); } });
    if (m && cp) trail.set(full[st.cur].pts.slice(0, idx(m, jd) + 1).concat([cp])); else if (m) trail.set([[0, 0, 0], [0, 0, 0]]);
    // HUD
    $('date').textContent = fd(jd);
    if (m) {
      if (cp) { const i = idx(m, jd), a = m.a[i], e = m.e[i];
        $('k1').textContent = m.name + ' · distance from Sun'; $('v1').textContent = nrm(cp).toFixed(2) + ' au';
        $('k2').textContent = 'Heliocentric speed'; $('v2').textContent = m.spd[i].toFixed(1) + ' km/s';
        $('st').textContent = (e < 1 ? 'a = ' + a.toFixed(2) + ' AU · ' : 'hyperbolic (escaping) · ') + 'e = ' + e.toFixed(3) + ' · i = ' + m.inc[i].toFixed(1) + '°';
        const nx = m.events.find(ev => ev.jd > jd); $('tm').textContent = nx ? (nx.body == 'Engine burn' ? 'Engine burn' : nx.body == 'Jupiter' && CAPTURED[m.id] ? 'Jupiter arrival' : nx.body + ' flyby') + ' in ' + Math.ceil(nx.jd - jd) + ' d' : 'Mission data ends here';
      } else { $('k1').textContent = m.name; $('v1').textContent = 'on the pad'; $('k2').textContent = 'Launch'; $('v2').textContent = m.launch; $('st').textContent = ''; $('tm').textContent = 'Launch in ' + Math.ceil(m.jd0 - jd) + ' d'; }
      while (st.ev < m.events.length && m.events[st.ev].jd <= jd) { const e = m.events[st.ev]; if (jd - e.jd < 15) { st.flashT = now + 2400; st.flashTx = e.body == 'Engine burn' ? '🔥 BURN' : e.body.toUpperCase(); } st.ev++; }
      while (st.ev > 0 && m.events[st.ev - 1].jd > jd) st.ev--;
      if (st.evShown !== st.ev && $('hev')) { st.evShown = st.ev; $('hev').querySelectorAll('.ev').forEach((d, j) => d.classList.toggle('past', j < st.ev)); }
    } else {
      const fl = M.filter(mm => jd >= mm.jd0 && jd <= mEnd(mm)).map(mm => mm.name);
      $('k1').textContent = 'In flight'; $('v1').textContent = fl.join(', ') || '-'; $('k2').textContent = 'Earth–Jupiter distance'; $('v2').textContent = nrm(sub(pJ, pE)).toFixed(2) + ' au';
      $('st').textContent = 'Real trajectories · JPL Horizons'; $('tm').textContent = '';
    }
    const fls = now < st.flashT; $('flash').style.opacity = fls ? 1 : 0; if (fls) $('flash').textContent = st.flashTx;
    drawDial(pE, pJ, ang(pE, pJ), D.hohmann.real.phi); drawCh(m);
    { const sep = cp ? nrm(sub(pJ, cp)) : 9, sh = sep < .15; $('jc').style.display = sh ? 'block' : 'none';
      if (sh) { drawJC(cp, pJ, jd, Math.max(sep * AU, 8e5)); $('jct').innerHTML = '<b>Jupiter from ' + m.name + '</b><br>' + Math.round(sep * AU).toLocaleString() + ' km'; } }
    R.render(sc, cam);
  }

  // ---------- story scenes (inserted before "Thank you")
  const byName = n => M.findIndex(m => m.name == n), mm = n => M[byName(n)], ev = (n, b) => mm(n).events.find(e => e.body == b);
  const P10 = mm('Pioneer 10'), V1 = mm('Voyager 1'), GAL = mm('Galileo'), JUNO = mm('Juno'), ULY = mm('Ulysses');
  const jE = ev('Juno', 'Earth'), gE = GAL.events.filter(e => e.body == 'Earth'), vJ = ev('Voyager 1', 'Jupiter'), uJ = ev('Ulysses', 'Jupiter');
  const burns = JUNO.events.filter(e => e.body == 'Engine burn');
  const HS = [
    ['How NASA actually did it', `Nine real missions, with trajectories from JPL Horizons. The chart on the right puts each launch next to the direct Lambert transfer line. Anything below the red line was impossible without help.`,
      () => enter(-1, { rate: 365, d: 28 })],
    ['Pioneer 10, 1972', `First spacecraft to Jupiter, on a direct shot. JPL data gives a launch C3 of <b>${P10.c3.toFixed(1)}</b>. My Lambert solver says the cheapest direct shot that day was <b>${P10.direct_c3.toFixed(1)}</b>. My solver reproduces NASA's 1972 launch.`,
      () => enter(byName('Pioneer 10'), { rate: 90 })],
    ['Voyager 1: thrown out of the Solar System', `Watch the energy line. Before Jupiter it is negative, so the spacecraft is bound to the Sun. After the ${vJ.date} flyby it jumps to <b>+${vJ.en1.toFixed(0)} km²/s²</b> and the orbit becomes hyperbolic (e = ${vJ.e1.toFixed(2)}). Jupiter threw it out of the Solar System.`,
      () => enter(byName('Voyager 1'), { t: vJ.jd - 120, rate: 30 })],
    ['Galileo: Venus, Earth, Earth', `Too heavy for a direct shot: launch C3 was only <b>${GAL.c3.toFixed(0)}</b>, against <b>${GAL.direct_c3.toFixed(0)}</b> for a direct one. Each flyby raises a: ${gE[0].a0.toFixed(2)} → ${gE[0].a1.toFixed(2)} → <b>${gE[1].a1.toFixed(2)} AU</b>. That last value is my Hohmann ellipse.`,
      () => enter(byName('Galileo'), { rate: 90 })],
    ['Juno: launch weak, on purpose', `Launch C3 was <b>${JUNO.c3.toFixed(0)}</b>, a third of a direct shot. Juno's own engine fired in deep space (🔥 ${burns.map(b => b.date).join(', ')}, found by my code from the energy jumps alone) to send it back to Earth.`,
      () => enter(byName('Juno'), { t: burns.length ? burns[0].jd - 60 : JUNO.jd0, rate: 30 })],
    ['Juno: the Earth flyby', `On ${jE.date} Juno passed <b>${Math.round(jE.alt_km)} km</b> above Earth. Its a went from ${jE.a0.toFixed(2)} to <b>${jE.a1.toFixed(2)} AU</b>, which is my Hohmann transfer ellipse (a = 3.10), with Earth doing the burn for free.`,
      () => enter(byName('Juno'), { t: jE.jd - 30, rate: 5, follow: true })],
    ['Ulysses: off the planets\' plane', `Ulysses went to Jupiter to <b>leave</b> the plane of the planets. One flyby tilted its orbit from i = ${uJ.i0.toFixed(0)}° to <b>${uJ.i1.toFixed(0)}°</b>, so it could fly over the Sun's poles. No rocket could do that directly.`,
      () => enter(byName('Ulysses'), { t: uJ.jd - 200, rate: 30, ph: 1.35, d: 12 })],
  ];
  SC.splice(SC.length - 1, 0, ...HS);
  return { get on() { return st.on; }, base: () => st.base, enter, exit, update };
})();
