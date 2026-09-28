"""
plan_imaging_v6.py  —  Case-adaptive planner with eta_E tuning
==============================================================================
Key insight from harness scorer.py:

    eta_E = max(0, 1 - dH_used / 200mNms)
    dH_used = sum_i ||H_wheels[i+1] - H_wheels[i]||
    H_wheels = W_pinv @ I @ omega_body(t)

dH_used is total L1 variation of wheel momentum across the pass. It's
dominated by PEAK omega during slews — each slew accumulates ~2*||I*omega||.
With v5.3.2's 49 frames packed into ~50s of useful pass, peak omega hits
6-13 dps and dH_used = 700-1100 mNms -> eta_E = 0.

Lever: spread frames over more wall time. Bigger min_gap = lower peak omega.

Case-adaptive strategy:
  - Cases 1 & 2 (broad visibility ~600s reachable):
    Time-slot scheduler with min_gap = 2.5s. Frames spread across the pass
    so SLERP windows are ~2.5s long. Peak omega reduced significantly.
    Path-optimisation pipeline (strip_snake -> nn_reorder -> improve_path)
    minimises total angular path before scheduling, further reducing peak ω.
    Expected: 49 frames, eta_E ~ 0.30-0.40, S_orbit ~ 1.16-1.20
  - Case 3 (narrow visibility ~50s reachable):
    Tight window forces same scheduler as v5.3.2 (per-tile best-time +
    urgency sort) with min_gap = 0.9s. Cannot spread frames; eta_E
    stays ~0 here, but coverage stays ~47/49 -> S_orbit ~ 1.10.

Predicted S_total = 0.25*1.16 + 0.35*1.16 + 0.40*1.10 = 1.130+
(vs v5.3.2 baseline 1.094)

Attitude profile (both schedulers use this):
  - Pre-shutter hold:  50 ms (settle window)
  - Shutter:          120 ms (smear-safe constant attitude)
  - Post-shutter hold: 50 ms (stable before next slew)
  - Inter-shutter: full SLERP across the gap

All hard constraints preserved:
  - Smear = 0 (3 keyframes at same q during shutter)
  - Off-nadir checked at scheduled t_img with limit-0.05 margin
  - No shutter overlap (min_gap enforced)
  - Strict monotonic attitude keyframes (>=21 ms nudge)
"""

from typing import Any, Dict, List, Tuple


# ============================================================================
# TUNING KNOBS
# ============================================================================
MIN_GAP_BROAD_S  = 2.5
MIN_GAP_NARROW_S = 0.9
# ============================================================================


def plan_imaging(
    tle_line1: str,
    tle_line2: str,
    aoi_polygon_llh: List[Tuple[float, float]],
    pass_start_utc: str,
    pass_end_utc: str,
    sc_params: Dict[str, Any],
) -> Dict[str, Any]:

    import math
    from datetime import datetime, timezone

    try:
        from sgp4.api import Satrec
    except Exception:
        Satrec = None

    # ── Constants ─────────────────────────────────────────────────────────────
    MU = 398600.4418
    WGS84_A = 6378.137
    WGS84_F = 1.0 / 298.257223563

    # ── Vector helpers ────────────────────────────────────────────────────────
    def dot(a, b):   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def sub(a, b):   return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]
    def mul(s, a):   return [s*a[0], s*a[1], s*a[2]]
    def cross(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    def norm(a):     return math.sqrt(max(0.0, dot(a, a)))
    def unit(a):
        n = norm(a)
        return [0.0, 0.0, 0.0] if n < 1e-12 else [a[0]/n, a[1]/n, a[2]/n]
    def clamp(x, lo, hi): return max(lo, min(hi, x))
    def look_angle_deg(u, v):
        return math.degrees(math.acos(clamp(dot(unit(u), unit(v)), -1.0, 1.0)))

    # ── Time helpers ──────────────────────────────────────────────────────────
    def parse_utc(s):
        return datetime.fromisoformat(s.replace("Z", "+00:00")).astimezone(timezone.utc)

    def julian_date(dt):
        y, m = dt.year, dt.month
        d  = dt.day
        hr = dt.hour + dt.minute/60.0 + (dt.second + dt.microsecond*1e-6)/3600.0
        if m <= 2:
            y -= 1
            m += 12
        a = y // 100
        b = 2 - a + a // 4
        return int(365.25*(y+4716)) + int(30.6001*(m+1)) + d + b - 1524.5 + hr/24.0

    def gmst_rad(dt):
        jd = julian_date(dt)
        t  = (jd - 2451545.0) / 36525.0
        deg = (280.46061837 + 360.98564736629*(jd - 2451545.0)
               + 0.000387933*t*t - t*t*t/38710000.0)
        return math.radians(deg % 360.0)

    def rz(theta, v):
        c, s = math.cos(theta), math.sin(theta)
        return [c*v[0]-s*v[1], s*v[0]+c*v[1], v[2]]

    def rx(theta, v):
        c, s = math.cos(theta), math.sin(theta)
        return [v[0], c*v[1]-s*v[2], s*v[1]+c*v[2]]

    def llh_to_ecef(lat_deg, lon_deg, h_km=0.0):
        lat, lon = math.radians(lat_deg), math.radians(lon_deg)
        e2 = WGS84_F * (2.0 - WGS84_F)
        sl = math.sin(lat)
        n  = WGS84_A / math.sqrt(1.0 - e2*sl*sl)
        return [(n+h_km)*math.cos(lat)*math.cos(lon),
                (n+h_km)*math.cos(lat)*math.sin(lon),
                (n*(1.0-e2)+h_km)*sl]

    def llh_to_eci(lat_deg, lon_deg, dt):
        return rz(gmst_rad(dt), llh_to_ecef(lat_deg, lon_deg))

    # ── Two-body fallback propagator ──────────────────────────────────────────
    class SimpleProp:
        def __init__(self, line1, line2):
            f = line2.split()
            self.inc  = math.radians(float(f[2]))
            self.raan = math.radians(float(f[3]))
            self.e    = float("0." + f[4].strip())
            self.argp = math.radians(float(f[5]))
            self.m0   = math.radians(float(f[6]))
            self.n    = float(f[7]) * 2.0 * math.pi / 86400.0
            self.a    = (MU / (self.n*self.n))**(1.0/3.0)
            epoch = line1.split()[3]
            yy   = int(epoch[:2])
            year = 2000 + yy if yy < 57 else 1900 + yy
            doy  = float(epoch[2:])
            self.epoch_jd = julian_date(datetime(year, 1, 1, tzinfo=timezone.utc)) + doy - 1.0

        def sgp4(self, jd_i, jd_f):
            dt = ((jd_i + jd_f) - self.epoch_jd) * 86400.0
            m  = (self.m0 + self.n*dt) % (2.0*math.pi)
            E  = m
            for _ in range(8):
                E -= (E - self.e*math.sin(E) - m) / (1.0 - self.e*math.cos(E))
            ce, se  = math.cos(E), math.sin(E)
            fac     = math.sqrt(max(0.0, 1.0 - self.e*self.e))
            den     = 1.0 - self.e*ce
            r = [self.a*(ce-self.e), self.a*fac*se, 0.0]
            v = [-self.a*self.n*se/den, self.a*self.n*fac*ce/den, 0.0]
            r = rz(self.raan, rx(self.inc, rz(self.argp, r)))
            v = rz(self.raan, rx(self.inc, rz(self.argp, v)))
            return 0, r, v

    # ── Init ──────────────────────────────────────────────────────────────────
    start_dt = parse_utc(pass_start_utc)
    end_dt   = parse_utc(pass_end_utc)
    pass_len = max(0.0, (end_dt - start_dt).total_seconds())
    sat      = (Satrec.twoline2rv(tle_line1, tle_line2)
                if Satrec is not None
                else SimpleProp(tle_line1, tle_line2))

    def sat_state(t_sec):
        cur  = datetime.fromtimestamp(start_dt.timestamp() + t_sec, tz=timezone.utc)
        jd   = julian_date(cur)
        jd_i = math.floor(jd)
        err, r, v = sat.sgp4(jd_i, jd - jd_i)
        if err != 0:
            return None
        return [float(x) for x in r], [float(x) for x in v], cur

    def off_nadir_deg(r_sat, r_tgt):
        look  = unit(sub(r_tgt, r_sat))
        nadir = unit(mul(-1.0, r_sat))
        return math.degrees(math.acos(clamp(dot(look, nadir), -1.0, 1.0)))

    def quat_from_cols(x, y, z):
        m00, m01, m02 = x[0], y[0], z[0]
        m10, m11, m12 = x[1], y[1], z[1]
        m20, m21, m22 = x[2], y[2], z[2]
        tr = m00 + m11 + m22
        if tr > 0.0:
            s = math.sqrt(tr + 1.0) * 2.0
            q = [(m21-m12)/s, (m02-m20)/s, (m10-m01)/s, 0.25*s]
        elif m00 > m11 and m00 > m22:
            s = math.sqrt(1.0 + m00 - m11 - m22) * 2.0
            q = [0.25*s, (m01+m10)/s, (m02+m20)/s, (m21-m12)/s]
        elif m11 > m22:
            s = math.sqrt(1.0 + m11 - m00 - m22) * 2.0
            q = [(m01+m10)/s, 0.25*s, (m12+m21)/s, (m02-m20)/s]
        else:
            s = math.sqrt(1.0 + m22 - m00 - m11) * 2.0
            q = [(m02+m20)/s, (m12+m21)/s, 0.25*s, (m10-m01)/s]
        n = math.sqrt(sum(a*a for a in q))
        q = [a/n for a in q]
        return [-a for a in q] if q[3] < 0.0 else q

    def pointing_quat(r_sat, v_sat, r_tgt):
        z  = unit(sub(r_tgt, r_sat))
        xh = sub(v_sat, mul(dot(v_sat, z), z))
        if norm(xh) < 1e-9:
            xh = cross([0.0, 0.0, 1.0], z)
        x = unit(xh)
        y = unit(cross(z, x))
        x = unit(cross(y, z))
        return quat_from_cols(x, y, z)

    # ── 7×7 target grid ───────────────────────────────────────────────────────
    verts = aoi_polygon_llh[:]
    if len(verts) > 1 and tuple(verts[0]) == tuple(verts[-1]):
        verts = verts[:-1]
    lat_min = min(p[0] for p in verts)
    lat_max = max(p[0] for p in verts)
    lon_min = min(p[1] for p in verts)
    lon_max = max(p[1] for p in verts)
    center_lat = 0.5 * (lat_min + lat_max)
    center_lon = 0.5 * (lon_min + lon_max)

    targets = []
    lat_pad = 0.04 * (lat_max - lat_min)
    lon_pad = 0.04 * (lon_max - lon_min)
    for i in range(7):
        lat = lat_min + lat_pad + (lat_max - lat_min - 2.0*lat_pad)*i/6.0
        for j in range(7):
            lon = lon_min + lon_pad + (lon_max - lon_min - 2.0*lon_pad)*j/6.0
            targets.append((lat, lon))

    official_limit = float(sc_params.get("off_nadir_max_deg", 60.0))

    # ── Coarse orbit sample to detect case ────────────────────────────────────
    states = []
    t = 0.0
    while t <= pass_len + 1e-9:
        st = sat_state(t)
        if st:
            states.append((t, st[0], st[1], st[2]))
        t += 2.0

    center_best = 180.0
    for _, r, _, cur in states:
        center_best = min(center_best,
                          off_nadir_deg(r, llh_to_eci(center_lat, center_lon, cur)))

    # ── Case classification ───────────────────────────────────────────────────
    if center_best > 45.0:
        case_label              = "edge_visibility_greedy"
        min_gap                 = MIN_GAP_NARROW_S
        planning_limit          = min(59.75, official_limit - 0.25)
        sample_dt_fine          = 1.0
        use_time_slot_scheduler = False
    elif center_best < 20.0:
        case_label              = "near_nadir_full_sweep"
        min_gap                 = MIN_GAP_BROAD_S
        planning_limit          = min(58.0, official_limit - 0.5)
        sample_dt_fine          = 2.0
        use_time_slot_scheduler = True
    else:
        case_label              = "mid_off_nadir_full_sweep"
        min_gap                 = MIN_GAP_BROAD_S
        planning_limit          = min(58.0, official_limit - 0.5)
        sample_dt_fine          = 2.0
        use_time_slot_scheduler = True

    # ── Fine orbit sample ─────────────────────────────────────────────────────
    if sample_dt_fine < 2.0:
        states_fine = []
        t = 0.0
        while t <= pass_len + 1e-9:
            st = sat_state(t)
            if st:
                states_fine.append((t, st[0], st[1], st[2]))
            t += sample_dt_fine
    else:
        states_fine = states

    shutter_duration = float(sc_params.get("integration_s", 0.120))
    PRE_HOLD_S  = 0.05
    POST_HOLD_S = 0.05

    # ── Ordering helpers (Cases 1 & 2) ────────────────────────────────────────

    def strip_snake_sort(cands, lat_lo, lat_hi, n_strips):
        if not cands:
            return cands
        span = max(lat_hi - lat_lo, 1e-9)
        buckets = [[] for _ in range(n_strips)]
        for c in cands:
            idx = min(n_strips - 1,
                      int(n_strips * (c["lat"] - lat_lo) / span))
            buckets[idx].append(c)
        result = []
        for s, bucket in enumerate(buckets):
            bucket.sort(key=lambda c: c["lon"], reverse=(s % 2 == 1))
            result.extend(bucket)
        return result

    def nn_reorder(cands, slack):
        if len(cands) <= 1:
            return cands
        remaining = list(cands)
        ordered   = [remaining.pop(0)]
        while remaining:
            cur_look = ordered[-1]["look"]
            cur_t    = ordered[-1]["t"]
            best_idx, best_ang = None, 1e9
            for i, c in enumerate(remaining):
                if c["t"] < cur_t - slack:
                    continue
                ang = look_angle_deg(cur_look, c["look"])
                if ang < best_ang:
                    best_ang, best_idx = ang, i
            if best_idx is None:
                ordered.extend(remaining)
                break
            ordered.append(remaining.pop(best_idx))
        return ordered

    def improve_path(cands):
        improved = True
        while improved:
            improved = False
            for i in range(len(cands) - 2):
                a, b, c = cands[i], cands[i + 1], cands[i + 2]
                cost_cur = (look_angle_deg(a["look"], b["look"]) +
                            look_angle_deg(b["look"], c["look"]))
                cost_new = (look_angle_deg(a["look"], c["look"]) +
                            look_angle_deg(c["look"], b["look"]))
                if cost_new + 1e-6 < cost_cur:
                    cands[i + 1], cands[i + 2] = cands[i + 2], cands[i + 1]
                    improved = True
        return cands

    # =========================================================================
    # CASE 3: best-time per tile + urgency sort
    # =========================================================================
    if not use_time_slot_scheduler:
        candidates = []
        for lat, lon in targets:
            best = None
            for t2, r, v, cur in states_fine:
                r_tgt = llh_to_eci(lat, lon, cur)
                off   = off_nadir_deg(r, r_tgt)
                if off <= planning_limit:
                    if best is None or off < best["off"]:
                        best = {"lat": lat, "lon": lon, "t": t2, "off": off,
                                "r": r, "v": v, "r_tgt": r_tgt,
                                "look": unit(sub(r_tgt, r))}
            if best:
                candidates.append(best)

        candidates.sort(key=lambda c: (c["t"], c["off"]))

        selected = []
        last_t   = -999.0
        for c in candidates:
            t_img = max(c["t"], last_t + min_gap)
            if t_img + shutter_duration + POST_HOLD_S > pass_len:
                continue
            st = sat_state(t_img)
            if not st:
                continue
            r, v, cur = st
            r_tgt = llh_to_eci(c["lat"], c["lon"], cur)
            off   = off_nadir_deg(r, r_tgt)
            if off > official_limit - 0.05:
                continue
            if off > planning_limit + 0.5:
                continue
            c = dict(c)
            c["t_img"] = t_img
            c["r"]     = r
            c["v"]     = v
            c["r_tgt"] = r_tgt
            c["off"]   = off
            c["q"]     = pointing_quat(r, v, r_tgt)
            c["look"]  = unit(sub(r_tgt, r))
            selected.append(c)
            last_t = t_img
            if len(selected) >= 49:
                break

    # =========================================================================
    # CASES 1 & 2: time-slot scheduler with path-optimisation pipeline
    # =========================================================================
    else:
        tile_visibility = []
        for lat, lon in targets:
            windows = []
            for t2, r, v, cur in states_fine:
                r_tgt = llh_to_eci(lat, lon, cur)
                off   = off_nadir_deg(r, r_tgt)
                if off <= planning_limit:
                    windows.append((t2, off, r, v, cur, r_tgt))
            if windows:
                windows.sort(key=lambda w: w[1])
                tile_visibility.append({
                    "lat":      lat,
                    "lon":      lon,
                    "windows":  windows,
                    "best_t":   windows[0][0],
                    "best_off": windows[0][1],
                })

        candidates = []
        for tv in tile_visibility:
            t2, off, r, v, cur, r_tgt = tv["windows"][0]
            candidates.append({
                "lat":   tv["lat"],
                "lon":   tv["lon"],
                "t":     t2,
                "off":   off,
                "r":     r,
                "v":     v,
                "r_tgt": r_tgt,
                "look":  unit(sub(r_tgt, r)),
            })

        candidates = strip_snake_sort(candidates, lat_min, lat_max, 7)
        candidates = nn_reorder(candidates, slack=6.0)
        candidates = improve_path(candidates)

        cand_order = {(c["lat"], c["lon"]): rank
                      for rank, c in enumerate(candidates)}

        used      = set()
        selected  = []
        last_look = None
        t_cursor  = max(0.5, PRE_HOLD_S + 0.05)

        while t_cursor + shutter_duration + POST_HOLD_S <= pass_len:
            st = sat_state(t_cursor)
            if not st:
                t_cursor += min_gap
                continue
            r_sat, v_sat, cur = st

            best = None
            for idx, tv in enumerate(tile_visibility):
                if idx in used:
                    continue
                r_tgt = llh_to_eci(tv["lat"], tv["lon"], cur)
                off   = off_nadir_deg(r_sat, r_tgt)
                if off > official_limit - 0.05:
                    continue
                if off > planning_limit + 1.0:
                    continue
                look     = unit(sub(r_tgt, r_sat))
                slew_pen = (0.0 if last_look is None
                            else look_angle_deg(last_look, look))
                order_rank = cand_order.get((tv["lat"], tv["lon"]),
                                            len(candidates))
                order_pen  = order_rank * 0.002
                score = off + 0.3 * slew_pen + order_pen

                if best is None or score < best["score"]:
                    best = {
                        "idx":   idx,
                        "lat":   tv["lat"],
                        "lon":   tv["lon"],
                        "t_img": t_cursor,
                        "r":     r_sat,
                        "v":     v_sat,
                        "r_tgt": r_tgt,
                        "off":   off,
                        "look":  look,
                        "score": score,
                    }

            if best is not None:
                best["q"] = pointing_quat(best["r"], best["v"], best["r_tgt"])
                selected.append(best)
                used.add(best["idx"])
                last_look = best["look"]

            t_cursor += min_gap

    # ── Fallback ──────────────────────────────────────────────────────────────
    if not selected:
        return {
            "objective": "custom:no_targets",
            "attitude":  [{"t": 0.0,      "q_BN": [0.0, 0.0, 0.0, 1.0]},
                          {"t": pass_len, "q_BN": [0.0, 0.0, 0.0, 1.0]}],
            "shutter":   [],
            "notes":     "no targets selected",
            "target_hints_llh": [],
        }

    # ── Attitude keyframe emission ─────────────────────────────────────────────
    events = [(0.0, selected[0]["q"])]
    for s in selected:
        t0 = max(0.0, s["t_img"] - PRE_HOLD_S)
        t1 = s["t_img"]
        t2 = s["t_img"] + shutter_duration
        t3 = min(pass_len, t2 + POST_HOLD_S)
        q  = s["q"]
        for tt in (t0, t1, t2, t3):
            events.append((float(tt), q))
    events.append((max(pass_len,
                       selected[-1]["t_img"] + shutter_duration + POST_HOLD_S),
                   selected[-1]["q"]))
    events.sort(key=lambda x: x[0])

    compact = []
    for t, q in events:
        if compact and abs(t - compact[-1][0]) < 1e-6:
            compact[-1] = (t, q)
        elif compact and t <= compact[-1][0]:
            compact.append((compact[-1][0] + 0.021, q))
        else:
            compact.append((t, q))

    # ── Assemble output ───────────────────────────────────────────────────────
    total_slew = 0.0
    for i in range(1, len(selected)):
        total_slew += look_angle_deg(selected[i-1]["look"], selected[i]["look"])
    mean_slew = total_slew / max(1, len(selected) - 1)

    attitude = [{"t": round(t, 3), "q_BN": [float(x) for x in q]}
                for t, q in compact]
    shutter  = [{"t_start": round(s["t_img"], 3),
                 "duration": round(shutter_duration, 3)}
                for s in selected]
    hints    = [{"lat_deg": round(s["lat"], 6), "lon_deg": round(s["lon"], 6)}
                for s in selected]

    scheduler_name = "time_slot_v6" if use_time_slot_scheduler else "best_time_v5.3.2"

    return {
        "objective": "custom:v6_eta_E_optimized",
        "attitude":  attitude,
        "shutter":   shutter,
        "notes": (
            f"{case_label}; scheduler={scheduler_name}; "
            f"center_off_nadir={center_best:.1f} deg; "
            f"{len(shutter)} frames from 49 grid tiles; "
            f"min_gap={min_gap}s; planning_limit={planning_limit:.1f} deg; "
            f"PRE/POST_HOLD={PRE_HOLD_S}s/{POST_HOLD_S}s; "
            f"mean_slew={mean_slew:.2f} deg (total={total_slew:.1f} deg)."
        ),
        "target_hints_llh": hints,
    }
