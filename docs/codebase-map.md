# Codebase Map

Quick reference for implementation agents. Prevents 60k-token exploration sweeps.

## Module Layout

```
src/
├── lib.rs              — pub mod declarations (lines 1–21), re-exports (lines 23–30)
├── types.rs            — all domain types, enums, newtypes (778 lines)
├── constants.rs        — physical constants, epochs, unit conversions (199 lines)
├── corrections.rs      — relativistic corrections: meff (lookup), aberr_light (Lorentz), deflect_light (GR bending)
├── flags.rs            — bitflags! structs: CalcFlags, SiderealBits, etc. (146 lines)
├── error.rs            — Error enum
├── context.rs          — Ephemeris (effective_config — swisseph-rs/112: per-call
│                          ephemeris-source resolution from CalcFlags EPHMASK bits, clamped to loaded
│                          backends with Jpl→Swiss→Moshier fallback, adjusts tidal_acceleration via
│                          user_tidal_acceleration (pre-bake value captured in Ephemeris::new);
│                          Cow::Borrowed for matched combos, Cow::Owned only for mismatches;
│                          calc, calc_ut, calc_inner, calc_speed3, extract_for_body, fixstar2, fixstar2_with_config (swisseph-rs/79, mirrors calc_with_config — threads an explicit &EphemerisConfig into calc_fixstar/calc_fixstar_moshier/_sweph/_jpl so TOPOCTR gets a per-call topographic override instead of always reading self.config; fixstar2 delegates to it with &self.config), fixstar2_ut, fixstar2_mag, calc_fixstar, houses/houses_ex/houses_ex2 — UT-based house wrappers: ARMC+eps+nutation setup, Sun-declination resolution for Sunshine, traditional-sidereal dispatch, RADIANS conversion, gauquelin_sector_geometric — swe_gauquelin_sector imeth 0/1 port: own ARMC/eps/nutation setup using the caller's flags directly (NOT houses_ex2's forced TIDAL_DEFAULT — swe_gauquelin_sector's C deltaT call genuinely carries the caller's ephemeris-source iflag), calc() or fixstar2() per starname, house_pos with hsys=Gauquelin (swisseph-rs/97: starname threaded through); gauquelin_sector — full swe_gauquelin_sector dispatcher (swisseph-rs/89, PNOC 8): imeth 0/1 → gauquelin_sector_geometric, imeth 2–5 → gauquelin_sector_risetrans (rise/set-based: finds bracketing rise+set via Ephemeris::rise_trans with DISC_CENTER/NO_REFRACTION per imeth, §7 bracket+re-search logic, linear interpolation into sectors 1–36; circumpolar bodies propagate Err)), EphemerisConfig, CalcResult; stars: StarCatalog field on Ephemeris
├── math.rs             — pure math functions: normalize, chebyshev, cartpol, cotrans, poly_eval
├── date.rs             — Julian Day ↔ calendar conversion, delta-T, UTC
├── obliquity.rs        — swi_epsiln port: all 11 obliquity models
├── bias.rs             — swi_bias port: GCRS↔J2000 frame rotation; icrs2fk5 (RB matrix), fk4_fk5 (B1950 RA correction)
├── precession.rs       — swi_precess port: 3 algorithm families, 11 models, JPLHOR paths
├── deltat/
│   ├── mod.rs          — calc_deltat dispatcher, 5 historical models, Bessel interpolation, future extrapolation, tidal correction
│   └── data.rs         — static tables: DT (409 entries 1620–2028), DT97 (43), DT2 (27), DTCF16 (54×6 spline)
├── nutation/
│   ├── mod.rs          — router + 5 algorithms: IAU 1980, Herring 1987, IAU 2000A/B, Woolard
│   └── data.rs         — generated nutation term tables (IAU 2000A, 2000B, 1980)
├── sidereal_time.rs    — swe_sidtime0/swe_sidtime port: 4 GMST models, 33-term EoE, long-term model
│                          (sidtime_long_term, pre-1850/post-2050 dates) — its 2 internal deltaT
│                          calls force TIDAL_DEFAULT via a local deltat_config override, matching
│                          C's swe_deltat_ex(tjd, -1, NULL) sentinel (swephlib.c:3291,3301), which
│                          always resolves SE_TIDAL_DEFAULT regardless of the actually-configured
│                          ephemeris backend — the same deltaT/tid_acc-inconsistency pattern as
│                          houses_ex2's; discovered via a gauquelin_sector golden-test mismatch at
│                          1600/1800 AD epochs (swisseph-rs/66)
├── calc.rs             — requested_source (swisseph-rs/112: extracts EPHMASK from CalcFlags with
│                          C precedence MOSEPH > JPLEPH > SWIEPH, returns Option<EphemerisSource>);
│                          calc_planet, calc_sun, calc_moon, calc_mean_node, calc_mean_apogee, calc_ecl_nut, extract_output, extract_ecl_nut, plaus_iflag, speed3_interval, denormalize_positions, calc_speed_3point: light-time, retarded velocity, aberration, deflection pipeline + mean element pipeline + SPEED3 helpers; apparent_planet/apparent_sun/apparent_moon (generic over PositionProvider, used by the sweph/JPL backends); topo_offset helper — computes the observer offset (zero vector when TOPOCTR isn't set) and threads it through both the Moshier pipeline (calc_planet/calc_sun/calc_moon, added to earth_helio directly) and the generic pipeline (added to earth_bary/earth_helio as `xobs`/`xobs_helio`) in place of the plain geocenter wherever it functions as "the observer" (light-time, parallax, aberration, deflection); SEFLG_HELCTR support (added for phenomena, swisseph-rs/83): finish_helctr (shared bias->precess/nutation->app_pos_rest tail) plus a HELCTR early-branch in each of calc_planet (Moshier, niter=0 analytic), apparent_planet (Swiss/JPL, niter=1 + re-eval; the light-time loop's initial dx is heliocentric but its extrapolation base xx0 is BARYCENTRIC, a literal sweph.c:2513-2594 quirk), calc_moon (Moshier) and apparent_moon (Swiss/JPL) — the Moon computes light-time dt ONCE from the heliocentric distance (no loop, sweph.c:4147-4152) and, for Swiss/JPL, re-evaluates moon_geo+earth_bary at t-dt minus the ORIGINAL-epoch sun_bary; plaus_iflag forces NOABERR|NOGDEFL for HELCTR so the heliocentric position is purely geometric+light-time. calc_sun/apparent_sun now handle both Body::Sun and Body::Earth via
│                          an `is_earth: bool` parameter (swisseph-rs/96): heliocentric Sun is the
│                          origin → all-zero (short-circuited in calc_inner); heliocentric/barycentric
│                          Earth threads through the Sun pipeline with conditional frame construction
│                          (HELCTR: earth_bary - sun_bary; BARYCTR: earth_bary directly), a light-time
│                          loop that re-evaluates Earth at t-dt (niter=1, 2 iterations, entering even
│                          for Moshier via swi_moshplan re-eval), and a conditional sign flip (skipped
│                          for HELCTR/BARYCTR). C's Swiss sweplan updates both xearth and xsun to
│                          retarded-time values while JPL's swi_pleph only updates xearth — a C-internal
│                          backend inconsistency (~5e-8 AU); Rust uses retarded sun_bary (matching Swiss
│                          bitwise, JPL within 5e-6°). BARYCTR Earth supported for Swiss/JPL (Moshier
│                          rejects BARYCTR globally). SEFLG_BARYCTR for planets + Moon
│                          (swisseph-rs/129): BARYCTR early-branches in apparent_planet (planet_bary
│                          directly, same 2-pass light-time as HELCTR but no sun_bary subtraction)
│                          and apparent_moon (moon_geo + earth_bary, single dt, no sun_bary); both
│                          share finish_helctr tail. Sun BARYCTR: early-return in apparent_sun
│                          (C's app_pos_etc_sbar, sweph.c:4254 — separate from app_pos_etc_sun
│                          because SEI_SUN==SEI_EARTH in C) returns sun_bary with single analytic
│                          light-time retardation, via finish_helctr. BARYCTR now fully implemented
│                          for all bodies on Swiss/JPL; Moshier rejects globally. Osculating lunar node/apogee (swisseph-rs/84, PNOC 3): plan_for_osc_elem (swi_plan_for_osc_elem port — rotates a raw pre-bias geocentric moon 6-vector into ecliptic-of-date via bias->precess->nutation matrix->ecliptic->ecliptic-nutation; ALWAYS precesses to date + uses obliquity-of-date because the SEFLG_J2000/SIDEREAL skips live inside the dead `#ifdef SID_TNODE_FROM_ECL_T0` — the ref doc's Part C pseudocode is WRONG on this, verified against sweph.c:5758; speed is a PURE rotation, precessed via `precess` not `precess_speed` and nutated with nutv=None); lunar_osc_elem (sweph.c:5168, D.2-D.4: node tangent-line direction + osculating-ellipse apogee computed together, node/apogee speeds via plain central difference over speed_intv — the D.2 quadratic node speed is dead code overwritten by D.3; osc_output_frame does the xreturn[24] assembly incl. the SEFLG_J2000 re-projection elif that precesses the of-date node/apogee to J2000, position via precess/speed via precess_speed, AND returns the x2000 J2000-equatorial vector the SEFLG_SIDEREAL ECL_T0/SSY_PLANE rigorous branches need — built by removing the full nutation matrix [nutate(..., backward=true)] from the equatorial-of-date vector + precessing to J2000, per sweph.c:5527-5540; calc_inner threads this x2000 to apply_sidereal so ECL_T0/SSY_PLANE no longer silently fall back to traditional ayanamsa subtraction, swisseph-rs/84 review fix); raw_osc_moon_moshier/sweph/jpl (per-backend raw geocentric-equ-J2000 pre-bias moon pos+vel, reusing compute()/SwephProvider::moon_geo/JplProvider::moon_geo). CRITICAL: the C MOSEPH branch (sweph.c:5336-5354) has NO light-time correction — only JPL/SWIEPH re-evaluate at t-dt; Moshier node/apogee speed carries a documented stateless artifact (< 4e-6 deg/day, see CLAUDE.md <stateless_tolerance> §3).
│                          Planetocentric positions (swisseph-rs/90, PNOC 9): pctr_light_time
│                          (§3a-§3c light-time iteration with center as observer, returns retarded
│                          epoch + dt + speed correction xxsp), pctr_pipeline (§4-§9: planetocentric
│                          subtraction, deflection with Earth-observer geometry via earth_bary/sun_bary,
│                          aberration with center-body velocity, frame bias at retarded t, precession
│                          at tjd, nutation at §1 priming epoch incl. nutv for speed). Deflection uses
│                          earth_helio = earth_bary - sun_bary and planet_for_defl = xx + earth_helio
│                          (self-consistent with the standard pipeline's convention); Moshier returns
│                          Err (BARYCTR unsupported). Ephemeris::calc_pctr in context.rs orchestrates
│                          (§0-§2 validation + provider-level barycentric fetches via pctr_bary_state,
│                          §3d re-eval at t, delegates to pctr_pipeline for §4-§9, sidereal tail)
│                          Asteroid calc (swisseph-rs/101): normalize_asteroid_aliases
│                          (Pluto-134340 + Ceres..Vesta(1..4) alias, called at top of
│                          calc_inner); AsteroidProvider<P: PositionProvider> wraps an inner
│                          provider (SwephProvider/JplProvider/MoshierEarthProvider), evaluates
│                          asteroid from .se1 file, unconditionally adds sun_bary for helio→bary
│                          (seas files don't set SEI_FLG_HELIO, C checks slot index not flag);
│                          MoshierEarthProvider returns Earth via Moshier + all-zero sun_bary
│                          (fresh-process MOSEPH semantics); calc_asteroid_sweph/jpl/moshier
│                          construct the provider stack and delegate to apparent_planet
├── topocentric.rs      — get_observer: swi_get_observer port (NONUT-forced mean-frame path only, docs/c-ref-topocentric.md §3), geodetic→geocentric flattening + diurnal rotation + precession to J2000, returns observer position+velocity offset (AU/AU-day) from the geocenter
├── moshier/
│   ├── mod.rs          — PlantTbl struct, PLANETS array re-export, element-count tests
│   ├── backend.rs      — compute() public API, compute_pipeline() for calc.rs, embofs_mosh, planet/earth velocity helpers, Body dispatch
│   ├── moon.rs         — moshmoon2() lunar series evaluator: MeanElements, mean_elements(), chewm(), moon1–4, mean_node(), mean_apogee(), correction interpolation
│   ├── moon_tables.rs  — generated const arrays: LR/MB/LRT/BT/LRT2/BT2 + z[25] + MEAN_NODE_CORR[304] + MEAN_APSIS_CORR[304]
│   ├── planets.rs      — moshplan2() series evaluator, sscc() harmonic recurrence, fundamental argument constants
│   └── tables.rs       — generated const arrays: 9 planet tables (do not hand-edit, see scripts/gen_moshier_tables.py)
├── jpl/
│   ├── mod.rs          — JplFile (mmap + JplHeader), JplFile::open, byte_order/header/bytes accessors. Re-exports ByteOrder, JplHeader. J_* body index constants. pub fn jpl_pleph (body assembly entry point).
│   ├── header.rs       — ByteOrder enum + Reader cursor, detect_byte_order (plausibility of ss[2]), parse_header (record-0 offsets), compute_ksize (ipt[] algorithm), validate_file_length, JplHeader struct
│   └── interp.rs       — read_record (mmap→Vec<f64>), interp (JPL forward-recurrence Chebyshev eval + sub-interval selection), state (record selection + body interpolation loop)
├── sweph_file/
│   ├── mod.rs          — SwissEphFile (mmap-based .se1 reader), body_file_id(Body → ipl value), evaluate_body re-export
│   ├── types.rs        — FileHeader, PlanetFileData, FileType, ByteOrder, SEI_*/SE_* body constants, SEI_FLG_* flags
│   ├── parse.rs        — binary format parser: Reader cursor, detect_byte_order, parse_file (header + per-planet metadata)
│   ├── segment.rs      — Chebyshev coefficient unpacking from mmap'd bytes: 6 packing modes (4/3/2/1-byte, nibble, quarter-byte)
│   └── evaluate.rs     — rot_back (orbital-plane→ecliptic/equatorial transform), evaluate_body (public API: file + body_id + jd → [x,y,z,vx,vy,vz])
├── houses.rs           — AscMc, HouseResult (public types); houses_armc driver (swe_houses_armc_ex2 port);
│                          calc_h (CalcH core; systems implemented: A/D/N/V/W equal-family,
│                          O/S/X/M/F quadrant-arithmetic, R/C/T/H/J great-circle/pole-height,
│                          P/K/G Newton-iteration (Placidus/Koch/Gauquelin-36),
│                          U/Y/L/Q closed-form misc (Krusinski cotrans chain, APC apc_sector,
│                          Pullen SD/SR), I/i Sunshine (Treindl + Makransky, sunshine_init,
│                          finite-diff cusp speeds, Porphyry fallback on Makransky ERR) —
│                          B still stubbed Err);
│                          Asc1/Asc2/AscDash core trig, fix_asc_polar, mc_like (shared MC/equasc
│                          projection), polar_shift_subset (shared C/H/J/R polar-circle 180° flip),
│                          NewtonCusp/placidus_newton_cusp (shared P/G Newton-iteration skeleton),
│                          apc_sector (radians-domain helper for Y);
│                          sidereal_houses_trad (traditional sidereal: houses_armc at tropical
│                          armc/eps then subtract ayanamsa from cusps/ascmc except armc; W routed
│                          through Equal + re-fixed to 30° multiples; N re-fixed to (i-1)*30
│                          unconditionally) — ayanamsa is passed in, kept pure;
│                          sidereal_houses_ecl_t0/sidereal_houses_ssypl (geometric sidereal
│                          projections onto the ecliptic-of-t0/solar-system-plane: build a vernal-
│                          point-like moving vector, precess/nutate it to tjde's true equator,
│                          derive an auxiliary obliquity+vernal-point from its orbital-plane normal
│                          via cross_prod/dot_prod_unit, compute houses at the resulting
│                          armcx/epsx, then subtract dvpxe+ayan_t0(+x00 for ssypl) from every
│                          cusp/ascmc except armc) — t0/ayan_t0 passed in raw (unlike
│                          ayanamsa::resolve_t0's callers, C's own sidereal_houses_ecl_t0/ssypl
│                          never apply the t0_is_UT deltaT adjustment); shared helpers
│                          rotate_to_true_equator, sidereal_houses_geom_core, apply_sidereal_shift;
│                          house_pos (swe_house_pos port: planet ecl.[lon,lat] → continuous house
│                          position 1.0..13.0, all 25 HouseSystem variants incl. Alcabitius which
│                          calc_h doesn't implement yet — house_pos's Alcabitius branch is
│                          self-contained; Koch returns Err on circumpolar dfac-out-of-range;
│                          Sunshine requires Some(sundec) in range or Err) — helpers:
│                          armc_to_mc_house_pos (degnorm-before-and-after-+180 variant, distinct
│                          from crate::math::armc_to_mc and mc_like — c-ref-houses.md §12.1),
│                          mc_transform_raw (Morinus's un-normalized tand/cose transform),
│                          bracket_interpolate_12 (shared J/default fallback), koch_house_pos,
│                          topocentric_house_pos, sunshine_apc_house_pos (shared I/i/Y formula)
├── eclipse.rs          — solar eclipse shadow geometry + shared eclipse/occultation helpers
│                          (swisseph-rs/72): EclipseWhere (geopos[0..1]+dcore[0..6]+EclipseFlags,
│                          mirrors C's swe_sol_eclipse_where output minus attr[], which needs the
│                          not-yet-ported eclipse_how); body_radius_au (drad lookup, shared by
│                          eclipse_where/eclipse_how, reuses constants::PLANETARY_DIAMETERS;
│                          named-asteroid diameter not yet threaded through, returns 0.0);
│                          calc_planet_star (shared body/star position dispatch: self.calc vs
│                          self.fixstar2, reused by lunar-eclipse/occultation modules);
│                          eclipse_where (swi_polcart-rebuilt rm/rs from polar lm/ls, NOT the
│                          direct swe_calc cartesian output — literal FP-fidelity hazard per
│                          docs/c-ref-eclipse-solar.md §3.2; two-pass earthobl ellipsoid
│                          refinement, `for niter in 0..2` with `continue` standing in for C's
│                          `goto iter_where`); sol_eclipse_where (public wrapper, pins ipl=Sun —
│                          masks ifl to calc::EPHMASK before calling eclipse_where, matching C's
│                          swe_sol_eclipse_where: this strips NONUT/TOPOCTR/etc. so eclipse_where's
│                          own NONUT branch is unreachable through this entry point — the bug
│                          that cost a debugging session before the C source was read directly).
│                          Ephemeris::sol_eclipse_where in context.rs delegates (same wrapper
│                          pattern as azalt/rise_trans). EclipseHow (attr[0..10]: magnitude,
│                          diameter_ratio, obscuration, core_diameter_km, azimuth, true_altitude,
│                          apparent_altitude, elongation, nasa_magnitude, saros_series,
│                          saros_member, flags); calc_planet_star_topo (calc_planet_star variant
│                          threading a per-call topographic config override through the planet
│                          branch only — stars have no TOPOCTR path yet, matching riseset.rs);
│                          eclipse_how (swisseph-rs/73, local circumstances at an observer:
│                          builds a topo_config clone with `config.topographic` set from
│                          geolon/geolat/geohgt rather than mutating `self` — same pattern as
│                          riseset.rs's `topo_config`; obscuration via the circular-segment
│                          lens-area formula; Saros lookup against the 181-entry
│                          SAROS_DATA_SOLAR table (swecl.c:107-298, ported verbatim as a
│                          `[(i32, f64); 181]` const) with the literal j/j+1 boundary-scan loop,
│                          not a simplified nearest-cycle formula); sol_eclipse_how (public
│                          wrapper: geopos[2] altitude-range validation, CENTRAL/NONCENTRAL merge
│                          from a second eclipse_where call, then a redundant topocentric az/alt
│                          recompute whose apparent-altitude<=0 gate can zero the whole
│                          classification — and attr[0..3]/attr[8..10] with it — while leaving
│                          azimuth/true_altitude/apparent_altitude/elongation populated, matching
│                          swe_sol_eclipse_how's own layered visibility override). Ephemeris::
│                          sol_eclipse_how in context.rs delegates. SolarEclipseGlobal
│                          (swisseph-rs/74, tret[0..7]: time_maximum, time_ra_conjunction,
│                          time_begin/end, time_totality_begin/end, time_centerline_begin/end —
│                          tret[8..9] omitted, unimplemented upstream); contact_dc (shared
│                          contact-time sample formula for n=0 eclipse begin/end, n=1
│                          totality/annularity begin/end, n=2 center-line begin/end — n=0
│                          literally divides by cos_umbra_half_angle not cos_penumbra_half_angle,
│                          a verbatim C quirk per docs/c-ref-eclipse-solar.md §5.5);
│                          sol_eclipse_when_glob (swe_sol_eclipse_when_glob port: `'next_try: loop`
│                          + `continue 'next_try` standing in for C's `goto next_try`; Meeus
│                          lunation stepping with the [21°,159°] F-argument rejection band;
│                          find_maximum-refined minimum-separation search then a 3-pass
│                          fixed-point ET→UT deltaT conversion; eclipse_where+eclipse_how
│                          confirmation; ifltype bit-cascade rejection/retry; find_zero-refined
│                          contact times with a 3-pass Newton polish; annular-total (hybrid)
│                          detection via a core-shadow sign change; secant-refined RA-conjunction
│                          instant (tret[1]); diverges from C only in that a negative-discriminant
│                          find_zero leaves tret[i1]/tret[i2] at 0.0 instead of C's stale-value
│                          carryover — unreachable for well-conditioned real eclipses). Ephemeris::
│                          sol_eclipse_when_glob in context.rs delegates.
│                          meeus_new_moon_estimate (swisseph-rs/75: the Meeus lunation-estimate +
│                          F-argument filter factored out of sol_eclipse_when_glob's inline block
│                          and reused by eclipse_when_loc, per the C ref doc's "factor this into
│                          one shared function" porting note — no behavior change to
│                          sol_eclipse_when_glob). SolarEclipseLocal (tret[0..7]: time_maximum,
│                          time_first_contact..time_fourth_contact (1st/4th = penumbra, 2nd/3rd =
│                          umbra — DIFFERENT tret[] index semantics than SolarEclipseGlobal, see
│                          c-ref-eclipse-solar.md §6.3), time_sunrise/time_sunset, attr:
│                          EclipseHow, flags); eclipse_when_loc (swe_sol_eclipse_when_loc's worker
│                          port, swisseph-rs/75: topocentric main convergence loop with a 2-then-3
│                          dtdiv step-size schedule (distinct dtstart/dt-floor constants from the
│                          glob search — do not unify); topo_angular_separation (shared
│                          overlap-gap sample helper, manually re-normalizes cartesian distance
│                          via sqrt-of-squares rather than reusing a polar distance component,
│                          matching the C source's own dead-polar-fetch pattern); contacts 2/3
│                          (umbra ingress/egress) apply an asymmetric 0.99916 rmoon correction
│                          (flanking samples corrected, the reused center dc[1] is not) and a
│                          SEFLG_SPEED-based secant refinement that linearly extrapolates a
│                          second sample from one calc call's velocity components rather than
│                          calling calc twice; contacts 1/4 (penumbra ingress/egress) mirror that
│                          shape with no 0.99916 correction and an asymmetric fabs(rsplusrm) (only
│                          in the secant refinement, not the initial sample); visibility scan
│                          (descending i=4..=0 so the i=0/max eclipse_how write survives last,
│                          matching C's shared-attr[]-clobbering behavior) computed via raw ifl
│                          (not the topocentric iflag); sunrise/sunset re-anchor via
│                          Ephemeris::rise_trans (Error::CircumpolarBody short-circuits the
│                          function, matching C's retc==-2 early return) called with the
│                          topocentric iflag — a literal C inconsistency vs. every eclipse_how
│                          call in the same function using raw ifl, preserved not "fixed".
│                          sol_eclipse_when_loc (public wrapper: geopos[2] altitude-range
│                          validation matching sol_eclipse_how, merges only NONCENTRAL — not
│                          CENTRAL — from a geocentric eclipse_where call at time_maximum, unlike
│                          sol_eclipse_how's CENTRAL|NONCENTRAL merge). Ephemeris::
│                          sol_eclipse_when_loc in context.rs delegates.
│                          saros_lookup (swisseph-rs/76: shared Saros series/member scan,
│                          factored out of eclipse_how's solar-only branch and reused by
│                          lun_eclipse_how — identical algorithm over two distinct 180/181-entry
│                          tables, SAROS_DATA_SOLAR/SAROS_DATA_LUNAR); LunarEclipseCore/
│                          lun_eclipse_how (swecl.c:3248-3372, selenocentric Earth-shadow-cone
│                          core: change of origin to the selenocentric frame, umbra/penumbra
│                          half-angles, shadow diameters on the fundamental plane with the
│                          `(1+1/50)` atmospheric enlargement and doubled cosf1/cosf2 division —
│                          literal C, not simplified — plus the 0.99405/0.98813 NASA-agreement
│                          factors, phase/umbral+penumbral magnitude, distance from opposition,
│                          Saros lookup; now also carries `dcore[0..4]` (r0/d0/D0/cosf1/cosf2,
│                          fields `r0`/`d0`/`cap_d0`/`cosf1`/`cosf2`) that `lun_eclipse_when`'s
│                          contact-time search needs — exposed on `LunarEclipseCore` since RSE
│                          10/12, swisseph-rs/77); LunarEclipseHow/swe_lun_eclipse_how (public
│                          wrapper, swecl.c:3190-3239: adds the Moon's topocentric az/alt at
│                          `geopos` via calc_ut_with_config+azalt, forces `flags` empty when
│                          apparent altitude ≤0 while leaving magnitudes/Saros populated —
│                          matching sol_eclipse_how's analogous horizon-visibility gate; unlike C,
│                          `geopos` is always required here rather than nullable — callers needing
│                          the geocentric-only core call lun_eclipse_how directly). Ephemeris::
│                          lun_eclipse_how in context.rs delegates.
│                          lun_contact_dc (swisseph-rs/77, shared contact-time sample formula for
│                          n=0 penumbral/n=1 partial(umbra)/n=2 totality begin-end, mirrors solar's
│                          contact_dc); LunarEclipseGlobal (tret[0,2..7]: time_maximum,
│                          time_partial_begin/end, time_totality_begin/end,
│                          time_penumbral_begin/end — tret[1] omitted, unused for lunar);
│                          lun_eclipse_when (swe_lun_eclipse_when port: `'next_try: loop` +
│                          `continue 'next_try`; Meeus lunation stepping with the (21°,159°)
│                          F-argument rejection band (same structure as solar's but a distinct
│                          17-term full-moon periodic series, not factored into a shared helper —
│                          F1's `sin(Om)` is computed with `Om` already in radians, a literal C
│                          quirk preserved exactly); find_maximum-refined minimum-separation search
│                          (dtdiv=4 fixed, unlike solar's 2-then-3 local-search schedule) then a
│                          3-pass fixed-point ET→UT deltaT conversion; confirmation via the
│                          internal `lun_eclipse_how` core directly (equivalent to C's public-
│                          wrapper-with-geopos=NULL call, skipping the topocentric gate);
│                          ifltype bit-cascade rejection/retry; find_zero-refined contact times
│                          (coarse bracket + 3-round 2-point secant refinement, dt/=2 each round —
│                          a different divisor than solar's dt/=3) with the same intentional
│                          divergence as sol_eclipse_when_glob on `find_zero` failure: slots left
│                          at 0.0 instead of C's stale-value carryover, unreachable for a confirmed
│                          eclipse). Ephemeris::lun_eclipse_when in context.rs delegates.
│                          LunarEclipseLocal (tret[0,2..9]: same index semantics as
│                          LunarEclipseGlobal plus tret[8]/tret[9] for moonrise/moonset — NOT
│                          SolarEclipseLocal's differing 1st/2nd/3rd/4th-contact layout);
│                          lun_eclipse_when_loc (swe_lun_eclipse_when_loc port: visibility scan
│                          (descending i=7..=0, skip i==1) via the public swe_lun_eclipse_how;
│                          moonrise/moonset clipping via Ephemeris::rise_trans, both searches
│                          anchored at `tret[6]-0.001` — Error::CircumpolarBody from either call
│                          skips clipping entirely (matches C's `retc>=0` guard checking the LAST-
│                          assigned retc, not each call's own), other Err propagates; the second
│                          moonset-clipping `if` reads `tret[6]`/`tret[7]` which may already have
│                          been mutated by the first moonrise-clipping `if` in the same iteration —
│                          literal C hazard (ref doc §5), preserved via sequential mutation of one
│                          `tret` array rather than "fixed" to read pre-mutation values; on no
│                          visible phase / degenerate rise-set ordering / final-instant no-eclipse,
│                          retries with `tjd_start = tret[0] ± 25` days). Ephemeris::
│                          lun_eclipse_when_loc in context.rs delegates.
│                          normalize_occulted_body (asteroid-134340->Pluto aliasing, shared by
│                          lun_occult_where/lun_occult_when_glob, swisseph-rs/78); lun_occult_where
│                          (swe_lun_occult_where port: thin wrapper over eclipse_where threading
│                          ipl/starname in place of the Sun — no attr[]/eclipse_how call, same
│                          "geometry only" scope as sol_eclipse_where); OccultGlobal (tret[0..7],
│                          same slot layout as SolarEclipseGlobal but tret[1] is the occulted
│                          body's transit instant, not specifically the Sun's); lun_occult_when_glob
│                          (swe_lun_occult_when_glob port: generic `dl/13` Newton-style Moon-body
│                          ecliptic-longitude bracketing in place of solar's Meeus lunation-number
│                          estimate — occultation search must work for any sidereal period incl. a
│                          fixed star's zero proper motion; reuses contact_dc/find_zero/find_maximum
│                          verbatim from the solar port; two C `eclipse_where` calls with identical
│                          arguments collapsed into one reused result, since eclipse_where is pure
│                          and `tjd` is unchanged between them — makes C's dead "extremely small
│                          percentage" fallback branch provably unreachable; `dtb` NOT divided by 3
│                          here unlike solar's contact-time refinement — literal C divergence;
│                          ifltype ANNULAR/HYBRID validity gated on `ipl==Body::Sun` — ported via
│                          `starname` presence rather than replicating C's `ipl<0`->0 sentinel
│                          clamp; no `SE_ECL_ONE_TRY` support, matching sol_eclipse_when_glob's
│                          bool-only `backward`). Ephemeris::lun_occult_where/lun_occult_when_glob
│                          in context.rs delegate. OccultLocal (tret[0..6]: same slot semantics as
│                          SolarEclipseLocal; for a fixed star, contacts 1/4 alias contacts 2/3);
│                          occult_when_loc (swisseph-rs/79, swe_lun_occult_when_loc's worker port:
│                          topocentric `iflag = TOPOCTR | ifl` WITHOUT EQUATORIAL — confirmed
│                          against swecl.c directly since c-ref-occultation.md's phrasing reads
│                          ambiguously here, a genuine literal divergence from solar's
│                          eclipse_when_loc; rough-conjunction dl/13 bracket shared in spirit with
│                          lun_occult_when_glob; contacts 2/3 always attempted, contacts 1/4
│                          branch on `starname` — planet: independent find_zero refine, star:
│                          `tret[1]=tret[2]`/`tret[4]=tret[3]` alias, matching swecl.c:2696-2699;
│                          occultation-only rise/set block (no solar equivalent): occulted body's
│                          own rise/set anchored at `tret[1]-0.1` (NOT solar's `-0.001`) fills
│                          tret[5]/tret[6], plus two independent Sun rise/set pairs at tret[1]/
│                          tret[4] set OCC_BEG_DAYLIGHT/OCC_END_DAYLIGHT — both quirks confirmed
│                          against swecl.c directly, not documented precisely enough in the ref
│                          doc's step 12 prose to port from the doc alone); lun_occult_when_loc
│                          (public wrapper: occultation's own geoalt error string, distinct from
│                          solar's). Ephemeris::lun_occult_when_loc in context.rs delegates.
│                          **Arc complete** (Phase 11: Rise/Set, Eclipses, Occultations, 12/12).
│                          Uncovered a pre-existing gap while implementing this: `eclipse_how`'s
│                          topocentric az/alt for a fixed star silently fell back to geocentric
│                          (fixstar2 had no per-call topographic override) — no golden test had
│                          ever exercised topocentric fixstar before occ_when_loc's Aldebaran
│                          case. Fixed by threading `config: &EphemerisConfig` through
│                          `calc_fixstar`/`calc_fixstar_moshier`/`_sweph`/`_jpl` (new
│                          `fixstar2_with_config`, mirrors `calc_with_config`) and adding the
│                          topocentric observer offset (`calc::topo_offset`, widened to
│                          `pub(crate)`) to `xobs`/`xobs_dt` when TOPOCTR is set — see
│                          docs/c-ref-fixstar.md step 6. `calc_planet_star_topo` (this file) now
│                          calls `fixstar2_with_config` instead of the geocentric-only `fixstar2`.
├── fictitious.rs       — fictitious planets element layer (swisseph-rs/122, 1/2):
│                          FictitiousCatalog (built-in 15-row Neely table + seorbel.txt
│                          parser with check_t_terms polynomial-in-T evaluator);
│                          load_fictitious_catalog (file path with built-in fallback,
│                          same pattern as stars.rs load_catalog); resolve_elements
│                          (T-term eval, degnorm→radians, mano epoch override);
│                          kepler (Kepler equation solver: fixed-point e<0.4,
│                          Newton e≥0.4, 1e-12 convergence); osc_el_plan
│                          (elements→J2000-equatorial-barycentric 6-vector: Gaussian
│                          PQR rotation, obliquity at tequ, precession to J2000,
│                          xearth/xsun barycentric shift; FICT_GEO flag switches
│                          dmot/K/anchor for geocentric bodies). No calc dispatch
│                          yet — wiring into Ephemeris::calc is 2/2's job.
├── ayanamsa.rs         — EMPTY stub
├── azalt.rs            — atmospheric refraction + horizontal coordinates: refrac (swe_refrac,
│                          Meeus true<->apparent, sea-level/no-dip), refrac_extended (swe_refrac_
│                          extended, Sinclair calc_astronomical_refr + calc_dip horizon-dip, 5-
│                          iteration Newton inversion for TrueToApp), azalt/azalt_rev (pure
│                          geometry cores, take precomputed armc+eps_true — Ephemeris::azalt/
│                          azalt_rev in context.rs resolve ARMC/eps_true/deltaT — forced
│                          TIDAL_DEFAULT via azalt_armc_eps, same pattern as houses_ex2 — and
│                          delegate); RefracDir/AzAltDir/HorDir direction enums (AzAltDir/HorDir
│                          kept distinct since C's SE_ECL2HOR==SE_HOR2ECL==0 collide)
├── riseset.rs          — rise/set/meridian-transit full algorithm (swisseph-rs/70): RiseSetResult
│                          (single JD, UT); rise_trans_true_hor (swe_rise_trans_true_hor port) —
│                          15-point culmination pre-pass (find_maximum-refined, 6-iteration
│                          shrinking window) + mesh insertion + 20-iteration bisection
│                          zero-crossing (sign change + RISE/SET direction match), circumpolar →
│                          Error::CircumpolarBody; calc_mer_trans (4 fixed Newton-like
│                          iterations, 361°/day rate constant) for MTRANSIT/ITRANSIT; shared
│                          closures resolve_xc/rdi_of/sample/refine_sample capture the search's
│                          fixed inputs (star position computed once and reused, matching C's
│                          "stars don't move over 28h" optimization — fixed stars have NO
│                          TOPOCTR support yet, untested by golden data which only covers
│                          Sun/Moon); Ephemeris::rise_trans_true_hor in context.rs is the public
│                          entry point and delegates here (same wrapper pattern as azalt/
│                          azalt_rev). Depends on Ephemeris::calc_ut_with_config/calc_with_config
│                          (context.rs) — calc/calc_ut refactored to thread an explicit
│                          `&EphemerisConfig` through calc_inner/calc_body_*/calc_speed3 so a
│                          caller-supplied `geopos` can override TOPOCTR's observer position
│                          without needing to match Ephemeris's own configured topographic
│                          position (mirrors C's per-call swe_set_topo, but stateless);
│                          azalt_armc_eps widened to pub(crate) for calc_mer_trans's ARMC.
│                          rise_set_fast (swisseph-rs/71, swe_rise_trans_fast port): semi-diurnal-
│                          arc estimate (armc/decl snapshot, sda clamped to 10°/180° for
│                          never-rises/never-sets rather than signaling circumpolar — this path
│                          NEVER returns Error::CircumpolarBody) + 2 Newton iterations (4 for
│                          Moon) via finite-difference azalt slope (0.001-day step), at most one
│                          tjd_ut+0.5 retry if the estimate lands before the input time;
│                          get_sun_rad_plus_refr computes the target altitude offset (disc radius
│                          reusing disc_diameter_m/disc_radius_deg + once-computed horizon
│                          refraction). Ephemeris::rise_trans (context.rs, swe_rise_trans
│                          dispatcher port) picks fast vs rise_trans_true_hor(horhgt=0.0) by the
│                          §4 eligibility gate (not fixstar, RISE/SET only, !FORCE_SLOW, no
│                          twilight, body in Sun..=TrueNode, |geolat|≤60 or Sun≤65).
├── heliacal.rs         — **Arc complete** (Phase 13: Heliacal visibility, swisseph-rs/104–111).
│                          Port of swehel.c: atmospheric extinction/optics model, visibility-limit
│                          magnitude, heliacal phenomena, and event search (both vis_lim and arc_vis
│                          strategies). Public API: Ephemeris::heliacal_ut (swe_heliacal_ut port),
│                          Ephemeris::heliacal_pheno_ut, Ephemeris::vis_limit_mag, heliacal_angle.
│                          **Types**: HeliacalEventType (MorningFirst..AcronymchalSetting),
│                          HeliacalEvent (start_visible/optimum_visibility/end_visible),
│                          HeliacalPheno (28 fields), VisLimitResult (dret[0..7] + vision flags),
│                          HeliacalAngleResult (optimal_altitude/arcus_visionis/sun_altitude_diff).
│                          **Layer 1 — atmosphere/optics** (1/8): default_heliacal_parameters,
│                          extinction (kw/koz/kr/ka/kt/deltam), airmass, optics (cva/pupil_dia/
│                          optic_factor), refraction (topo_alt_from_app_alt/app_alt_from_topo_alt).
│                          **Layer 2 — object location** (3/8): object_loc (7 angle types),
│                          azalt_cart, magnitude, calc_rise_and_set/my_rise_trans/rise_set.
│                          **Layer 3 — vis_limit_mag** (4/8): vis_lim_magn + vis_limit_mag.
│                          **Layer 4 — heliacal_pheno_ut** (6/8): Moon crescent geometry
│                          (width_moon/length_moon/q_yallop/yallop_grade), deter_tav,
│                          heliacal_pheno_ut (visibility-window search + DeterTAV sampling).
│                          **Layer 5 — event search infrastructure** (7/8): get_synodic_period,
│                          TCON table, find_conjunct_sun (Pluto → Err, not C's OOB read),
│                          get_asc_obl/get_asc_obl_diff/get_asc_obl_with_sun (oblique-ascension
│                          bracket+bisection), get_heliacal_day (adaptive day/minute stepping),
│                          get_acronychal_day (photopic-forced convergence),
│                          time_optimum_visibility/time_limit_invisible/get_heliacal_details.
│                          **Layer 6 — event drivers** (8/8, swisseph-rs/111):
│                          vis_lim path: heliacal_ut_vis_lim (conjunction/oblique-ascension seed →
│                          get_heliacal_day/get_acronychal_day → get_heliacal_details),
│                          moon_event_vis_lim (Moon-specific: find_conjunct_sun → get_heliacal_day
│                          → optimum/boundary → sunset/sunrise clamp → TypeEvent==4 reorder);
│                          arc_vis path: heliacal_ut_arc_vis (adaptive day-step halving search with
│                          HeliacalAngle self-adjusting sunsangle feedback, AVKIND_VR per-minute
│                          TAV-minimization via x2min parabola vertex, AVKIND_PTO symmetric-crossing
│                          averaging, AVKIND_MIN7/MIN9 fixed-depth overrides; AltM=-1/AziM=0 always
│                          — Moon interference never factored, matching C's dead #if 0 block;
│                          topo_config threaded via calc_with_config/fixstar2_with_config since C
│                          relies on global swe_set_topo), moon_event_arc_vis (new-moon anchor via
│                          pheno_ut phase-angle walk, per-minute DeterTAV minimization, AVKIND_VR
│                          only for Moon); top-level: heliacal_ut (§7 swe_heliacal_ut port:
│                          Sun/Moon rejection, TypeEvent validation, acronychal 5/6→3/4 remapping
│                          for arc_vis, synodic-period retry loop with MAX_COUNT_SYNPER=5 or
│                          LONG_SEARCH=1M, SEARCH_1_PERIOD post-hoc rejection).
│                          **Deviations from C**: find_conjunct_sun returns Err for Pluto (C has
│                          latent TCON OOB read); arc_vis dret[1..2] explicitly zeroed (C leaves
│                          uninitialized); arc_vis sanity bound uses tjd_start not stale loop var;
│                          moon_event_arc_vis strips TOPOCTR from pheno_ut flags (geocentric
│                          quantity, C uses global swe_set_topo).
│                          **Quirks preserved**: all C algorithmic quirks ported literally for
│                          golden-test parity (see docs/c-ref-heliacal-search.md §8)
├── phenomena.rs        — swe_pheno/swe_pheno_ut port (swisseph-rs/83): Phenomena output struct
│                          (phase_angle, phase, elongation, apparent_diameter, apparent_magnitude,
│                          horizontal_parallax = attr[0..5]); MAG_ELEM[21][4] table + EULER/
│                          EULER_SATURN literals (kept distinct from f64::consts::E for FP
│                          fidelity); pheno core (§1 body remap incl. asteroid-134340->Pluto and
│                          Ceres..Vesta offset; two masked flag copies iflag/iflagp with iflagp
│                          forcing HELCTR; §3 light-time-lagged heliocentric via calc(HELCTR);
│                          magnitude branch cascade 5a-5j — Bowell §5k asteroid branch present but
│                          returns Err(EphemerisNotAvailable) for numbered asteroids, same
│                          swed.ast_H/G gap as eclipse::body_radius_au); pheno_ut deltaT re-call
│                          wrapper. Returns (Phenomena, CalcFlags) — the second is C's "flags
│                          actually used". Everything goes through Ephemeris::calc (constraint
│                          app-uses-calc-not-backends:phenomena). Ephemeris::pheno/pheno_ut
│                          delegates in context.rs; Phenomena re-exported in lib.rs.
├── nodaps.rs           — swe_nod_aps / swe_nod_aps_ut, mean (swisseph-rs/85, PNOC 4) + osculating
│                          (swisseph-rs/86, PNOC 5) branches: NodApsMethod bitflags
│                          (MEAN/OSCU/OSCU_BAR/FOPOINT), NodesApsides output (asc/desc/peri/aphe
│                          [f64;6]), the 5 VSOP mean-equinox-of-date element tables
│                          EL_NODE/PERI/INCL/ECCE/SEMA[8][4]; nod_aps (A.1 remap+reject, A.2 setup,
│                          dispatches to mean_branch or osculating_branch by A.3's eligibility gate);
│                          mean_branch (A.3: Moon via calc::mean_lunar_elements else 4-term
│                          polynomials + cotrans orbital->ecliptic + eccentric-anomaly node distance);
│                          osculating_branch (A.4: instantaneous two-body/angular-momentum ellipse —
│                          A.4.1 reference distance/Gmsm/dt/ellipse_is_bary via
│                          Ephemeris::nodaps_osc_body_j2000, A.4.2 up to 3 samples rotated via
│                          calc::plan_for_osc_elem, A.4.3-A.4.4 per-sample node-tangent-line +
│                          angular-momentum ellipse elements [uu/sema/ecce/ny] producing all FOUR
│                          points [xq perihelion, xa aphelion-or-2nd-focal-point, xn/xs
│                          ellipse-corrected asc/desc node] — the retrograde cosincl flip IS present
│                          here (absent in lunar_osc_elem's D.3, Moon is never retrograde), A.4.5
│                          central-difference assembly); transform_nodaps_output (A.5 shared
│                          pipeline, reused unchanged by both branches with is_true_nodaps=true for
│                          osculating: ecl->equ, precess to J2000, barycenter/observer,
│                          deflect_light+aberr_light, precess back, app_pos_rest tail,
│                          apply_sidereal, extract_output). Goes through calc/context, not backends
│                          (app-uses-calc-not-backends:nodaps): Moon elements via
│                          calc::mean_lunar_elements re-export; observer vectors via
│                          Ephemeris::nodaps_observer (all 3 backends: Moshier xear=heliocentric
│                          Earth/sun_bary=0, Swiss/JPL xear=real barycentric Earth via
│                          SwephProvider/JplProvider.positions(Sun,...) dummy-body call — earth_bary/
│                          sun_bary/earth_helio are always populated regardless of the queried body);
│                          Ephemeris::nodaps_osc_body_j2000 (A.4.1/A.4.2 body-position helper: TRUEPOS
│                          J2000-equatorial pos+speed in the requested helio/bary frame across all 3
│                          backends, Moon always geocentric via the existing raw_moon_at, Earth gets
│                          the EMB correction added back via the same raw_moon_at — Moshier rejects
│                          want_bary=true with Error::UnsupportedFlags(BARYCTR), matching
│                          calc_body_moshier's BARYCTR guard: Moshier has no real SSB).
│                          CONSTRAINT-BOUND PRECISION: mean branch's apparent DESCENDING NODE is
│                          ill-conditioned (docs/swisseph-c-potential-bugs.md §7, ~3.5e-4°pos/
│                          1.2e-2°/day speed, tolerances in tests/golden/nodaps.rs::tolerance);
│                          osculating branch's node-direction `fac=z/ż` division similarly amplifies
│                          backend FP noise into BOTH asc/desc nodes (§8, verified via
│                          identical-input formula check — feeding C's own dumped xpos[1] reproduces
│                          C's uu/cosnode/sinnode to ~12 digits, so the divergence is 100% upstream
│                          backend noise, not a formula bug; worst observed ~1.2e-3°pos/1.9e-2°/day
│                          for Jupiter OSCU_BAR descending node), tiered tolerances in
│                          tests/golden/nodaps.rs::osc_tolerance (peri/aphe tightest, asc looser, desc
│                          loosest). NodApsMethod/NodesApsides re-exported in lib.rs.
│                          PLMASS/IPL_TO_ELEM/OSCU_BAR_DISTANCE_THRESHOLD_AU/NODE_CALC_INTV in
│                          constants.rs (orbit.rs PNOC 6 needs them too). swi_mean_lunar_elements
│                          ported in moshier/moon.rs (mean_lunar_elements; note mean_elements_t2
│                          stale-T2 quirk). SwephProvider/JplProvider (calc.rs) widened to pub(crate)
│                          for context.rs's nodaps_observer/nodaps_osc_body_j2000 reuse.
│                          REVIEW FIX (swisseph-rs/85, /86 codex review): `transform_nodaps_output`
│                          originally ignored SEFLG_HELCTR/BARYCTR entirely — `nodaps_observer`'s
│                          `xobs` always subtracted Earth regardless of flags, so HELCTR/BARYCTR
│                          requests silently returned geocentric output (order-AU error, untested by
│                          any golden case). Fixed per swecl.c:5401-5436's A.5.1 observer-selection
│                          logic (`select_xobs` in nodaps.rs: HELCTR real-ephemerides ->
│                          sun_bary, BARYCTR/HELCTR-on-Moshier -> unchanged/topo-only, else -> +xear).
│                          `ObsFrame.xobs` renamed to `topo` (raw topocentric offset only) so
│                          `select_xobs` can choose per-flag rather than baking Earth in upfront.
│                          Separately, `swi_deflect_light` (sweph.c:3743) ALWAYS reads the true
│                          Earth/Sun-bary globals for its geometry, ignoring whatever the caller's
│                          observer-frame variable holds — only `swi_aberr_light` honors the
│                          HELCTR/BARYCTR-selected xobs. Missing this split initially blew up the
│                          deflection correction under HELCTR (feeding it `sun_bary`, a
│                          near-zero-magnitude vector relative to real Earth-Sun distance, as its
│                          "observer"). Fixed by computing a separate always-true
│                          `earth_helio_true = xear - sun_bary` for deflect_light specifically,
│                          leaving aberr_light/position-shift on the flag-aware `xobs`. Also added:
│                          `Ephemeris::nod_aps`/`nod_aps_ut` now reject `SEFLG_TOPOCTR` without
│                          `config.topographic` (previously silently returned geocentric, matching
│                          `calc_with_config`'s existing guard). Golden coverage added:
│                          "helctr_bary_mean"/"helctr_bary_osc" keys in nodaps.json (24 cases,
│                          Mercury/Jupiter/Pluto × SWIEPH_HELCTR/SWIEPH_BARYCTR/MOSEPH_HELCTR); peri/
│                          aphe longitude gets a 5e-6° floor there (`helctr_bary_tolerance`) for the
│                          2nd-order retarded-Sun-position term this port omits in deflection's
│                          "planethel" vector (order-of-magnitude verified via dt·xsun_speed, matches
│                          calc.rs's own established omission elsewhere). Plus a direct unit test for
│                          the TOPOCTR-without-config rejection.
├── orbit.rs            — swe_get_orbital_elements / swe_orbit_max_min_true_distance (swisseph-rs/87,
│                          PNOC 6): OrbitalElements struct (17 named dret slots + as_array());
│                          get_orbital_elements (swecl.c:5783-5971 — J2000/TRUEPOS/NONUT/XYZ state
│                          vector → Kepler elements: node via r_z/v_z projection, inclination via
│                          r×v, vis-viva sema, anomaly chain, node/apsis refinement, periods);
│                          get_gmsm (swecl.c:5687-5742 — central-body GM, ports the PLMASS/IPL_TO_ELEM
│                          Pluto quirk + ORBEL_AA Mercury double-count literally); osc_get_orbit_
│                          constants/osc_get_ecl_pos/osc_iterate_max_dist/osc_iterate_min_dist
│                          (ellipse sampling + coordinate-descent search — the iterate functions
│                          ALWAYS restart ean=0 and leave `xa` at the overshoot position, both
│                          verified against the C source, not the ref doc); orbit_max_min_true_
│                          distance (swecl.c:6170-6287 — heliocentric single-ellipse branch, else
│                          geocentric two-ellipse grid scan with the literal asymmetric loop-bound
│                          quirk (outer j*2°→362°, inner i*1°→181°) + 300-iter refinement).
│                          SEFLG_ORBEL_AA is bit-aliased onto CalcFlags::TOPOCTR (mass-summation
│                          method, NOT topocentric — the bit never reaches eph.calc). Goes through
│                          Ephemeris::calc/context (app-uses-calc-not-backends:orbit), NOT backends.
│                          HELCTR Earth now uses calc(Earth, HELCTR) directly (swisseph-rs/96
│                          removed the -(geocentric Sun) workaround). BARYCTR for planets now implemented
│                          (swisseph-rs/129); the r>6-AU barycentric branch in
│                          orbit_max_min_true_distance is now reachable on Swiss/JPL
│                          (Moshier still errors — C-Moshier errors there too). Ephemeris::
│                          get_orbital_elements/orbit_max_min_true_distance in context.rs delegate.
│                          OrbitalElements re-exported in lib.rs. PLMASS/IPL_TO_ELEM/
│                          OSCU_BAR_DISTANCE_THRESHOLD_AU in constants.rs (shared with nodaps.rs).
├── crossings.rs        — swe_solcross / mooncross / mooncross_node / helio_cross + UT variants
│                          (swisseph-rs/88, PNOC 7): MoonCrossing struct (jd + longitude + latitude);
│                          solcross/solcross_ut (Sun longitude crossing via Newton iteration,
│                          mean-speed 360/365.24 initial estimate); mooncross/mooncross_ut (Moon
│                          longitude crossing, mean-speed 360/27.32); mooncross_node/mooncross_node_ut
│                          (Moon latitude zero-crossing: day-stepping bracket with fixed-reference
│                          xlat comparison, then Newton refinement jd -= lat/lat_speed);
│                          helio_cross/helio_cross_ut (heliocentric longitude crossing for any planet
│                          including Earth (swisseph-rs/96) except Sun/Moon/nodes/apsides, Chiron uses hardcoded 0.01971 mean-speed
│                          for initial estimate only, dir≥0 forward / dir<0 backward). All go through
│                          Ephemeris::calc/calc_ut (app-uses-calc-not-backends:crossings), NOT backends.
│                          Convergence threshold CROSS_PRECISION = 1 milliarcsecond (1/3600000°).
│                          MoonCrossing re-exported in lib.rs. Ephemeris delegates in context.rs.
└── stars.rs            — StarCatalog, Star, load_catalog, builtin_star (8 ayanamsa ref stars), search, parse

tests/
├── golden/
│   ├── main.rs         — test harness: golden_data_path(), assert_f64_exact(), assert_f64_eps()
│   ├── asteroid.rs    — golden tests for asteroid calc (swisseph-rs/101: 300 cases, 6 bodies
│                         {Chiron..Vesta} × 5 epochs × 10 flag combos {SWIEPH×9, JPLEPH×1},
│                         positions 1e-9 / speeds 1e-7 for SWIEPH; JPLEPH widened to 2e-6 pos /
│                         1e-5 speed (JPL vs SWIEPH Earth/Sun source diff); TopoPosition
│                         configured; retflag checked for SWIEPH only (C returns SWIEPH for
│                         JPLEPH+asteroid since the asteroid file is .se1);
│                         swisseph-rs/102: golden_asteroid_numbered — 72 cases, 4 numbered
│                         asteroids {433 Eros, 7066 Nessus, 136199 Eris (>99999 s%06d naming),
│                         2060 Chiron-as-numbered (SEI_FILE_ANY_AST path)} × 3 epochs × 6 flags
│                         {SWIEPH×5, JPLEPH×1}, same tolerances as main battery;
│                         golden_asteroid_moseph — 27 cases, {Ceres, Vesta, Eros(433)} × 3
│                         epochs × 3 MOSEPH flags, process-isolated generator (separate binary,
│                         no SWIEPH/JPLEPH calls to keep sun_bary zero), positions 5e-7 / speeds
│                         1e-7 (stateless architecture tolerance ~0.4 mas); error/alias tests:
│                         Chiron/Pholus beyond-limits guard (Swiss + Moshier), asteroid not in
│                         asteroid_numbers, asteroid outside file range, alias identities
│                         134340↔Pluto and 1↔Ceres (bitwise-exact))
│   ├── azalt.rs        — golden tests for refraction/horizontal coords (swisseph-rs/69: refrac
│                          28 cases (7 inalt × 2 atpress × 2 dir, exact-or-1e-9 fallback);
│                          refrac_ext 56 cases (× 2 geoalt, out + dret[0..4], exact-or-1e-9);
│                          azalt/azalt_rev 8 cases each (2 tjd_ut × 2 geopos × 2 dir, via
│                          Ephemeris::azalt/azalt_rev, eps 1e-7 — compounds sidtime/obliquity)
│   ├── calc.rs        — golden tests for calc pipeline (1176 cases: 14 bodies × 7 epochs × 12 flag combos incl. SPEED3, no_speed)
│   ├── calc_helctr.rs — dedicated SEFLG_HELCTR/BARYCTR golden tests (swisseph-rs/94, Earth added
│                          swisseph-rs/96, full BARYCTR swisseph-rs/129, 1760 cases:
│                          3 backends {moshier,sweph,jpl} × Sun..Pluto+Moon+Earth (11) × 5 epochs
│                          × 8 HELCTR flag combos + 4 BARYCTR combos for all 11 bodies;
│                          BARYCTR cases for Swiss/JPL only (Moshier rejects). Sun BARYCTR
│                          returns sun_bary (C's app_pos_etc_sbar), NOT earth_bary.
│                          Heliocentric Sun = all-zero (origin), short-circuited in calc_inner.
│                          Positions eps 1e-9, speed eps 1e-7; JPL Earth HELCTR widened to 5e-6/1e-5
│                          for a C-internal backend inconsistency in the light-time loop's sun_bary
│                          handling (Swiss sweplan updates both xearth+xsun to retarded time; JPL
│                          swi_pleph only updates xearth, leaving xsun at original epoch — Rust uses
│                          retarded sun_bary matching Swiss bitwise). Epochs avoid sepl_18 file
│                          boundary. JPL rows skip if ephe/de441.eph absent)
│   ├── calc_topo.rs   — golden tests for SEFLG_TOPOCTR (170 cases across 3 sub-matrices, swisseph-rs/80: moshier — 90 cases, 3 observers × 5 bodies × 3 epochs incl. a SPEED3 file-boundary epoch × 2 flag shapes {speed, speed_noaberr}; sweph — 40 cases, 2 observers × 5 bodies × 2 epochs (incl. the sepl_18 SPEED3 file-boundary epoch, widened tolerance there per the documented C-state artifact) × 2 flag shapes; jpl — 40 cases, same shape as sweph; positions eps 1e-9/speeds eps 1e-7 except the sweph file-boundary widening and an OPEN-BUG widening for jpl epochs != J2000 (swisseph-rs/81 — JPL TOPOCTR diverges from C away from J2000, root cause unconfirmed) — TOPOCTR+SPEED+!NOABERR forces SPEED3 (calc.rs plaus_iflag) for the "speed" shape only; "speed_noaberr" exercises the non-SPEED3 analytic-speed path)
│   ├── fictitious_elements.rs — golden tests for fictitious-planet element layer
│                          (swisseph-rs/122): 114 cases (19 bodies ipl 40–58 × 6 epochs),
│                          element resolution bitwise-exact (assert_f64_exact for all 8
│                          fields: tjd0, tequ, mano, sema, ecce, parg, node, incl),
│                          osc_el_plan output eps 1e-9 pos / 1e-7 vel; C harness includes
│                          swemplan.c directly for read_elements_file/swi_osc_el_plan,
│                          emits xearth/xsun so Rust test feeds identical Earth/Sun state;
│                          bodies 55–58 without-file error case (unit test with builtin only)
│   ├── corrections.rs — golden tests for corrections (30 meff + 40 aberr + 15 pipeline)
│   ├── math.rs         — golden tests for math module
│   ├── date.rs         — golden tests for date module
│   ├── eclipse.rs     — golden tests for eclipse (sol_where: 4 cases — 1999/2021/2024 known
│                          central solar eclipses at their actual maximum-eclipse UT instants
│                          (CENTRAL|TOTAL, CENTRAL|ANNULAR ×2, one with NONUT set to confirm it's
│                          masked away same as C) + a plain-noon no-eclipse epoch; asserts
│                          central_longitude/central_latitude/core_diameter_km eps 1e-7 + exact
│                          retval flags bitmask; sol_how: 8 cases — the same 4 sol_where epochs
│                          (incl. the no-eclipse epoch, exercising the horizon-visibility
│                          clearing path) × 2 observers (8.55,47.37,500 near-central;
│                          -100.0,40.0,0 off-track), asserts all 11 attr[] fields
│                          (magnitude/diameter_ratio/obscuration/core_diameter_km/azimuth/
│                          true_altitude/apparent_altitude/elongation/nasa_magnitude/
│                          saros_series/saros_member) eps 1e-7 + exact retval flags bitmask;
│                          sol_when_glob: 4 cases — 2 tjd_start (2000/2020) × 2 backward,
│                          ifltype=0 (all types), asserts tret[0..7] (time_maximum,
│                          time_ra_conjunction, time_begin/end, time_totality_begin/end,
│                          time_centerline_begin/end) eps 1e-5 day + exact retval flags bitmask;
│                          sol_when_loc: 8 cases — 2 geopos (8.55,47.37,500 near-central for the
│                          sol_where set; -71.0,-33.0,500 Chile, near the 2019/2020 tracks) × 2
│                          tjd_start (2000/2019) × 2 backward, asserts tret[0..7] (time_maximum,
│                          time_first_contact..time_fourth_contact, time_sunrise, time_sunset —
│                          DIFFERENT tret[] index semantics than sol_when_glob, see
│                          c-ref-eclipse-solar.md §6.3) + all 11 attr[] fields eps 1e-5 day/eps +
│                          exact retval flags bitmask; passed on the first implementation attempt,
│                          no escape-hatch escalation needed; lun_how: 3 cases — 2001-01-09 (total,
│                          visible from Zurich), 2021-05-26 (total geocentrically but Moon below
│                          Zurich's horizon, retval==0 while attr[] stays populated — exercises the
│                          horizon-visibility clearing path), 2024-09-18 (small partial) at their
│                          actual maximum-eclipse UT instants × the same near-central Zurich
│                          observer, asserts attr[0,1,4..10] eps 1e-7 + exact retval flags bitmask;
│                          passed on the first implementation attempt; lun_when: 4 cases — 2
│                          tjd_start (2000/2019) × 2 backward, ifltype=0, asserts tret[0,2..7] eps
│                          1e-5 day + exact retval flags bitmask, incl. a partial-only case
│                          (tret[4..5]==0) and a penumbral-only case (tret[2..5]==0); lun_when_loc:
│                          4 cases — 1 geopos (near-central Zurich) × the same 2 tjd_start × 2
│                          backward, asserts tret[0,2..9] eps 1e-5 day + all attr[] fields except
│                          attr[8] (duplicate) + exact retval flags bitmask; both passed on the
│                          first implementation attempt, no escape-hatch escalation needed;
│                          occ_where: 3 cases — Venus/Mars (ipl) + Aldebaran (starname) all at
│                          tjd_ut=2458800.5, via a `make_eph()` with `ephe_path` set (unlike this
│                          file's other tests) so Aldebaran resolves through the fixstar catalog;
│                          asserts central_longitude/central_latitude + all 6 dcore[] shadow-cone
│                          fields (via the same swi_test_eclipse_where_dcore hook sol_where uses,
│                          already generic over ipl/starname) eps 1e-7 + exact retval flags
│                          bitmask; occ_when_glob: 6 cases — same 3 occulted bodies ×
│                          2 backward, tjd_start=2451545.0, ifltype=0, asserts tret[0..7] eps 1e-5
│                          day + exact retval flags bitmask; both passed on the first
│                          implementation attempt (after fixing gen_eclipse.c's star-name buffer:
│                          swe_fixstar strcpy()s into its `starname` argument in place, so a
│                          string-literal pointer segfaults — must copy into a local `char[]`
│                          first, and `swe_set_ephe_path` needs an explicit path for the same
│                          reason fixstar tests do), no escape-hatch escalation needed;
│                          occ_when_loc: 4 cases — Venus (planet) + Aldebaran (star) × 2 backward,
│                          geopos=(8.55,47.37,500) Zurich, tjd_start=2451545.0; asserts
│                          tret[0..6] eps 1e-5 day + attr[0..7] eps 1e-5 + exact retval flags
│                          bitmask (Venus/forward case exercises OCC_BEG_DAYLIGHT/
│                          OCC_END_DAYLIGHT bits 8192/16384; star cases exercise the
│                          contact-1/4-aliased-from-2/3 point-source branch); FIRST attempt
│                          revealed the pre-existing fixstar-TOPOCTR gap (see eclipse.rs's map
│                          entry) via ~4e-3° azimuth error on the star cases — root-caused by
│                          diffing intermediate tret[]/attr[] against gen_eclipse.c's raw output
│                          per-field, not by re-reading the C; fixed at the calc_fixstar level,
│                          not worked around in the test)
│   ├── pheno.rs       — golden tests for swe_pheno (swisseph-rs/83): 120 cases, Sun..Pluto (10) ×
│                          4 epochs (J2000/2024/1950/1800-Jan-5) × {MOSEPH, MOSEPH|TRUEPOS, SWIEPH},
│                          picks moshier vs sweph Ephemeris by the SWIEPH flag bit, asserts
│                          attr[0..5] eps 1e-9 (Moon magnitude 1e-8) + exact retflag. Exercises
│                          magnitude branches 5a-5j; the Bowell §5k asteroid branch is
│                          golden-uncovered (no backend computes asteroid positions yet — Ceres
│                          omitted). Transitively validates the new SEFLG_HELCTR calc paths
│                          (planets+Moon, Moshier+Swiss) at 1e-9. Pre-1900 epoch nudged to
│                          1800-Jan-5 (off sepl_18's tfstart) to avoid a documented C stateless
│                          file-boundary artifact in swe_pheno's elongation.
│   ├── heliacal_internals.rs — golden tests for heliacal atmospheric/optics layer
│                          (swisseph-rs/104, sub-task 1/8): extinction battery 48 cases
│                          (AltO×AltS×sunra rotated over Lat/HeightEye/datm, asserts Deltam + kt +
│                          kR + kOZ + kW + ka at 1e-12); airmass battery 14 cases (AppAltO×Press,
│                          asserts Airmass + Xext{rayleigh,water,aerosol} + Xlay_ozone); app_alt
│                          battery 36 cases (alt×TempE×PresE, asserts AppAltfromTopoAlt +
│                          TopoAltfromAppAlt at 1e-12); optic battery 20 cases (B×4 dobs configs
│                          {default,age60,binocular,optical_params}, asserts CVA + PupilDia +
│                          OpticFactor{intensity,background} at 1e-12); search battery 13 cases
│                          (swisseph-rs/110, sub-task 7/8): find_conjunct_sun 8 cases (Venus/Mars ×
│                          TypeEvent {1,2} × tjd_start {2453000,2451545}), get_heliacal_day 1 case
│                          (Venus morning first from conjunction seed), time_optimum_visibility 1
│                          case, time_limit_invisible 2 cases (dir=±1), get_heliacal_details 1 case
│                          (Venus te=1), all at Cairo observer, eps 2e-5 day.
│                          C harness includes swehel.c directly for access to static functions
│   ├── heliacal.rs       — golden tests for swe_vis_limit_mag (swisseph-rs/107, sub-task 4/8):
│                          33 cases — 4 objects {venus,sirius,moon,mercury} × per-object UT
│                          instants (daytime for planets, nighttime for stars/Moon) at Cairo
│                          observer (31.25°E, 30.1°N, 30m), + VISLIM_DARK/NOMOON/PHOTOPIC/
│                          SCOTOPIC flag variants, OPTICAL_PARAMS custom dobs, MOSEPH duplicates;
│                          14 photopic (retval 0), 14 scotopic (retval 1), 5 below-horizon
│                          (retval -2). Positions eps 1e-7, limiting magnitude eps 5e-7 (chain
│                          compounds azalt+refraction+brightness+extinction+optics). C harness
│                          gen_heliacal.c links libswe.a (public API, not swehel.c internals).
│                          Also: swe_heliacal_pheno_ut (swisseph-rs/109, sub-task 6/8):
│                          15 cases — Moon crescent × TypeEvent=3 at Mecca (2 young-crescent
│                          evenings + 1 HIGH_PRECISION), Venus × TypeEvent=1/2 at Cairo/Mecca,
│                          Sirius × TypeEvent=1 at Cairo (± HIGH_PRECISION), Mercury morning first,
│                          Moon morning last, Mars/Jupiter evening first (early-exit guard path),
│                          MOSEPH duplicates (Venus + Moon). Geometry slots eps 1e-7, time/duration
│                          slots (TfirstVR/TbVR/TlastVR/TbYallop/RiseO/RiseS/Lag/TvisVR) eps 1e-5.
│                          Also: swe_heliacal_ut (swisseph-rs/111, sub-task 8/8):
│                          12 event cases — Sirius/Venus/Moon at Cairo/Mecca/Athens observers,
│                          TypeEvent 1–4, LONG_SEARCH, AVKIND_VR ×2 (arc_vis path), HIGH_PRECISION;
│                          vis_lim cases eps 2e-5 day, arc_vis cases eps 1e-4 day (coarser search),
│                          arc_vis dret[1]/dret[2] asserted == 0.0
│   ├── obliquity_bias.rs — golden tests for obliquity + bias
│   ├── precession.rs  — golden tests for precession (374 cases)
│   ├── nutation.rs    — golden tests for nutation (80 cases + router tests)
│   ├── deltat.rs      — golden tests for delta-T (217 cases: 5 models × 43 epochs)
│   ├── sidereal_time.rs — golden tests for sidereal time (128 cases: 4 models × 32 epochs)
│   ├── mean_elements.rs — golden tests for mean node, mean apogee, ECL_NUT (231 cases: 165 tropical [3 bodies × 11 epochs × 5 flag combos] + 66 sidereal [2 lunar bodies × 11 epochs × 3 sid_modes {Lahiri, Lahiri|ECL_T0, Lahiri|SSY_PLANE}, SEFLG_SIDEREAL|SEFLG_SPEED] — regression guard for the swisseph-rs/84 review follow-up threading x2000 through mean_element_pipeline; tropical eps 1e-10, sidereal 1e-9 pos / 1e-7 speed)
│   ├── truenode.rs     — golden tests for osculating node/apogee (swisseph-rs/84: 252 cases, 2
│                          bodies {SE_TRUE_NODE, SE_OSCU_APOG} × {MOSEPH, SWIEPH}. 168 tropical: × 6
│                          flag combos {SPEED, SPEED|EQUATORIAL, SPEED|XYZ, SPEED|NONUT, SPEED|J2000,
│                          no_speed} × 7 gen_calc epochs. 84 sidereal (SEFLG_SIDEREAL|SEFLG_SPEED): ×
│                          3 sid_modes {Lahiri traditional, Lahiri|ECL_T0, Lahiri|SSY_PLANE} × 7
│                          epochs — regression guard for the swisseph-rs/84 review fix that threads
│                          x2000 so the ECL_T0/SSY_PLANE rigorous branches no longer silently fall
│                          back to traditional ayanamsa subtraction. Positions eps 1e-9, speeds 1e-7
│                          EXCEPT Moshier node/apogee speed relaxed to 5e-6 for the documented C
│                          global-cache finite-difference artifact — see CLAUDE.md
│                          <stateless_tolerance> §3)
│   ├── nodaps.rs        — golden tests for swe_nod_aps mean (swisseph-rs/85) + osculating
│                          (swisseph-rs/86) branches. mean: 200 Moshier cases — 10 bodies
│                          {Sun/Moon/Mercury..Neptune/Earth} × 4 epochs (incl. pre-1900) × 5 flags
│                          {SPEED, SPEED|EQUATORIAL, no_speed, SPEED|TRUEPOS, SPEED|EQUATORIAL|
│                          TRUEPOS}, asserting all four asc/desc/peri/aphe [f64;6]. Per-point/flag
│                          tolerance (see the `tolerance` fn): TRUEPOS geometry tight (1e-9 pos/1e-8
│                          speed, bit-exact), apparent asc/peri/aphe 1e-6, apparent descending node
│                          relaxed to 1e-3°/2e-2°/day for the C node-distance ill-conditioning
│                          (docs/swisseph-c-potential-bugs.md §7). oscu: 72 cases — 9 bodies
│                          {Moon/Mercury..Neptune/Pluto} × 4 epochs (pre-1900 nudged to 1800-Jan-5,
│                          `osc_epochs` in gen_nodaps.c, off the sepl_18 file boundary) × 2 backends
│                          {MOSEPH, SWIEPH} × SEFLG_SPEED, method=OSCU. oscu_bar: 8 SWIEPH cases —
│                          Jupiter/Saturn/Pluto (beyond 6 AU)/Mercury (inside it) × 2 epochs,
│                          method=OSCU_BAR (Moshier has no real barycenter, rejects with
│                          Error::UnsupportedFlags — untested by this SWIEPH-only battery). fopoint:
│                          6 MOSEPH cases — Moon/Mars/Jupiter × 2 epochs, method=OSCU|FOPOINT.
│                          Shared `check_case`/`ephemeris_for` (backend-keyed Ephemeris cache) helpers.
│                          Osculating tolerance is tiered by point (see `osc_tolerance` fn) —
│                          peri/aphe 5e-5°pos/1e-4°/day, asc 1e-3°/1e-4°/day, desc 2e-3°/3e-2°/day —
│                          for the `fac=z/ż` node-direction division amplifying backend FP noise into
│                          BOTH nodes (docs/swisseph-c-potential-bugs.md §8, verified via an
│                          identical-input formula check, not guessed))
│   ├── orbit.rs         — golden tests for swe_get_orbital_elements + swe_orbit_max_min_true_distance
│                          (swisseph-rs/87, PNOC 6). elements: 130 cases — Mercury..Pluto+Earth ×
│                          {Moshier default/HELCTR/ORBEL_AA(=TOPOCTR bit) × 4 epochs incl. pre-1900,
│                          Swiss default/HELCTR × 2 epochs, MOSEPH|BARYCTR<6-AU-fallback × 5 bodies ×
│                          2 epochs}. Asserts all 17 dret slots; every PLANET matches C bit-for-bit
│                          (residual 0.0), only Earth carries a ~4e-8° artifact from the -(Sun) helio
│                          derivation, so a uniform 1e-6 both accommodates Earth and guards the
│                          planets. Node(3)/arg-peri(4) relaxed to 1e-3 for near-planar orbits
│                          (incl<0.1° isolates Earth from Uranus's 0.77°): those two are
│                          ill-conditioned when the orbit is near-coplanar with the reference
│                          ecliptic, but their sum peri(5) stays tight (same class as nodaps'
│                          descending-node singularity). maxmin: 30 Moshier cases — Mercury/Venus/
│                          Mars/Jupiter/Pluto × 3 epochs × {geocentric two-ellipse, heliocentric};
│                          dmax/dmin 1e-8 (actual ~1e-10 after the 300-iter refine), dtrue 1e-9.
│   ├── pctr.rs           — golden tests for swe_calc_pctr (swisseph-rs/90, PNOC 9). 90 cases total:
│                          6 (ipl,iplctr) pairs {(Sun,Mars),(Earth,Mars),(Jupiter,Mars),(Moon,Venus),
│                          (Saturn,Jupiter),(Earth,Moon)} × 3 epochs {J2000,2460600,2415020.5} × 5
│                          flag combos {MOSEPH×4, SWIEPH×1}. 72 MOSEPH cases assert Err (BARYCTR
│                          unsupported); 18 SWIEPH cases assert positions eps 5e-8 / speeds eps 1e-7
│                          + exact retflag. Position tolerance wider than standard 1e-9 due to pctr's
│                          deflection geometry using earth_helio vs C's earth_bary (worst case 1.44e-8
│                          for Saturn-Jupiter, same architectural cause as the documented stateless
│                          deflection tolerance)
│   ├── crossings.rs      — golden tests for swe_solcross / mooncross / mooncross_node / helio_cross
│                          (swisseph-rs/88, PNOC 7). 66 Moshier cases total. solcross: 18 cases —
│                          x2cross {0,180,359.5} × jd_start {2451500,2440000,2460600} × {ET,UT}.
│                          mooncross: 18 cases — same shape. mooncross_node: 6 cases — jd_start
│                          {2451545,2440000,2460600} × {ET,UT}, asserts jd+lon+lat. helio_cross:
│                          24 cases — ipl {Mercury,Mars,Jupiter} × x2cross {0,120.5} × dir {1,-1}
│                          × {ET,UT}. Crossing times eps 1e-6 day, lon eps 1e-7, node lat eps 5e-9.
│   ├── moshier_backend.rs — golden tests for backend::compute (110 cases: 10 bodies × 11 epochs + Earth zero-check)
│   ├── moshier_moon.rs — golden tests for moshmoon2 (11 cases: Moon at 11 epochs)
│   ├── moshier_planet.rs — golden tests for moshplan2 (81 cases: 9 planets × 9 epochs)
│   ├── se1_header.rs  — golden tests for SE1 file parsing (11 planet metadata fields, byte-order detection on 84 files)
│   ├── sweph_eval.rs  — golden tests for evaluate_body (80 cases: 10 bodies × 8 epochs, bitwise-exact positions + velocities)
│   ├── jpl_pleph.rs   — golden tests for jpl_pleph (84 cases: 11 bodies × 7 epochs barycentric + 7 geocentric Moon, 1e-9 eps)
│   └── houses.rs      — golden tests for houses_armc (angles_special: 30 cases system-independent special points;
│                         equal_family: 150 cases, 5 systems A/D/N/V/W × 30 battery cases, bitwise-exact;
│                         quad_arith: 150 cases, 5 systems O/S/X/M/F × 30 battery cases, eps 1e-9/1e-7;
│                         great_circle: 150 cases, 5 systems R/C/T/H/J × 30 battery cases, eps 1e-9;
│                         iterative: 84 cases, 2 systems P/K × 6 armc × 7 geolat (incl. ±78 polar) ×
│                         1 eps, eps 1e-9 cusps/1e-7 speeds; gauquelin36: 42 cases, G × 6 armc × 7 geolat
│                         × 1 eps, cusps[1..36], same eps; closed_form_misc: 120 cases, 4 systems
│                         U/Y/L/Q × 30 battery cases, eps 1e-9 cusps/1e-7 speeds (U speeds eps 0 —
│                         stale pre-switch values, asserted exactly per c-ref-houses.md §4.2e);
│                         sunshine: 76 cases — 60: 2 systems I/i × 6 armc × 5 geolat (1 sundec per
│                         case, rotated through {-23,-10,0,10,23}), eps 1e-9 cusps (1e-8 for
│                         Makransky 'i')/1e-7 speeds (driver-level finite-diff); + 16: 2 systems ×
│                         2 armc × geolat {70,-70} × sundec {23,-23} (all four combinations satisfy
│                         |tand(geolat)·tand(sundec)|≥1, triggering Makransky's circumpolar ERR →
│                         Porphyry fallback; Treindl never short-circuits on it, included at the
│                         same combos for contrast), eps 1e-9 cusps/speeds for the 'i' fallback
│                         subset (fill_porphyry is closed-form, bitwise-exact elsewhere);
│                         sunshine_requires_sundec: negative test, Sunshine + sundec=None
│                         returns Err; ut_wrapper: 42 cases — 36: 6 (tjd_ut,geolat,geolon) triples
│                         × 6 systems P/K/C/R/W/I + 6: 1 triple × 6 systems with SEFLG_NONUT,
│                         eps 1e-7 cusps/ascmc, 1e-6 speeds via Ephemeris::houses_ex2 (compounds
│                         deltaT/obliquity/nutation/sidtime — looser than the pure-armc tests);
│                         sidereal_trad: 9 cases, systems P/W/E × 3 triples, SEFLG_SIDEREAL +
│                         Lahiri, same tolerances; sidereal_geom: 18 cases, systems P/C/W x 3
│                         triples x 2 sid_modes (Lahiri|SE_SIDBIT_ECL_T0, Lahiri|SE_SIDBIT_SSY_PLANE),
│                         same tolerances, exercises sidereal_houses_ecl_t0/ssypl; house_pos: 150
│                         cases, all 25 house-system chars
│                         × 2 armc/geolat/eps triples (1 temperate, 1 polar — armc=105/geolat=67
│                         chosen so 2/3 xpin succeed and 1/3 hits Koch's genuine hpos==0 circumpolar
│                         failure) × 3 xpin, eps 1e-7 hpos, "err" cases (1: Koch) assert Err instead
│                         (err is driven by hpos==0.0, not serr non-empty — several systems set an
│                         informational serr on a valid nonzero hpos); gauquelin_sector: 36 cases,
│                         6 (tjd_ut,ipl,imeth) combos × 3 bodies (Sun/Moon/Mars) × 2 imeth, eps 1e-6
│                         dgsect via Ephemeris::gauquelin_sector_geometric — this test caught a
│                         pre-existing sidtime_long_term deltaT/tid_acc bug at 1600/1800 AD epochs,
│                         fixed in sidereal_time.rs (see that file's codebase-map entry));
│                         gauquelin_riseset: 78 cases (swisseph-rs/89, PNOC 8), 72 planet: 6
│                         ut_triples × 3 bodies (Sun/Moon/Mars) × 4 imeth {2,3,4,5}, SEFLG_MOSEPH,
│                         atpress=0, attemp=0, eps 1e-6 dgsect via Ephemeris::gauquelin_sector; 1
│                         circumpolar case (Moon at -64° latitude, 1600-Jan-1) asserts Err; + 6
│                         star: Aldebaran × 2 ut_triples × imeth {0,1,2}, exercises the fixstar2
│                         geometric branch + rise/set star path (swisseph-rs/97); star rise/set
│                         eps 1e-5 for the known fixstar-TOPOCTR gap)
├── riseset.rs         — golden tests for rise_trans_true_hor + rise_trans (swisseph-rs/70,71;
│                         full: 36 cases, 3 geopos (Zurich/Null Island/Tromso) × 2 bodies
│                         (Sun/Moon) × 2 epochs × 3 rsmi (RISE/SET/MTRANSIT, all with FORCE_SLOW
│                         OR'd in for parity with the C harness though it's a no-op on this
│                         function), eps 1e-6 time (≈0.1s); 4 of the 36 (Tromso Sun RISE/SET at
│                         both epochs) are circumpolar (C retval -2) and assert
│                         Err(Error::CircumpolarBody) instead of a time; dip: 6 cases, horhgt=-100
│                         auto-dip sentinel × atpress ∈ {0, 1013.25} × 3 geopos, locks in that
│                         calc_dip receives atpress unmodified (not auto-estimated); mtrans_flags:
│                         12 cases, epheflag NONUT|TRUEPOS × 3 geopos × 2 bodies × MTRANSIT/
│                         ITRANSIT, locks in calc_mer_trans's narrower SEFLG_EPHMASK-only mask;
│                         fast: 24 cases via Ephemeris::rise_trans (swe_rise_trans dispatcher),
│                         3 geopos all |lat|≤60 (Zurich/Null Island/Tokyo) × 2 bodies × 2 epochs ×
│                         RISE/SET (no FORCE_SLOW — that's what selects the fast path), eps 1e-6
│                         vs C tret0 + a same-input fast-vs-full cross-check eps 1e-5 day)
│   ├── source_consistency.rs — unit tests for epheflag/EphemerisConfig source consistency
│                          (swisseph-rs/112): 5 cases — SWIEPH flags on Moshier config clamps to
│                          MOSEPH (flags_used), MOSEPH flags on Swiss config produces bitwise-
│                          identical Moshier output, MOSEPH-on-Swiss calc_ut at 1800 AD uses DE404
│                          tid_acc (bitwise match vs pure Moshier), matched-flags identity, no
│                          EPHMASK uses config default. No C golden data — pure Rust-side invariants
├── golden-data/
│   ├── asteroid.json   — C-generated reference data for asteroid calc (300 cases: 6 bodies ×
│                          5 epochs × 10 flag combos; see tests/golden/asteroid.rs)
│   ├── calc.json       — C-generated reference data for calc pipeline (swe_calc full pipeline)
│   ├── corrections.json — C-generated reference data for corrections (meff, aberr_light, pipeline)
│   ├── math.json       — C-generated reference data for math
│   ├── date.json       — C-generated reference data for date
│   ├── eclipse.json    — C-generated reference data for swe_sol_eclipse_where (sol_where key),
│                          swe_sol_eclipse_how (sol_how key), swe_sol_eclipse_when_glob
│                          (sol_when_glob key), swe_sol_eclipse_when_loc (sol_when_loc key),
│                          swe_lun_eclipse_how (lun_how key), swe_lun_eclipse_when (lun_when key),
│                          swe_lun_eclipse_when_loc (lun_when_loc key), swe_lun_occult_where
│                          (occ_where key), swe_lun_occult_when_glob (occ_when_glob key), and
│                          swe_lun_occult_when_loc (occ_when_loc key, swisseph-rs/79) — all keys
│                          this module needs now, final RSE arc task
│   ├── pctr.json       — C-generated reference data for swe_calc_pctr (90 cases; see tests/golden/pctr.rs)
│   ├── pheno.json      — C-generated reference data for swe_pheno (120 cases; see tests/golden/pheno.rs)
│   ├── obliquity_bias.json — C-generated reference data for obliquity/bias
│   ├── precession.json — C-generated reference data for precession
│   ├── nutation.json   — C-generated reference data for nutation
│   ├── deltat.json     — C-generated reference data for delta-T
│   ├── sidereal_time.json — C-generated reference data for sidereal time
│   ├── mean_elements.json — C-generated reference data for mean node, mean apogee, ECL_NUT (231 cases: 165 tropical + 66 sidereal)
│   ├── truenode.json   — C-generated reference data for swe_calc SE_TRUE_NODE/SE_OSCU_APOG (252 cases: 168 tropical + 84 sidereal)
│   ├── nodaps.json     — C-generated reference data for swe_nod_aps mean branch (key "mean", 200 cases; PNOC 5 adds oscu/oscu_bar/fopoint keys)
│   ├── moshier_backend.json — C-generated reference data for backend::compute (swe_calc with ICRS)
│   ├── moshier_moon.json — C-generated reference data for moshmoon2
│   ├── moshier_planet.json — C-generated reference data for moshplan2
│   ├── se1_header.json — C-generated reference data for SE1 file headers (sepl_18, semo_18)
│   ├── sweph_eval.json — C-generated reference data for evaluate_body (raw Chebyshev eval + rot_back + ecl→equ rotation)
│   ├── jpl_pleph.json  — C-generated reference data for jpl_pleph (84 cases via swi_pleph against de441.eph)
│   ├── fixstar.json    — C-generated reference data for swe_fixstar2 (196 position cases + 4 mag cases, 7 stars × 4 epochs × 7 flags)
│   ├── azalt.json      — C-generated reference data for swe_refrac/swe_refrac_extended/swe_azalt/swe_azalt_rev (refrac: 28, refrac_ext: 56, azalt: 8, azalt_rev: 8)
│   ├── houses.json     — C-generated reference data for swe_houses_armc_ex2 (battery: 6 armc × 5 geolat × 1 eps, reused across all houses sub-tasks; iterative/gauquelin36 keys add a 7th/8th polar geolat (±78) to exercise the Placidus/Koch/Gauquelin Porphyry fallback; closed_form_misc key reuses the standard 5-geolat battery for U/Y/L/Q; sunshine key reuses the standard 6 armc × 5 geolat battery for I/i, crossed with a rotated (not full cross-product) Sun-declination set {-23,-10,0,10,23}, plus a dedicated circumpolar-Sun sub-battery (geolat {70,-70} × sundec {23,-23}) to exercise Makransky's ERR→Porphyry fallback; ut_wrapper key: swe_houses_ex2 (UT-based) over 6 (tjd_ut,geolat,geolon) triples × 6 systems, + a SEFLG_NONUT variant at 1 triple; sidereal_trad key: swe_houses_ex2 with SEFLG_SIDEREAL + swe_set_sid_mode(SE_SIDM_LAHIRI) over 3 triples × 3 systems P/W/E; house_pos key: swe_house_pos over all 25 house-system chars × 2 (armc,geolat,eps) triples × 3 xpin, "err" field is hpos==0.0 (Koch's real failure sentinel), NOT serr-non-empty (P/G/J/L/Q/default set an informational serr on valid results) — the static sundec cache 'I'/'i' need is primed via a preceding swe_houses_armc_ex2(ascmc[9]=sundec) call; gauquelin_sector key: swe_gauquelin_sector imeth∈{0,1} over 6 ut_triples × 3 bodies (Sun/Moon/Mars); gauquelin_riseset key (swisseph-rs/89, PNOC 8): swe_gauquelin_sector imeth∈{2,3,4,5} over 6 ut_triples × 3 bodies × 4 imeth = 72 cases, retval recorded for circumpolar Err)
│   ├── riseset.json    — C-generated reference data for swe_rise_trans_true_hor + swe_rise_trans (full key: 36 cases, 3 geopos × 2 bodies × 2 epochs × 3 rsmi, retval recorded so circumpolar -2 cases assert Err; dip key: 6 cases, horhgt=-100 × atpress∈{0,1013.25} × 3 geopos; mtrans_flags key: 12 cases, NONUT|TRUEPOS × 3 geopos × 2 bodies × MTRANSIT/ITRANSIT; fast key: 24 cases via swe_rise_trans, 3 geopos all \|lat\|≤60 × 2 bodies × 2 epochs × RISE/SET, no FORCE_SLOW)
│   ├── heliacal.json   — C-generated reference data for swe_vis_limit_mag (vis_limit key, 33 cases),
│                          swe_topo_arcus_visionis (arcvis key), swe_heliacal_angle (helangle key),
│                          swe_heliacal_pheno_ut (pheno key, 15 cases; swisseph-rs/109),
│                          swe_heliacal_ut (events key, 12 cases; swisseph-rs/111);
│                          see tests/golden/heliacal.rs
│   └── crossings.json  — C-generated reference data for swe_solcross/mooncross/mooncross_node/helio_cross (66 cases: 18 solcross + 18 mooncross + 6 mooncross_node + 24 helio_cross, all Moshier)
└── c-gen/
    ├── gen_asteroid.c  — C harness to regenerate asteroid.json (6 bodies × 5 epochs × 10 flags
    │                       = 300 cases; swe_close+swe_set_ephe_path("ephe")+swe_set_topo per call;
    │                       no MOSEPH cases per decision 1; aborts on any swe_calc error)
    ├── gen_calc.c      — C harness to regenerate calc.json (full swe_calc pipeline, 14 bodies × 7 epochs × 12 flags, ECL_NUT cleanup)
    ├── gen_eclipse.c   — C harness to regenerate eclipse.json (swe_sol_eclipse_where: 3 known
    │                       central eclipses at their real maximum-eclipse UT instants + 1
    │                       no-eclipse epoch, one case with SEFLG_NONUT to confirm it's masked
    │                       away by swe_sol_eclipse_where's own `ifl &= SEFLG_EPHMASK`;
    │                       swe_sol_eclipse_how: the same 4 epochs × 2 observers (near-central
    │                       and off-track); swe_sol_eclipse_when_glob: 2 tjd_start × 2 backward,
    │                       ifltype=0; swe_sol_eclipse_when_loc: 2 geopos (near-central; Chile) ×
    │                       2 tjd_start × 2 backward; swe_lun_eclipse_how: 3 known lunar eclipses
    │                       (2001/2021/2024) at their real maximum-eclipse UT instants × 1 observer
    │                       (near-central Zurich); swe_lun_eclipse_when: 2 tjd_start (2000/2019) ×
    │                       2 backward, ifltype=0; swe_lun_eclipse_when_loc: 1 geopos (near-central
    │                       Zurich) × the same 2 tjd_start × 2 backward; swe_lun_occult_where +
    │                       swe_lun_occult_when_glob: 3 occulted bodies (Venus, Mars, Aldebaran
    │                       starname) at tjd_ut/tjd_start=2451545.0 (where) or 2458800.5
    │                       (when_glob) × 2 backward (when_glob only); swe_lun_occult_when_loc
    │                       (swisseph-rs/79): 2 of those 3 bodies (Venus, Aldebaran — skips Mars
    │                       to keep the battery small) × 1 geopos (near-central Zurich) × 1
    │                       tjd_start=2451545.0 × 2 backward)
    ├── gen_pheno.c     — C harness to regenerate pheno.json (swe_pheno: Sun..Pluto × 4 epochs ×
    │                       {MOSEPH, MOSEPH|TRUEPOS, SWIEPH}; asteroids omitted, boundary-safe 1800 epoch)
    ├── gen_mean_elements.c — C harness to regenerate mean_elements.json (mean node, mean apogee, ECL_NUT; 165 tropical + 66 sidereal [MeanNode/MeanApogee × 11 epochs × 3 sid_modes via swe_set_sid_mode])
    ├── gen_truenode.c  — C harness to regenerate truenode.json (swe_calc SE_TRUE_NODE/SE_OSCU_APOG, 2 bodies × {MOSEPH,SWIEPH}; 168 tropical [6 flags × 7 epochs] + 84 sidereal [3 sid_modes {Lahiri, Lahiri|ECL_T0, Lahiri|SSY_PLANE} × 7 epochs, via swe_set_sid_mode]; sets ephe path for SWIEPH)
    ├── gen_nodaps.c    — C harness to regenerate nodaps.json (swe_nod_aps, SE_NODBIT_MEAN: 10 bodies × 4 epochs × 5 flags = 200 cases; 2 TRUEPOS combos for tight geometry assertion)
    ├── gen_corrections.c — C harness to regenerate corrections.json (meff copied from sweph.c, swi_aberr_light direct, pipeline via swe_calc)
    ├── gen_obliquity_bias.c — C harness to regenerate obliquity_bias.json
    ├── gen_precession.c — C harness to regenerate precession.json
    ├── gen_nutation.c  — C harness to regenerate nutation.json
    ├── gen_deltat.c    — C harness to regenerate deltat.json
    ├── gen_sidereal_time.c — C harness to regenerate sidereal_time.json
    ├── gen_moshier_backend.c — C harness to regenerate moshier_backend.json (swe_calc with all corrections disabled + ICRS)
    ├── gen_moshier_moon.c — C harness to regenerate moshier_moon.json
    ├── gen_moshier_planet.c — C harness to regenerate moshier_planet.json
    ├── gen_sweph_eval.c — C harness to regenerate sweph_eval.json (raw SE1 Chebyshev eval via swed.pldat internals)
    ├── gen_se1_header.c — standalone binary parser, dumps header + planet metadata as JSON
    ├── gen_jpl_pleph.c  — C harness to regenerate jpl_pleph.json (swi_pleph direct calls against de441.eph)
    ├── gen_fixstar.c    — C harness to regenerate fixstar.json (swe_fixstar2: 7 stars × 4 epochs × 7 flags + 4 mag cases)
    ├── gen_azalt.c      — C harness to regenerate azalt.json (swe_refrac 7 inalt × 2 atpress × 2 dir; swe_refrac_extended × 2 geoalt; swe_azalt/swe_azalt_rev 2 tjd_ut × 2 geopos × 2 dir; swe_set_ephe_path(NULL))
    ├── gen_houses.c     — C harness to regenerate houses.json (swe_houses_armc_ex2: angles_special,
    │                       equal_family, quad_arith, great_circle, iterative, gauquelin36,
    │                       closed_form_misc, sunshine — sunshine key sets ascmc[9]=sundec before
    │                       calling, per c-ref-houses.md §11; ut_wrapper/sidereal_trad keys use
    │                       swe_houses_ex2 (UT-based) instead, over a 6-triple
    │                       (tjd_ut,geolat,geolon) battery; house_pos key: swe_house_pos over all
    │                       25 house-system chars, primes the 'I'/'i' static sundec cache via a
    │                       swe_houses_armc_ex2(ascmc[9]=sundec) call immediately before each
    │                       swe_house_pos call; gauquelin_sector key: swe_gauquelin_sector
    │                       imeth∈{0,1} reusing the ut_wrapper triples × 3 bodies;
    │                       gauquelin_riseset key (swisseph-rs/89): swe_gauquelin_sector
    │                       imeth∈{2,3,4,5} × same triples × 3 bodies = 72 cases)
    ├── gen_riseset.c    — C harness to regenerate riseset.json (swe_rise_trans_true_hor: full
    │                       key, 3 geopos × 2 bodies (Sun/Moon) × 2 epochs × 3 rsmi
    │                       (RISE/SET/MTRANSIT, | SE_BIT_FORCE_SLOW_METHOD), SEFLG_MOSEPH,
    │                       records retval so circumpolar -2 cases assert Err; dip key: horhgt=
    │                       -100 × atpress∈{0,1013.25} × 3 geopos; mtrans_flags key: NONUT|
    │                       TRUEPOS × 3 geopos × 2 bodies × MTRANSIT/ITRANSIT; fast key:
    │                       swe_rise_trans (swisseph-rs/71), 3 geopos all |lat|≤60 (Zurich/Null
    │                       Island/Tokyo) × 2 bodies × 2 epochs × RISE/SET, no FORCE_SLOW)
    ├── gen_heliacal.c   — C harness to regenerate heliacal.json (swe_vis_limit_mag: 4 objects ×
    │                       per-object UT instants at Cairo observer, + VISLIM_DARK/NOMOON/PHOTOPIC/
    │                       SCOTOPIC flag variants, OPTICAL_PARAMS, MOSEPH duplicates = 33 cases;
    │                       swe_heliacal_pheno_ut (swisseph-rs/109): 15 cases — Moon crescent at
    │                       Mecca, Venus/Sirius/Mercury/Mars/Jupiter at Cairo, MOSEPH duplicates;
    │                       swe_heliacal_ut (swisseph-rs/111): 12 event cases — Sirius/Venus/Moon
    │                       at Cairo/Mecca/Athens, TypeEvent 1–4, LONG_SEARCH, AVKIND_VR ×2,
    │                       HIGH_PRECISION; links libswe.a normally, does NOT #include swehel.c)
    └── gen_cross.c      — C harness to regenerate crossings.json (66 cases: solcross 18,
                            mooncross 18, mooncross_node 6, helio_cross 24, all SEFLG_MOSEPH)
```

## Key Types in types.rs

### Astronomical model enums (lines 517–597)

| Type | Lines | Variants | repr |
|---|---|---|---|
| `PrecessionModel` | 522–534 | IAU1976=1..Newcomb=11 (11 variants) | i32 |
| `NutationModel` | 538–544 | IAU1980=1..Woolard=5 | i32 |
| `DeltaTModel` | 548–554 | 5 variants | i32 |
| `SiderealTimeModel` | 558–563 | IAU1976=1..Longterm=4 | i32 |
| `BiasModel` | 567–571 | None=1, IAU2000=2, IAU2006=3 | i32 |
| `JplHorMode` | 575–577 | LongAgreement=1 | i32 |
| `JplHoraMode` | 581–585 | V1=1, V2=2, V3=3 | i32 |
| `AstroModels` | 588–597 | 8 fields: delta_t, prec_longterm, prec_shortterm, nutation, bias, jplhor_mode, jplhora_mode, sidereal_time |
| `Default` impl | 680–693 | longterm=shortterm=Vondrak2011, nutation=IAU2000B, bias=IAU2006, jplhora=V3, sidereal=Longterm |

### Frame/obliquity types (after AstroModels, before JdTt)

| Type | Lines | Notes |
|---|---|---|
| `FrameTransform` | ~599 | J2000ToGcrs, GcrsToJ2000 |
| `PrecessionDirection` | ~603 | J2000ToDate, DateToJ2000 |
| `Epsilon` | ~607 | eps, sin_eps, cos_eps + `Epsilon::new(eps_rad)` |
| `Nutation` | ~628 | dpsi, deps (radians) |

### Julian Day newtypes

| Type | Line | Notes |
|---|---|---|
| `JdTt` | ~624 | newtype `(pub f64)` |
| `JdUt1` | ~627 | newtype `(pub f64)` |

### Other key types

| Type | Lines | Notes |
|---|---|---|
| `Body` | 82–111 | enum: Sun..Vesta + Fictitious/Asteroid/PlanetMoon/Comet |
| `HouseSystem` | 272–298 | 22 house system variants |
| `CalendarType` | 373–376 | Julian, Gregorian |
| `SiderealMode` | 396–445 | 42 sidereal mode variants |
| `EphemerisSource` | 510–514 | Jpl, Swisseph, Moshier |
| `UtcComponents` | 640–647 | year, month, day, hour, min, sec |
| `DeltaT` | 659–661 | trait |
| `DegreeParts` | 668–674 | degrees, minutes, seconds, second_fraction, sign |

## Flags (src/flags.rs)

Six `bitflags!` structs. Most relevant:

**CalcFlags (u32)** — lines 3–31:
- JPLEPH=1, SWIEPH=2, MOSEPH=4, HELCTR=8, TRUEPOS=16, J2000=32
- NONUT=64, SPEED3=128, SPEED=256, NOGDEFL=512, NOABERR=1024
- EQUATORIAL=2048, XYZ=4096, RADIANS=8192, BARYCTR=16384
- TOPOCTR=32768, SIDEREAL=65536, ICRS=131072
- **DPSIDEPS_1980=262144** (C: SEFLG_JPLHOR)
- **JPLHOR_APPROX=524288**
- CENTER_BODY=1048576

## Constants (src/constants.rs)

Key constants for quick reference:

| Name | Value | Line |
|---|---|---|
| DEGTORAD | π/180 | 29 |
| STR | 4.8481368e-6 (arcsec→rad) | 38 |
| TWOPI | 2π | 34 |
| J2000 | 2451545.0 | 59 |
| B1950 | 2433282.42345905 | 60 |
| J1900 | 2415020.0 | 61 |
| B1850 | 2396758.2035810 | 62 |
| PREC_IAU_1976_CTIES | 2.0 | ~64 |
| PREC_IAU_2000_CTIES | 2.0 | ~65 |
| PREC_IAU_2006_CTIES | 75.0 | ~66 |
| DPSI_DEPS_IAU1980_TJD0_HORIZONS | 2437684.5 | ~70 |

## Math Functions (src/math.rs)

All `pub fn`. Key functions and their line ranges:

| Function | Lines | Signature |
|---|---|---|
| normalize_degrees | 11–20 | (f64) → f64 |
| normalize_radians | 22–31 | (f64) → f64 |
| mod_2pi | 33–39 | (f64) → f64 |
| mods3600 | 41–43 | (f64) → f64 — arcsec modulo 1296000 |
| diff_degrees_norm | 49–51 | (f64, f64) → f64 |
| diff_degrees | 49–52 | (f64, f64) → f64 |
| diff_radians | 54–57 | (f64, f64) → f64 |
| midpoint_degrees | 63–66 | (f64, f64) → f64 |
| midpoint_radians | 68–70 | (f64, f64) → f64 |
| csnorm | 76–93 | (i32) → i32 |
| difcsn | 95–97 | (i32, i32) → i32 |
| difcs2n | 99–102 | (i32, i32) → i32 |
| d2l | 108–114 | (f64) → i32 |
| chebyshev_eval | 120–131 | (f64, &[f64]) → f64 |
| chebyshev_deriv | 133–156 | (f64, &[f64]) → f64 |
| cross_prod | 167–173 | ([f64;3], [f64;3]) → [f64;3] — swi_cross_prod |
| dot_prod_unit | 177–185 | ([f64;3], [f64;3]) → f64 — swi_dot_prod_unit, clamped to [-1,1] |
| rotate_x | 190–194 | ([f64;3], f64) → [f64;3] |
| rotate_x_sincos | 196–202 | ([f64;3], f64, f64) → [f64;3] |
| cartesian_to_polar | 204–221 | ([f64;3]) → [f64;3] |
| polar_to_cartesian | 223–230 | ([f64;3]) → [f64;3] |
| cartesian_to_polar_with_speed | 236–265 | ([f64;6]) → [f64;6] |
| polar_to_cartesian_with_speed | 267–289 | ([f64;6]) → [f64;6] |
| cotrans | 295–302 | ([f64;3], f64) → [f64;3] |
| cotrans_with_speed | 304–318 | ([f64;6], f64) → [f64;6] |
| split_degrees | 332–392 | (f64, SplitDegFlags) → DegreeParts |
| poly_eval | 394–396 | (&[f64], f64) → f64 — Horner's method |
| OWEN_T0S | 402 | [f64; 5] — Owen interval boundaries |
| owen_t0_icof | 404 | (f64) → (f64, usize) — Owen interval + index |
| owen_chebyshev_basis | 418 | (f64) → (usize, [f64; 10]) — shared by obliquity + precession |
| find_maximum | ~427 | (f64,f64,f64,f64) → (f64,f64) — parabola extremum; offset relative to the rightmost (`y2`) sample, not the middle one; shared by riseset.rs + future eclipse contact-time refinement |
| find_zero | ~446 | (f64,f64,f64,f64) → Option<(f64,f64)> — parabola root(s), same offset convention; `None` on negative discriminant |
| **unit tests** | 438+ | |

## Golden Test Pattern

1. JSON data in `tests/golden-data/<name>.json` — top-level object, keys are test groups, values are arrays of case structs
2. Rust test file in `tests/golden/<name>.rs`:
   - `#[derive(Deserialize)]` case structs matching JSON fields
   - Top-level struct aggregating all case vectors
   - `fn load() -> TopStruct` using `super::golden_data_path()`
   - `#[test]` per group calling `super::assert_f64_exact()` or `super::assert_f64_eps()`
3. Add `mod <name>;` to `tests/golden/main.rs`
4. C generator in `tests/c-gen/` compiled against `../swisseph/libswe.a`

## C Source Reference

- C repo: `../swisseph/`
- C catalogues: `../swisseph/claude/catalogue-{public,internal}.md`
- Ephemeris data: `../swisseph/ephe/`
- Pre-built library: `../swisseph/libswe.a`
- Key headers: `swephexp.h` (public API, model enums), `swephlib.h` (internal functions, PREC constants), `sweph.h` (internal constants)
- `swe_set_astro_models(char *samod, ...)` — comma-separated string: "delta_t,prec_long,prec_short,nutation,bias,jplhor_mode,jplhora_mode,sidt"

## Floating-Point Fidelity Notes

- C polynomial models compute `result * DEGTORAD / 3600` (two runtime ops). Rust must match: `poly_eval(...) * DEGTORAD / 3600.0` — NOT `* STR`.
- C `eps *= DEGTORAD/3600.0` (Laskar, Vondrák) folds to single multiply. Rust: `* (DEGTORAD / 3600.0)` with parens.
- Owen 1990 returns degrees (not arcsec): multiply by `DEGTORAD`, no `/3600`.
- Vondrák `swi_ldp_peps` returns radians directly via `* AS2R` = `* (DEGTORAD / 3600.0)`.
- **`+=` vs `= x +`**: C's `L = L + a + b + c` accumulates left-to-right with L in each step. Rust's `l += a + b + c` evaluates `a + b + c` first, then adds to l. When L is large (~481k) and corrections are small (~6), the different accumulation order produces ULP-level rounding differences that propagate through backward-difference velocity (÷1e-4) and deflection speed (÷5e-7). Always use `l = l + ...` to match C's evaluation order.
- **Multiplication order matters**: `2.0 * x * DEGTORAD` ≠ `2.0 * DEGTORAD * x` due to FP non-associativity. Match C's grouping exactly.

## Insertion Points for New Modules

| What | Where | After |
|---|---|---|
| New model enums | src/types.rs | ~585 (after JplHoraMode) |
| New shared types | src/types.rs | after AstroModels block |
| New AstroModels field | src/types.rs | line 596 + update Default at 681 |
| New constants | src/constants.rs | after existing epoch block |
| New `pub mod` | src/lib.rs | alphabetical in lines 7–21 |
| New re-export | src/lib.rs | inside `pub use types::{...}` block |
| New golden test mod | tests/golden/main.rs | after existing mod declarations |
