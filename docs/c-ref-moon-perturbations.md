# C Reference: Moon Perturbations — `swemmoon.c` (DE404 branch)

**Source file**: `../swisseph/swemmoon.c`  
**Branch**: `#else` of `#ifdef MOSH_MOON_200` — this is the DE404 version used by Swiss Ephemeris (MOSH_MOON_200 is never defined).  
**Goal**: Complete, bitwise-exact porting reference for `moon1()`, `moon2()`, `mean_elements()`, `mean_elements_pl()`.

---

## Thread-Local Variable Map

All C thread-locals (`TLS`) are declared at file scope (~lines 811–843):

```c
static TLS double ss[5][8];   // sin lookup: ss[angle_index][harmonic-1]
static TLS double cc[5][8];   // cos lookup: cc[angle_index][harmonic-1]
static TLS double l;          // Moon ecliptic longitude accumulator (arcsec×10^-5)
static TLS double B;          // Ecliptic latitude accumulator
static TLS double moonpol[3]; // [lon, lat, rad] accumulators
static TLS double SWELP;      // Mean longitude of moon (arcsec)
static TLS double M;          // Mean anomaly of sun (arcsec)
static TLS double MP;         // Mean anomaly of moon = l (arcsec)
static TLS double D;          // Mean elongation of moon (arcsec)
static TLS double NF;         // Mean distance from ascending node = F (arcsec)
static TLS double T;          // Julian centuries from J2000
static TLS double T2;         // T²
static TLS double T3;         // T³ (set externally, not used in moon1/2)
static TLS double T4;         // T⁴ (set externally, not used in moon1/2)
static TLS double f;          // local temporary (planetary combination, persists into moon2!)
static TLS double g;          // current angle (radians)
static TLS double Ve;         // Venus mean longitude (arcsec)
static TLS double Ea;         // Earth mean longitude (arcsec)
static TLS double Ma;         // Mars mean longitude (arcsec)
static TLS double Ju;         // Jupiter mean longitude (arcsec)
static TLS double Sa;         // Saturn mean longitude (arcsec)
static TLS double cg;         // cos(g) temporary
static TLS double sg;         // sin(g) temporary
static TLS double l1;         // t^1 longitude accumulator
static TLS double l2;         // t^2 longitude accumulator
static TLS double l3;         // t^3 longitude accumulator
static TLS double l4;         // t^4 longitude accumulator
```

**Rust mapping**:
- `T`, `T2`, `T4` → state fields
- `SWELP`, `M`, `MP`, `D`, `NF` → `MeanElements` struct (all arcsec)
- `Ve`, `Ea`, `Ma`, `Ju`, `Sa` → `PlanetaryElements` struct (all arcsec)
- `l`, `l1`, `l2`, `l3`, `l4` → longitude accumulator fields
- `B` → latitude accumulator field
- `moonpol[3]` → `[lon, lat, rad]` accumulator array
- `ss[5][8]`, `cc[5][8]` → sin/cos lookup tables
- `f`, `g`, `cg`, `sg` → local temporaries
- `STR = 4.8481368110953599359e-6` (arcseconds to radians)

**Critical**: `f = 18*Ve - 16*Ea` is set in `moon1()` and **reused** in `moon2()` — it is a file-scope thread-local, not a local variable. The Rust port must preserve this across function boundaries.

---

## Call sequence (from `swi_moshmoon2`, line 848)

```c
T = (J - J2000) / 36525.0;
T2 = T * T;
mean_elements();
mean_elements_pl();
moon1();
moon2();
moon3();   // applies all accumulators: l += Horner(l4..l1)*T*1e-5, assembles moonpol
moon4();   // nutation, not in scope here
```

### `moon3()` — how the accumulators are consumed (line 1444)

```c
static void moon3(void) {
    moonpol[0] = 0.0;
    chewm( LR, NLR, 4, 1, moonpol );   // main longitude/radius series (T^0)
    chewm( MB, NMB, 4, 3, moonpol );   // main latitude series (T^0)
    // Horner polynomial evaluation:
    l += (((l4 * T + l3) * T + l2) * T + l1) * T * 1.0e-5;
    moonpol[0] = SWELP + l + 1.0e-4 * moonpol[0];   // arcsec
    moonpol[1] = 1.0e-4 * moonpol[1] + B;           // arcsec
    moonpol[2] = 1.0e-4 * moonpol[2] + 385000.52899; // km
}
```

---

## DE404 `z[]` Array (lines 284–313)

25 elements, z[0]..z[24]. Fitted from DE404 data, -3000 to +3000.

```c
static const double z[] = {
/* Replacements for higher-degree secular terms in mean elements.
   Units: arc seconds; time in Julian centuries. */
-1.312045233711e+01, /* z[0]  F (NF),  t^2 coefficient */
-1.138215912580e-03, /* z[1]  F (NF),  t^3 coefficient */
-9.646018347184e-06, /* z[2]  F (NF),  t^4 coefficient */
 3.146734198839e+01, /* z[3]  l (MP),  t^2 coefficient */
 4.768357585780e-02, /* z[4]  l (MP),  t^3 coefficient */
-3.421689790404e-04, /* z[5]  l (MP),  t^4 coefficient */
-6.847070905410e+00, /* z[6]  D,       t^2 coefficient */
-5.834100476561e-03, /* z[7]  D,       t^3 coefficient */
-2.905334122698e-04, /* z[8]  D,       t^4 coefficient */
-5.663161722088e+00, /* z[9]  L (SWELP), t^2 coefficient */
 5.722859298199e-03, /* z[10] L (SWELP), t^3 coefficient */
-8.466472828815e-05, /* z[11] L (SWELP), t^4 coefficient */
/* Planetary perturbation t^2 amplitudes. Units: arcsec × 10^5 */
-8.429817796435e+01, /* z[12] t^2 cos(18V - 16E - l) */
-2.072552484689e+02, /* z[13] t^2 sin(18V - 16E - l) */
 7.876842214863e+00, /* z[14] t^2 cos(10V - 3E - l)  */
 1.836463749022e+00, /* z[15] t^2 sin(10V - 3E - l)  */
-1.557471855361e+01, /* z[16] t^2 cos(8V - 13E)      */
-2.006969124724e+01, /* z[17] t^2 sin(8V - 13E)      */
 2.152670284757e+01, /* z[18] t^2 cos(4E - 8M + 3J)  */
-6.179946916139e+00, /* z[19] t^2 sin(4E - 8M + 3J)  */
-9.070028191196e-01, /* z[20] t^2 cos(18V - 16E)     */
-1.270848233038e+01, /* z[21] t^2 sin(18V - 16E)     */
-2.145589319058e+00, /* z[22] t^2 cos(2J - 5S)       */
 1.381936399935e+01, /* z[23] t^2 sin(2J - 5S)       */
/* T^3 terms */
-1.999840061168e+00, /* z[24] t^3 sin(l') = t^3 sin(M) */
};
```

**Note**: The MOSH_MOON_200 z[] has 71 elements; the DE404 z[] has only 25. The index mapping is completely different between the two branches. Never mix them.

---

## `sscc()` — Sin/Cos Lookup Table Builder (line 1696)

```c
static void sscc(int k, double arg, int n) {
    double cu, su, cv, sv, s;
    int i;
    su = sin(arg);
    cu = cos(arg);
    ss[k][0] = su;    // sin(1·arg)
    cc[k][0] = cu;    // cos(1·arg)
    sv = 2.0 * su * cu;
    cv = cu * cu - su * su;
    ss[k][1] = sv;    // sin(2·arg)
    cc[k][1] = cv;    // cos(2·arg)
    for (i = 2; i < n; i++) {
        s  = su * cv + cu * sv;
        cv = cu * cv - su * sv;
        sv = s;
        ss[k][i] = sv;  // sin((i+1)·arg)
        cc[k][i] = cv;  // cos((i+1)·arg)
    }
}
```

`ss[k][j]` = sin((j+1) · arg), `cc[k][j]` = cos((j+1) · arg).

---

## `chewm()` — Perturbation Table Evaluator (line 1628)

```c
static void chewm(const short *pt, int nlines, int nangles, int typflg, double *ans) {
    // For each row in the table:
    //   Read nangles integer multipliers; combine sin/cos via trig addition formulas
    //   using ss[m][k-1], cc[m][k-1] for the k-th harmonic of angle m
    // Then accumulate into ans[] based on typflg:
    //   typflg=1 (LR, large):  ans[0] += (10000*j+k)*sv; ans[2] += (10000*j+k)*cv
    //   typflg=2 (LRT2, small): ans[0] += j*sv; ans[2] += k*cv
    //   typflg=3 (MB, large lat): ans[1] += (10000*j+k)*sv
    //   typflg=4 (BT/BT2, lat):  ans[1] += j*sv
}
```

Table index mapping for `sscc`/`chewm` angle indices (0-based):
- Index 0 = D (elongation), harmonics up to 6
- Index 1 = M (solar anomaly), harmonics up to 4
- Index 2 = MP (lunar anomaly = l), harmonics up to 4
- Index 3 = NF (node distance = F), harmonics up to 4

---

## `moon1()` — DE404 Version (lines 1182–1364)

This is the `#else` branch of `#ifdef MOSH_MOON_200`.

### Initialization

```c
static void moon1(void) {
double a;

/* Bug fix (Bhanu Pinnamaneni, 17-aug-2009): zero ss/cc before use.
   Without this, random values could leak in from prior calls. */
int i, j;
for (i = 0; i < 5; i++) {
    for (j = 0; j < 8; j++) {
        ss[i][j] = 0;
        cc[i][j] = 0;
    }
}

/* Build sin/cos lookup tables */
sscc( 0, STR*D,  6 );   // ss[0][0..5], cc[0][0..5] = sin/cos(kD), k=1..6
sscc( 1, STR*M,  4 );   // ss[1][0..3], cc[1][0..3] = sin/cos(kM), k=1..4
sscc( 2, STR*MP, 4 );   // ss[2][0..3], cc[2][0..3] = sin/cos(kMP), k=1..4
sscc( 3, STR*NF, 4 );   // ss[3][0..3], cc[3][0..3] = sin/cos(kNF), k=1..4

moonpol[0] = 0.0;  // longitude accumulator (T^2 scale, units: arcsec×10^-5)
moonpol[1] = 0.0;  // latitude accumulator  (T^2 scale, units: arcsec×10^-5)
moonpol[2] = 0.0;  // radius accumulator    (T^2 scale, units: km×10^-5)
```

### Phase A: T² Series (chewm + planetary perturbations)

```c
/* Tabulated T^2 terms, scale 1.0 = 10^-5" or km */
chewm( LRT2, NLRT2, 4, 2, moonpol ); // 25 longitude+radius rows → moonpol[0], moonpol[2]
chewm( BT2,  NBT2,  4, 4, moonpol ); // 12 latitude rows          → moonpol[1]
```

#### Planetary perturbation terms — Phase A

Units: `l`, `l1` are in arcsec×10^-5; `l2` is in arcsec×10^5 (DE404 z[] scale).  
`moonpol[2]` radius corrections are in km×10^-5.

```c
f = 18 * Ve - 16 * Ea;   // combination used repeatedly; PERSISTS INTO moon2()

/* --- Term 1: angle = 18V - 16E - l --- */
g = STR * (f - MP);
cg = cos(g);
sg = sin(g);
l  = 6.367278 * cg + 12.747036 * sg;   // t^0 longitude (ASSIGNMENT, not +=)
l1 = 23123.70 * cg - 10570.02 * sg;    // t^1 longitude (ASSIGNMENT, not +=)
l2 = z[12] * cg + z[13] * sg;          // t^2 longitude (ASSIGNMENT, not +=)
moonpol[2] += 5.01 * cg + 2.72 * sg;  // radius correction (×10^-5 km)

/* --- Term 2: angle = 10V - 3E - l --- */
g = STR * (10.*Ve - 3.*Ea - MP);
cg = cos(g);
sg = sin(g);
l  += -0.253102 * cg + 0.503359 * sg;
l1 += 1258.46 * cg + 707.29 * sg;
l2 += z[14] * cg + z[15] * sg;

/* --- Term 3: angle = 8V - 13E --- */
g = STR * (8.*Ve - 13.*Ea);
cg = cos(g);
sg = sin(g);
l  += -0.187231 * cg - 0.127481 * sg;
l1 += -319.87 * cg - 18.34 * sg;
l2 += z[16] * cg + z[17] * sg;

a = 4.0*Ea - 8.0*Ma + 3.0*Ju;   // reused for terms 4, 5, and later in Phase B

/* --- Term 4: angle = 4E - 8M + 3J --- */
g = STR * a;
cg = cos(g);
sg = sin(g);
l  += -0.866287 * cg + 0.248192 * sg;
l1 += 41.87 * cg + 1053.97 * sg;
l2 += z[18] * cg + z[19] * sg;

/* --- Term 5: angle = 4E - 8M + 3J - l --- */
g = STR * (a - MP);
cg = cos(g);
sg = sin(g);
l  += -0.165009 * cg + 0.044176 * sg;
l1 += 4.67 * cg + 201.55 * sg;
/* No l2 term */

/* --- Term 6: angle = 18V - 16E --- */
g = STR * f;
cg = cos(g);
sg = sin(g);
l  += 0.330401 * cg + 0.661362 * sg;
l1 += 1202.67 * cg - 555.59 * sg;
l2 += z[20] * cg + z[21] * sg;

/* --- Term 7: angle = 18V - 16E - 2l --- */
g = STR * (f - 2.0*MP);
cg = cos(g);
sg = sin(g);
l  += 0.352185 * cg + 0.705041 * sg;
l1 += 1283.59 * cg - 586.43 * sg;
/* No l2 term in DE404 (MOSH_MOON_200 has z[48]/z[49] here — not present in DE404!) */

/* --- Term 8: angle = 2J - 5S --- */
g = STR * (2.0*Ju - 5.0*Sa);
cg = cos(g);
sg = sin(g);
l  += -0.034700 * cg + 0.160041 * sg;
l2 += z[22] * cg + z[23] * sg;
/* No l1 term */

/* --- Term 9: angle = L - F = SWELP - NF --- */
g = STR * (SWELP - NF);
cg = cos(g);
sg = sin(g);
l  += 0.000116 * cg + 7.063040 * sg;
l1 +=  298.8 * sg;   // only sin term; no cos term for l1
/* No l2 term */
```

#### T³ radius and longitude corrections

```c
/* T^3 terms — note: DE404 has far fewer than MOSH_MOON_200 */

sg = sin( STR * M );
/* CRITICAL: this is ASSIGNMENT (=), not accumulation (+=)!
   Original Moshier code had +=, which was a bug (l3 uninitialized).
   Fixed in DE404 branch: */
l3 = z[24] * sg;   // z[24] = -1.999840061168e+00 (t^3 sin(M))
l4 = 0;            // explicitly zeroed; no t^4 planetary terms in DE404

/* Radius corrections scaled by T (so they end up as T^3 relative to the T^2 base) */
g = STR * (2.0*D - M);
sg = sin(g);
cg = cos(g);
moonpol[2] +=  -0.2655 * cg * T;

g = STR * (M - MP);
moonpol[2] +=  -0.1568 * cos(g) * T;

g = STR * (M + MP);
moonpol[2] +=   0.1309 * cos(g) * T;

g = STR * (2.0*(D + M) - MP);
sg = sin(g);
cg = cos(g);
moonpol[2] +=   0.5568 * cg * T;

/* Accumulate chewm LRT2 longitude result into l2 */
l2 += moonpol[0];   // moonpol[0] now holds sum of LRT2 table longitude terms

g = STR * (2.0*D - M - MP);
moonpol[2] +=  -0.1910 * cos(g) * T;

/* Scale: multiply accumulated T^2 terms by T so they become T^3 relative to caller */
moonpol[1] *= T;   // latitude  × T
moonpol[2] *= T;   // radius    × T
```

### Phase B: T¹ Series (chewm + planetary perturbations)

```c
/* Reset longitude accumulator — moonpol[1], [2] still hold T^2×T scaled values */
moonpol[0] = 0.0;

chewm( BT,  NBT,  4, 4, moonpol );  // 16 latitude rows  → moonpol[1]
chewm( LRT, NLRT, 4, 1, moonpol );  // 41 longitude+radius rows → moonpol[0], moonpol[2]
```

#### Latitude perturbations — Phase B (units: arcsec×10^-5)

```c
g = STR * (f - MP - NF - 2355767.6);  /* 18V - 16E - l - F */
moonpol[1] +=  -1127. * sin(g);

g = STR * (f - MP + NF - 235353.6);   /* 18V - 16E - l + F */
moonpol[1] +=  -1123. * sin(g);

g = STR * (Ea + D + 51987.6);
moonpol[1] +=   1303. * sin(g);

g = STR * SWELP;
moonpol[1] +=    342. * sin(g);
```

#### Longitude+speed perturbations — Phase B

```c
/* --- Term 10: angle = 2V - 3E --- */
g = STR * (2.*Ve - 3.*Ea);
cg = cos(g);
sg = sin(g);
l  +=  -0.343550 * cg - 0.000276 * sg;
l1 +=   105.90  * cg + 336.53 * sg;

/* --- Term 11: angle = 18V - 16E - 2D --- */
g = STR * (f - 2.*D);
cg = cos(g);
sg = sin(g);
l  += 0.074668 * cg + 0.149501 * sg;
l1 += 271.77   * cg - 124.20  * sg;

/* --- Term 12: angle = 18V - 16E - 2D - l --- */
g = STR * (f - 2.*D - MP);
cg = cos(g);
sg = sin(g);
l  += 0.073444 * cg + 0.147094 * sg;
l1 += 265.24   * cg - 121.16  * sg;

/* --- Term 13: angle = 18V - 16E + 2D - l --- */
g = STR * (f + 2.*D - MP);
cg = cos(g);
sg = sin(g);
l  += 0.072844 * cg + 0.145829 * sg;
l1 += 265.18   * cg - 121.29  * sg;

/* --- Term 14: angle = 18V - 16E + 2(D - l) --- */
g = STR * (f + 2.*(D - MP));
cg = cos(g);
sg = sin(g);
l  += 0.070201 * cg + 0.140542 * sg;
l1 += 255.36   * cg - 116.79  * sg;

/* --- Term 15: angle = E + D - F --- */
g = STR * (Ea + D - NF);
cg = cos(g);
sg = sin(g);
l  +=  0.288209 * cg - 0.025901 * sg;
l1 += -63.51   * cg - 240.14  * sg;

/* --- Term 16: angle = 2E - 3J + 2D - l --- */
g = STR * (2.*Ea - 3.*Ju + 2.*D - MP);
cg = cos(g);
sg = sin(g);
l  += 0.077865 * cg + 0.438460 * sg;
l1 += 210.57   * cg + 124.84  * sg;

/* --- Term 17: angle = E - 2M (Mars) --- */
g = STR * (Ea - 2.*Ma);
cg = cos(g);
sg = sin(g);
l  += -0.216579 * cg + 0.241702 * sg;
l1 +=  197.67  * cg + 125.23  * sg;

/* a = 4E - 8M + 3J is still in scope from Phase A */

/* --- Term 18: angle = 4E - 8M + 3J + l --- */
g = STR * (a + MP);
cg = cos(g);
sg = sin(g);
l  += -0.165009 * cg + 0.044176 * sg;
l1 +=   4.67   * cg + 201.55  * sg;

/* --- Term 19: angle = 4E - 8M + 3J + 2D - l --- */
g = STR * (a + 2.*D - MP);
cg = cos(g);
sg = sin(g);
l  += -0.133533 * cg + 0.041116 * sg;
l1 +=   6.95   * cg + 187.07  * sg;

/* --- Term 20: angle = 4E - 8M + 3J - 2D + l --- */
g = STR * (a - 2.*D + MP);
cg = cos(g);
sg = sin(g);
l  += -0.133430 * cg + 0.041079 * sg;
l1 +=   6.28   * cg + 169.08  * sg;

/* --- Term 21: angle = 3V - 4E --- */
g = STR * (3.*Ve - 4.*Ea);
cg = cos(g);
sg = sin(g);
l  += -0.175074 * cg + 0.003035 * sg;
l1 +=   49.17  * cg + 150.57  * sg;

/* --- Term 22: angle = 2(E + D - l) - 3J + 213534" --- */
/* (only l1, no l term) */
g = STR * (2.*(Ea + D - MP) - 3.*Ju + 213534.);
l1 +=  158.4 * sin(g);

/* Accumulate chewm LRT longitude result into l1 */
l1 += moonpol[0];   // moonpol[0] holds sum of LRT table longitude terms

/* Scale Phase B latitude and radius by 0.1*T (so final unit is 10^-4 arcsec) */
a = 0.1 * T;
moonpol[1] *= a;   // latitude
moonpol[2] *= a;   // radius
```

### End of `moon1()`

After `moon1()` returns:
- `l`, `l1`, `l2`, `l3`, `l4` hold polynomial coefficients for the Horner series
- `moonpol[1]` = latitude T-series contribution (combined T^2 and T^1 terms, scaled)
- `moonpol[2]` = radius T-series contribution (combined, scaled)
- `f = 18*Ve - 16*Ea` remains live in thread-local storage for use by `moon2()`

---

## `moon2()` — DE404 Version (lines 1367–1442)

T^0 planetary perturbation terms. Adds directly to `l` (longitude, arcsec×10^-5) and `B` (latitude, arcsec).

**Important**: `f = 18*Ve - 16*Ea` from `moon1()` is used in two of the latitude terms below.

### Longitude terms (28 terms, arcsec×10^-5)

```c
static void moon2(void) {
/* terms in T^0 */

g = STR * (2*(Ea - Ju + D) - MP + 648431.172);
l += 1.14307 * sin(g);

g = STR * (Ve - Ea + 648035.568);
l += 0.82155 * sin(g);

g = STR * (3*(Ve - Ea) + 2*D - MP + 647933.184);
l += 0.64371 * sin(g);

g = STR * (Ea - Ju + 4424.04);
l += 0.63880 * sin(g);

g = STR * (SWELP + MP - NF + 4.68);
l += 0.49331 * sin(g);

g = STR * (SWELP - MP - NF + 4.68);
l += 0.4914 * sin(g);

g = STR * (SWELP + NF + 2.52);
l += 0.36061 * sin(g);

g = STR * (2.*Ve - 2.*Ea + 736.2);
l += 0.30154 * sin(g);

g = STR * (2.*Ea - 3.*Ju + 2.*D - 2.*MP + 36138.2);
l += 0.28282 * sin(g);

g = STR * (2.*Ea - 2.*Ju + 2.*D - 2.*MP + 311.0);
l += 0.24516 * sin(g);

g = STR * (Ea - Ju - 2.*D + MP + 6275.88);
l += 0.21117 * sin(g);

g = STR * (2.*(Ea - Ma) - 846.36);
l += 0.19444 * sin(g);

g = STR * (2.*(Ea - Ju) + 1569.96);
l -= 0.18457 * sin(g);   /* NOTE: subtraction (-=) */

g = STR * (2.*(Ea - Ju) - MP - 55.8);
l += 0.18256 * sin(g);

g = STR * (Ea - Ju - 2.*D + 6490.08);
l += 0.16499 * sin(g);

g = STR * (Ea - 2.*Ju - 212378.4);
l += 0.16427 * sin(g);

g = STR * (2.*(Ve - Ea - D) + MP + 1122.48);
l += 0.16088 * sin(g);

g = STR * (Ve - Ea - MP + 32.04);
l -= 0.15350 * sin(g);   /* NOTE: subtraction (-=) */

g = STR * (Ea - Ju - MP + 4488.88);
l += 0.14346 * sin(g);

g = STR * (2.*(Ve - Ea + D) - MP - 8.64);
l += 0.13594 * sin(g);

g = STR * (2.*(Ve - Ea - D) + 1319.76);
l += 0.13432 * sin(g);

g = STR * (Ve - Ea - 2.*D + MP - 56.16);
l -= 0.13122 * sin(g);   /* NOTE: subtraction (-=) */

g = STR * (Ve - Ea + MP + 54.36);
l -= 0.12722 * sin(g);   /* NOTE: subtraction (-=) */

g = STR * (3.*(Ve - Ea) - MP + 433.8);
l += 0.12539 * sin(g);

g = STR * (Ea - Ju + MP + 4002.12);
l += 0.10994 * sin(g);

g = STR * (20.*Ve - 21.*Ea - 2.*D + MP - 317511.72);
l += 0.10652 * sin(g);

g = STR * (26.*Ve - 29.*Ea - MP + 270002.52);
l += 0.10490 * sin(g);

g = STR * (3.*Ve - 4.*Ea + D - MP - 322765.56);
l += 0.10386 * sin(g);
```

### Latitude terms (8 terms, arcsec)

**Note**: `B` is ASSIGNED (=) on the first term, then accumulated (+=). Also uses `f` from `moon1()`.

```c
g = STR * (SWELP + 648002.556);
B =  8.04508 * sin(g);   /* ASSIGNMENT — B is initialized here */

g = STR * (Ea + D + 996048.252);
B += 1.51021 * sin(g);

/* f = 18*Ve - 16*Ea from moon1() — still live in TLS */
g = STR * (f - MP + NF + 95554.332);
B += 0.63037 * sin(g);

g = STR * (f - MP - NF + 95553.792);
B += 0.63014 * sin(g);

g = STR * (SWELP - MP + 2.9);
B +=  0.45587 * sin(g);

g = STR * (SWELP + MP + 2.5);
B += -0.41573 * sin(g);   /* NOTE: negative amplitude */

g = STR * (SWELP - 2.0*NF + 3.2);
B +=  0.32623 * sin(g);

g = STR * (SWELP - 2.0*D + 2.5);
B +=  0.29855 * sin(g);
}
```

---

## `mean_elements()` — DE404 Version (lines 1763–1818)

```c
static void mean_elements(void) {
double fracT = fmod(T, 1);

/* Mean anomaly of sun = l' (J. Laskar)
   Split into integer + fractional multipliers to avoid precision loss: */
M = mods3600(129600000.0 * fracT - 3418.961646 * T + 1287104.76154);
M += ((((((((
   1.62e-20    * T
 - 1.0390e-17  ) * T
 - 3.83508e-15 ) * T
 + 4.237343e-13) * T
 + 8.8555011e-11) * T
 - 4.77258489e-8 ) * T
 - 1.1297037031e-5) * T
 + 1.4732069041e-4 ) * T
 - 0.552891801772  ) * T2;

/* DE404 branch (else of #ifdef MOSH_MOON_200): */

/* Mean distance of moon from its ascending node = F */
/* Effective: (1739527263.0983 - 2.079419901760e-01) * T + 335779.55755 */
NF = mods3600(1739232000.0 * fracT + 295263.0983 * T - 2.079419901760e-01 * T + 335779.55755);

/* Mean anomaly of moon = l */
/* Effective: (1717915923.4728 - 2.035946368532e-01) * T + 485868.28096 */
MP = mods3600(1717200000.0 * fracT + 715923.4728 * T - 2.035946368532e-01 * T + 485868.28096);

/* Mean elongation of moon = D */
/* Effective: (1602961601.4603 + 3.962893294503e-01) * T + 1072260.73512 */
D = mods3600(1601856000.0 * fracT + 1105601.4603 * T + 3.962893294503e-01 * T + 1072260.73512);

/* Mean longitude of moon, referred to mean ecliptic and equinox of date */
/* Effective: (1732564372.83264 - 6.784914260953e-01) * T + 785939.95571 */
SWELP = mods3600(1731456000.0 * fracT + 1108372.83264 * T - 6.784914260953e-01 * T + 785939.95571);

/* Higher-degree secular correction using DE404 z[] (quadratic in T, times T^2):
   Each of the 4 elements gets a 3-term correction: z[i]*T^2 + z[i+1]*T^3 + z[i+2]*T^4 */
NF    += ((z[2]*T + z[1])*T + z[0])*T2;   /* F:  z[0..2] */
MP    += ((z[5]*T + z[4])*T + z[3])*T2;   /* l:  z[3..5] */
D     += ((z[8]*T + z[7])*T + z[6])*T2;   /* D:  z[6..8] */
SWELP += ((z[11]*T + z[10])*T + z[9])*T2; /* L:  z[9..11] */
}
```

**Unit**: All of `M`, `NF`, `MP`, `D`, `SWELP` are in arcseconds after `mods3600` (range ≈ ±648000).

**`mods3600`** (line 1734): reduces x modulo 360 degrees expressed in arcseconds (i.e., modulo 1296000"):
```c
static double mods3600(double x) {
    x = fmod(x, 1296000.0);
    if (x < 0.0) x += 1296000.0;
    return x;
}
```

---

## `mean_elements_pl()` — Planetary Mean Longitudes (lines 1820–1850)

All elements in arcseconds, evaluated as polynomials in T with `mods3600` reduction.

```c
void mean_elements_pl(void) {
/* Mean longitudes of planets (Laskar, Bretagnon) */

Ve = mods3600( 210664136.4335482 * T + 655127.283046 );
Ve += ((((((((
  -9.36e-023 * T
 - 1.95e-20  ) * T
 + 6.097e-18 ) * T
 + 4.43201e-15) * T
 + 2.509418e-13) * T
 - 3.0622898e-10) * T
 - 2.26602516e-9 ) * T
 - 1.4244812531e-5) * T
 + 0.005871373088 ) * T2;

Ea = mods3600( 129597742.26669231 * T + 361679.214649 );
Ea += (((((((( -1.16e-22 * T
 + 2.976e-19 ) * T
 + 2.8460e-17 ) * T
 - 1.08402e-14 ) * T
 - 1.226182e-12 ) * T
 + 1.7228268e-10 ) * T
 + 1.515912254e-7 ) * T
 + 8.863982531e-6 ) * T
 - 2.0199859001e-2 ) * T2;

Ma = mods3600( 68905077.59284 * T + 1279559.78866 );
Ma += (-1.043e-5*T + 9.38012e-3)*T2;

Ju = mods3600( 10925660.428608 * T + 123665.342120 );
Ju += (1.543273e-5*T - 3.06037836351e-1)*T2;

Sa = mods3600( 4399609.65932 * T + 180278.89694 );
Sa += (( 4.475946e-8*T - 6.874806E-5 ) * T + 7.56161437443E-1)*T2;
}
```

---

## Key Differences: DE404 vs MOSH_MOON_200

| Feature | MOSH_MOON_200 | DE404 (used by Swiss Ephemeris) |
|---------|--------------|--------------------------------|
| `z[]` size | 71 elements | 25 elements |
| Index mapping | z[0..5]=F,l,D,L each 6 terms | z[0..11]=F,l,D,L each 3 terms |
| `l3` initialization | `l3 += ...` multiple terms | `l3 = z[24]*sin(M)` (bug fix) |
| `l4` | multiple terms | `l4 = 0` |
| Term 7 (18V-16E-2l) l2 | z[48]/z[49] | none |
| T^3 l3 series (2D-M, etc.) | present | absent |
| ss/cc zeroing | not done | explicit double-loop |
| `mean_elements` secular | 6-term polynomial per element | 3-term polynomial per element |

---

## Amplitude and Scale Notes

### `l` accumulator family

`l`, `l1`, `l2`, `l3`, `l4` are all `double` fields with implicit scales:

- `l` accumulates in **arcsec** (the T^0 planetary longitude corrections, e.g., `6.367278 * cg + 12.747036 * sg` gives ~±14 arcsec max).
- `l1`, `l2`, `l3`, `l4` are polynomial coefficients. The Horner in `moon3()` is:
  ```
  l += (((l4·T + l3)·T + l2)·T + l1)·T · 1e-5
  ```
  So `l1*T*1e-5` is added to `l` in arcsec. For `l1 = 23123.70`, that gives 0.231·T arcsec/century — a slow secular drift.

- The DE404 `z[]` comment says these longitude terms are "arc seconds times 10^5". This means the values (e.g., z[12] = -84.3) represent coefficients in the Horner *before* the 1e-5 factor; the physical T^2 contribution is z[12]*cos(angle)*T^2*1e-5 arcsec. For T=1 century, ≈ -0.000843 arcsec.

**Do not second-guess the scales.** Copy coefficients exactly and apply the Horner with the 1e-5 factor as written.

### `moonpol[]` accumulator

| Stage | moonpol[0] (lon) | moonpol[1] (lat) | moonpol[2] (rad) |
|-------|-----------------|-----------------|-----------------|
| After chewm T^2 calls | arcsec×10^-5 | arcsec×10^-5 | km×10^-5 |
| After `l2 += moonpol[0]` | — (l2 now holds this) | — | — |
| After `moonpol[1] *= T; moonpol[2] *= T` | — | ×T | ×T |
| After `moonpol[0] = 0.0` (Phase B start) | reset | scaled | scaled |
| After chewm T^1 calls | arcsec×10^-5 (new) | accumulated | accumulated |
| After `l1 += moonpol[0]` | — (l1 now holds this) | — | — |
| After `a = 0.1*T; moonpol[1] *= a; moonpol[2] *= a` | — | ×0.1T more | ×0.1T more |
| In moon3: `moonpol[0]` after chewm LR/MB | arcsec×10^-5 | arcsec×10^-5 | km×10^-5 |
| Final: `moonpol[1] = 1e-4*moonpol[1] + B` | — | arcsec | — |
| Final: `moonpol[2] = 1e-4*moonpol[2] + 385000.52899` | — | — | km |

The final longitude: `moonpol[0] = SWELP + l + 1e-4 * moonpol[0]` in arcsec.
