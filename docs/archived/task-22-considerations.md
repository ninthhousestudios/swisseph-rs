## Developer Blindspots & Considerations: Delta-T ($\Delta T$) Computation

### 1. Array Boundary Vulnerabilities in Bessel Interpolation

* **The 4-Point Stencil Risk:** The 4th-order Bessel interpolation requires a stencil of four consecutive points from the tabulated data: $y_{-1}, y_0, y_1, y_2$. For a given fractional year $Y$, if $Y$ falls within the very first interval ($1620 \le Y < 1621$) or the very last interval of the table (e.g., $Y_{\text{max}-1} \le Y < Y_{\text{max}}$), a naive array indexing implementation will cause a runtime panic (underflow on index `-1` or out-of-bounds on index `+2`).
* *Claude Watchout:* Claude will likely write a generic interpolation helper assuming the stencil is always valid, neglecting the edge fallback logic present in `swephlib.c`.
* *Correction:* Force Claude to explicitly handle boundary mitigation. For the first and last slots, the C implementation drops down to a lower-order (linear or quadratic) interpolation or alters the stencil window selection. Ensure the code uses safe slicing and explicit fallback paths without out-of-bounds risks.



### 2. Time Parameter Inconsistencies ($TJD$ vs. Decimal Year)

* **The Interconversion Pitfall:** The entry point takes `tjd` (Julian Day), but the router and the individual model polynomials switch between Julian centuries ($T$), Millenniums ($M$), and decimal calendar years ($Y$).
* *Claude Watchout:* Claude might use inconsistent conversions for the decimal year across different models. For instance, a naive `Y = 2000.0 + (tjd - 2451545.0) / 365.25` introduces phase errors compared to the exact calendar date conversion routines utilized in the original C code for table indexing.
* *Correction:* The plan must specify a unified, deterministic internal helper to compute the exact decimal year fraction from `tjd` that matches the C library's assumptions for table lookups versus polynomial evaluations.



### 3. State Infiltration via Ephemeris Configuration

* **Tidal Acceleration Defaults:** The tidal acceleration adjustment formula modifies the output based on a reference acceleration ($\dot{n}$):

$$\text{correction} = -0.000091 \times (\text{tid\_acc} - \text{tid\_acc\_ref}) \times (\text{Year} - 1955)^2$$



In C, `tid_acc_ref` changes implicitly depending on which ephemeris file backend is currently compiled or active in global state.
* *Claude Watchout:* Claude might hardcode a single reference value or introduce a match block on a global enum, breaking isolation.
* *Correction:* This calculation must remain completely stateless. The `tid_acc` and its corresponding `tid_acc_ref` must be derived cleanly from the passed `&EphemerisConfig` struct. If the config dictates a default behavior, it must map explicitly to the selected ephemeris type variant encapsulated in that configuration object.



### 4. Verification of Leap Seconds vs. UT1

* **Semantic Disambiguation:** As noted, $\Delta T = TT - UT1$. While UTC depends on discrete atomic leap seconds, $UT1$ is driven by the Earth's variable rotation.
* *Claude Watchout:* Claude may mistakenly cross-reference or integrate the `Ephemeris.leap_seconds()` table when performing extrapolation for the near future, blending UTC step-adjustments into a continuous $UT1$ trend curve.
* *Correction:* Verify that Claude treats the tabulated data array up to the present day as the absolute truth source for $UT1$ variations. Leap second steps should play absolutely no role in this mathematical pipeline; future extrapolation must rely purely on the specified long-term linear/parabolic drift models.



### 5. Table Archiving Strategy

* **Compact Storage:** The $\Delta T$ table spans ~400 entries from 1620 onwards.
* *Correction:* Ensure Claude exposes this array as a `static` or `const` slice of `f64` values packed directly inside the `deltat` internal module. Do not allow it to implement dynamic file reading or runtime heap allocations for this static dataset.



---

How does Claude plan to handle the boundary conditions of the 4th-order Bessel interpolation for dates immediately adjacent to the 1620 and ~2024 table thresholds?
