## Developer Blindspots & Considerations: Sidereal Time

### 1. The Time Scale Cross-Contamination Trap (Critical)

* **The Conflict:** `swe_sidtime` takes a Universal Time parameter (`tjd_ut`). Greenwich Mean Sidereal Time (GMST) is calculated as a polynomial function directly of $T_{\text{UT1}}$ (derived from `tjd_ut`). However, the Equation of the Equinoxes (EoE) depends on Nutation ($\Delta\psi$) and Obliquity ($\varepsilon$), which are driven by gravitational dynamics and must be evaluated using Terrestrial Time ($T_{\text{TT}}$ / Ephemeris Time `tjd_et`).
* *Claude Watchout:* Claude will likely pass `tjd_ut` blindly into the internal calls for nutation and obliquity. This will cause an immediate precision failure because it skips the $\Delta T$ time-shift step ($tjd\_et = tjd\_ut + \Delta T$).
* *Correction:* Ensure Claude's plan explicitly routes `tjd_ut` through the newly built `DeltaT` module to generate the correct ephemeris time before triggering any underlying nutation or obliquity functions.

### 2. Unit Asymmetry and C-API Signature Drift

* **The Discrepancy:** In the original C library, `swe_sidtime0` expects its input parameters `eps` (obliquity) and `nut` (nutation in longitude) to be passed in **degrees**. However, internally, our pure Rust rewrite aims to normalize modern angular measurements to **radians** to eliminate constant runtime step conversions.
* *Claude Watchout:* Claude will easily mix up these conventions. If it replicates the exact code from `swe_sidtime0`, it will expect degrees, while the rest of our ported pipeline outputs radians—leading to a silent, severe scaling corruption of the math.
* *Correction:* Enforce complete internal unit homogeneity. The ported Rust core function must take `eps` and `nut` in radians. The final transformation of the accumulated polynomial terms must scale uniformly to hours at the very exit gate:

$$\text{Hours} = \text{Radians} \times \left(\frac{12.0}{\pi}\right)$$



### 3. Complementary Terms and Fundamental Arguments Linkage

* **The 33-Term Fourier Array:** The Equation of the Equinoxes for advanced models (IAU 2000+) uses a 33-term series of complementary terms to adjust for Earth's non-rigidity. This series requires evaluating Delaunay fundamental arguments.
* *Claude Watchout:* As highlighted in the nutation task, there are multiple historical expressions for Delaunay arguments (IAU 1980 vs. IAU 2000). Claude might pull the wrong Delaunay helper function or attempt to recycle arguments evaluated under an older model context.
* *Correction:* Ensure the 33 complementary terms are explicitly linked to the **IAU 2000 / IERS 2003** fundamental argument polynomial expressions, irrespective of the overarching model configuration chosen by the caller.

### 4. Floating-Point Range Reduction and Boundary Collisions

* **The 24-Hour Wrap:** Sidereal time must strictly fall within the $[0.0, 24.0)$ hour range. C often handles this via primitive loop-based accumulation: `while (st >= 24.0) st -= 24.0;`.
* *Claude Watchout:* Claude will look to optimize this using Rust's `rem_euclid(24.0)`. While mathematically sound, floating-point inaccuracies around numbers close to $24.0$ (e.g., $23.999999999999996$) can snap exactly to $24.0$ due to rounding behavior under certain hardware configurations.
* *Correction:* Direct Claude to implement an explicit boundary sanitization block right before returning the value: if the output equals or rounds up to `24.0`, it must explicitly collapse back to `0.0`.

### 5. Config Penetration for the Convenience Wrapper

* **The Signature Design:** The convenience function `swe_sidtime(tjd_ut)` cannot remain a single-parameter function in a stateless, pure architecture. It needs to know *which* default astronomical models are configured, and it needs access to the underlying ephemeris files or math rules.
* *Claude Watchout:* Claude might try to inject default fallback configurations inside the function body or pull them from an implicit global scope to preserve the single-argument signature of C.
* *Correction:* The ported function signature must be pure, explicit, and pass configuration context down:
```rust
pub fn sidereal_time(tjd_ut: f64, config: &EphemerisConfig) -> f64

```



---

Are you ready to audit Claude's structural definitions and polynomial data arrays for the 4 GMST models?
