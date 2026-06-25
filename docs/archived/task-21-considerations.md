## Developer Blindspots & Considerations: Nutation Algorithms

### 1. The Cache Temptation & State Isolation

* **The C Interpolation Trap:** In `swephlib.c`, `swi_nutation` implements a 3-point quadratic interpolation cache using static arrays (`t_cache`, `dpsi_cache`, `deps_cache`) to avoid re-evaluating the heavy IAU 2000 series for nearby time steps. Because the task explicitly forbids this interpolation cache, Claude must be watched closely.
* *Claude Watchout:* Claude may look at the performance regression of dropping the cache and attempt to sneak inner mutability or a hidden global state manager into the `Ephemeris` structure or module scope.
* *Correction:* Enforce that `swi_nutation` is a completely stateless, pure function matching the signature `fn nutation(tjd: f64, model: NutationModel) -> Nutation`. Performance optimizations must be pushed downstream to user-space architecture.



### 2. Disambiguating Fundamental Arguments (Delaunay Expressions)

* **IAU 1980 vs. IAU 2000 Coefficients:** A critical pitfall in the C code is that the 5 Delaunay arguments ($l, l', F, D, \Omega$) do not use the same polynomial coefficients across all models. IAU 1980 uses the expressions from Seidelmann (1982), while IAU 2000/2006 uses the IERS 1996/2003 expressions, which include higher-degree terms and different constants.
* *Claude Watchout:* Claude will likely try to create a single, unified `calc_delaunay_arguments(t: f64)` function to deduplicate code. This will silently corrupt the numerical fidelity of one of the families.
* *Correction:* Force Claude to explicitly decouple or parameterize the fundamental argument generation based on the target astronomical standard (e.g., `FundamentalArgs::compute_iau1980(t)` vs. `FundamentalArgs::compute_iau2000(t)`).



### 3. Unit Scale and Precision Inversions

* **Arcseconds, Milliarcseconds, and Radians:** The underlying C implementation stores nutation table coefficients as integers or compact floats scaled to specific units (e.g., ten-thousandths of an arcsecond or milliarcseconds) to prevent precision loss during C structural packing. The final output must be in radians.
* *Claude Watchout:* Claude frequently mixes up conversion factors when porting raw values from C arrays. For example, it might apply a global `arcsec_to_rad` multiplier to coefficients that are already scaled down by an order of magnitude, or apply it prematurely before the time-dependent polynomial multiplication ($T \times \text{coeff}$).
* *Correction:* All raw data tables must explicitly document their base units in Rust code comments, and the conversion to radians must happen uniformly at the very end of the summation loop:

$$\Delta\psi_{\text{rad}} = \Delta\psi_{\text{units}} \times \text{SCALE\_FACTOR}$$





### 4. Memory Footprint and Compilation Bottlenecks

* **Large Table Compilation Stack:** The IAU 2000A model contains 678 Luni-solar terms and 687 Planetary terms. Storing these as unoptimized, deeply nested structural literals in code blocks can cause Rust compiler memory usage to skyrocket or blow past stack limits during evaluation if handled incorrectly.
* *Claude Watchout:* Claude might model these rows using overly complex structural layouts or, worse, attempt to load them via safe initialization routines at runtime.
* *Correction:* Tables must be hardcoded as flat, contiguous `const` arrays of primitive types or simple, trivial `Copy` structs (e.g., `const IAU2000A_LUNI_SOLAR: [NutationTerm; 678] = [...];`).



### 5. Loop Mechanics & Safe Iteration

* **Pointer Arithmetic Translation:** The C implementation steps through these massive tables using pointer arithmetic offsets (`struct nut *np`).
* *Claude Watchout:* Claude might introduce explicit index trackers (`table[i]`) and trigger bound-checking penalties on every single term iteration, or attempt to use `unsafe` blocks to match the C pointer iteration.
* *Correction:* Enforce clean, idiomatic, and bounds-check-optimized iterators. The compiler will auto-vectorize loops structured like this:
```rust
for term in IAU2000A_LUNI_SOLAR.iter() {
    // Accumulate pure math step
}

```





### 6. Derived Outputs Optimization

* **The Matrix/Sine Dilemma:** The C struct provides `snut` ($\sin(\Delta\psi)$) and a full nutation matrix.
* *Correction:* Do not allow Claude to blow up the output type of the core math function. The pure function should exclusively return `Nutation { dpsi: f64, deps: f64 }`. Provide a separate helper implementation block or utility function on that struct to allow downstream callers to compute the sine value or the transformation matrix from those raw angles if needed.



Are you ready to review Claude's architectural plan for compiling these large coefficient arrays and loop structures?
