Porting `swephlib.c` is where we lay the raw mathematical tracks for the entire ephemeris engine. Because these functions are called millions of times during a dense chart or ephemeris generation run, any subtle drift in precision or performance overhead here will cascade heavily into later modules.

Here are the critical traps, optimizations, and structural details to look out for when drafting the final implementation plan.

---

## 1. Ergonomics vs. Type Safety for Vectors

The plan suggests using raw arrays like `[f64; 3]` and `[f64; 6]`. While this maps nicely to C pointers, it leaves room for subtle bugs—like accidentally feeding a position vector into a velocity calculation, or mixing up a 6-element state vector's components.

* **Watch out for:** Using unstructured arrays loses the semantic clarity of what the data represents.
* **Recommendation:** Introduce thin, zero-cost semantic wrapper types or aliases early in the plan:
```rust
pub struct Vector3D(pub [f64; 3]);
pub struct StateVector(pub [f64; 6]); // or struct { pos: Vector3D, vel: Vector3D }

```


This keeps the memory layout identical to C (`#[repr(transparent)]` or simple naming conventions) while forcing compile-time correctness.

## 2. The Polar Singularity Trap ($r_{xy} \to 0$)

In `cartesian_to_polar_with_velocity`, the Jacobian matrix transformation involves dividing by $r_{xy} = \sqrt{x^2 + y^2}$. When a body is exactly at or incredibly close to the geographic/ecliptic pole, $r_{xy}$ approaches zero.

* **Watch out for:** C code often implicitly relies on platform-specific floating-point behavior or explicit `if (rxy < 1e-20)` guards to prevent a `NaN` explosion during velocity projection.
* **Recommendation:** Ensure the plan mandates a strict audit of how the original C handles division when $r_{xy} \approx 0$, and explicitly test the Rust translation with coordinates like `[0.0, 0.0, 1.0]` to ensure velocities don't pan out to `NaN` or panic on division.

## 3. Do Not Skip the Centisecond (`csnorm`) Functions

The preliminary plan marks centisecond variants as *"maybe skip or make optional"*. **This is a dangerous trap.**

* **Watch out for:** The core Swiss Ephemeris engine (especially code handling house systems, aspect grids, and fast eclipses) heavily uses 32-bit integer centiseconds (`int32` or `centisec`) for rapid, exact angular comparisons and internal caching.
* **Recommendation:** Do not skip them. They should be ported alongside their degree/radian counterparts. They are vital for avoiding cumulative float rounding errors in downstream modules.

## 4. `f64::rem_euclid` vs. Literal C Expressions

For angle normalization (`normalize_degrees`, `normalize_radians`), Rust developers naturally reach for `x.rem_euclid(360.0)`.

* **Watch out for:** The original C code typically utilizes explicit flooring expressions, such as:
```c
x -= floor(x / 360.0) * 360.0;

```


While mathematically identical to `rem_euclid`, floating-point operations at extreme values (e.g., coordinates input at $10^{12}$ degrees or values exactly on the $360.0$ boundary) can yield sub-bit discrepancies due to compiler-specific FMA (Fused Multiply-Add) instructions.
* **Recommendation:** To honor the bitwise golden-data mandate from `CLAUDE.md`, the plan should dictate writing the Rust math exactly as written in C first, rather than relying on idiomatic Rust abstraction shortcuts, unless differential tests prove `rem_euclid` is identical down to the last bit.

## 5. Chebyshev Slice Invariants and Clenshaw Performance

The Chebyshev evaluation routines (`swi_echeb` and `swi_edcheb`) rely on a slice of coefficients (`coeffs: &[f64]`).

* **Watch out for:** In a tight loop, indexing into a slice (`coeffs[i]`) introduces Rust bounds checking overhead. Furthermore, the recurrence relation requires safe access up to the polynomial degree specification.
* **Recommendation:** The plan should ensure that the recurrence loops utilize iterators, or that the function establishes a clear precondition check on the slice length at the entry point so the compiler can safely optimize away interior bounds checks.

## 6. Strict Semantic Mapping of Float-to-Integer Rounding

The function `swe_d2l(x)` converts a `double` to a long integer with rounding.

* **Watch out for:** Rust's standard `as i32` cast **truncates** toward zero. C's rounding routines in `swephlib.c` often implement a specific banker's rounding, round-half-away-from-zero, or custom offset adjustments to ensure consistency across negative and positive coordinates.
* **Recommendation:** Explicitly document the exact rounding direction expected by `swe_d2l` (e.g., using `.round()`, `.floor()`, or `.ceil()` before casting) so historical BCE charts don't wind up offset by a day or a minute of arc due to truncation errors.

---

Should we structure the vector-handling math to explicitly use standard nested tuples/arrays, or would you like to introduce a lightweight geometry/vector utility internal to the math module to make the coordinate transformations more readable?
