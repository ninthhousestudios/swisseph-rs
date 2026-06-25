## Developer Blindspots & Considerations: Precession Router & Algorithms

### 1. Hidden C State & Caching Traps (Critical)

* **The Static Cache Illusion:** In `swephlib.c`, `swi_precess` and its underlying algorithms heavily rely on internal `static` variables to cache the last evaluated time (`tjd`), model, and calculated rotation matrices. This allows successive calls for the same date to skip computation.
* *Claude Watchout:* Claude will be tempted to introduce `lazy_static`, `thread_local`, or inner mutability (`RefCell`/`Mutex`) to preserve this performance optimization. **This violates Invariant 1 (Stateless Design).** * *Correction:* The Rust implementation must be a pure, stateless function: `fn precess(x: &[f64; 3], tjd: f64, model: AstroModel, direction: Direction) -> [f64; 3]`. Optimization via caching must be entirely omitted from this layer; every call must compute from scratch.



### 2. Numerical Fidelity & Polynomial Evaluation

* **Horner's Method vs. Naive Expansion:** `precess_2` evaluates polynomials up to $T^{10}$ (11 terms for long-term models like Laskar 1986).
* *Claude Watchout:* Claude might naively use `f64::powf` or write out the polynomial expansion explicitly (e.g., `a + b*T + c*T*T...`). This introduces floating-point accumulation errors and destroys numerical fidelity over the required $\pm 5000$-year testing range.
* *Correction:* Force Claude to evaluate all polynomial arrays using an explicit, loop-based Horner’s method to guarantee bit-for-bit parity with the C implementation's precision characteristics.


* **Angle Reduction & Phase Drift:** `precess_3` (Vondrák 2011) utilizes a mix of 4th-degree polynomials and 24 periodic terms (sine/cosine arguments). At extreme dates ($\pm 5000$ years), large values of $T$ passed into trigonometric functions will cause phase drift if angle reduction is not identical to C.
* *Correction:* Ensure Claude does not use standard Rust `rem_euclid` or `f64::sin` blindly without verifying how `swephlib.c` handles argument reduction for periodic terms (check for custom wrapping constants or explicit modulo operations against $2\pi$).



### 3. API Patterns & Type Safety

* **Direction Semantics:** The C implementation passes `direction` as an integer (`1` for forward, `-1` or `0` for backward) and often uses it mathematically as a sign multiplier for angles, or to invert the rotation matrix.
* *Claude Watchout:* Claude might pass raw `i32` or `isize` for the direction parameter to preserve this math shortcut.
* *Correction:* Define a explicit enum `Direction { Forward, Backward }`. The conversion or sign application must be explicitly matched inside the function, preserving type safety without leaking raw magic numbers into the API.


* **Model Coverage:** The `AstroModels` enum in `types.rs` must exhaustively map `SEMOD_PREC_IAU_1976` through `SEMOD_PREC_NEWCOMB`. Ensure Claude implements `TryFrom<i32>` for this enum rather than using `unsafe` casting or raw integer matching in the router switch block.

### 4. Algorithmic Matrix Transposition Traps

* **Matrix Multiplication Order:** `precess_1` constructs Euler angle rotations $R_z(Z) \times R_y(-\Theta) \times R_z(z)$.
* *Claude Watchout:* When handling the inverse direction (backward precession), C often transposes the final matrix or reverses the sequence of angle applications. Claude is prone to mixing up row-major versus column-major ordering when converting C array indexing (e.g., `A[i][j]`) to Rust arrays.
* *Correction:* Ensure matrix multiplications are verified explicitly against the directional logic of C. Do not allow `unsafe` pointer offsets or raw structural casts for coordinate manipulation. Use structural `[f64; 3]` for vectors and clear, safe multi-dimensional arrays or dedicated math types for the $3 \times 3$ transformation matrices.



### 5. Strategy for Claude's Implementation Plan

Before Claude begins generating code, ensure its plan explicitly addresses:

1. The exact mathematical definition and execution order of Horner's method for polynomial evaluations.
2. An exhaustive mapping table for all 11 models within a clean `match` statement in the `swi_precess` router.
3. A definition of how the forward/inverse direction modifies the matrix operations without mutating external arguments.
