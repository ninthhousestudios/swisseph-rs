Here are the critical areas that require particular attention when Opus drafts the final implementation plan for the Date/Time module (`swedate.c`). These points ensure we strictly adhere to the architecture invariants in `CLAUDE.md` and avoid common pitfalls when porting legacy C astronomical math to Rust.

---

## 1. Stateless Design & File I/O for Leap Seconds

The preliminary plan mentions: *"loading from file can be a method on `Ephemeris`"*.

According to `CLAUDE.md`, `Ephemeris` must hold **only read-only configuration** and have **no mutable cache or internal state**.

* **Attention Item:** The final plan must explicitly define how `seleapsec.txt` is handled. It should be loaded and parsed during the construction of `Ephemeris` via `EphemerisConfig` (or provided as an optional path/buffer in the config).
* The functions inside the date module should accept the parsed leap-second data as a read-only reference argument, ensuring the pipeline remains entirely pure.

## 2. Strong Typing for Time Scales (TT vs. UT1 vs. UTC)

The plan proposes returning a struct `{ tt: f64, ut1: f64 }` for `utc_to_jd`, and notes considering a `JulianDay` newtype.

* **Attention Item:** In astronomy, mixing up Ephemeris Time (TT/TDB) and Universal Time (UT1/UTC) is a frequent source of severe bugs. Instead of a single generic `JulianDay` newtype, the actual plan should evaluate using **distinct newtypes** for different time scales (e.g., `JulianDayTt`, `JulianDayUt1`).
* If distinct newtypes are used, the plan needs to detail how ergonomics are preserved (e.g., implementing `Add<f64>`, `Sub`, or explicit conversion methods) so the mathematical code doesn't become overly verbose.

## 3. Numerical Fidelity: `floor` and Integer Truncation

The Julian Day formula involves several `floor` operations and mixed-type arithmetic:

```
if Gregorian: u2 = floor(|u|/100) − floor(|u|/400); JD -= u2 + 2

```

* **Attention Item:** C and Rust handle casting and rounding differently, especially with negative numbers. C's integer division truncates toward zero, whereas `floor` moves toward negative infinity.
* The plan must explicitly map out how Rust types (`f64`, `i32`) will interact. Every instance of C integer truncation—such as `(long)` casts in `swedate.c`—must be mapped to explicit Rust operations (like `.floor()` or `as i32`) rather than blind translation, as this directly impacts negative year (BCE) handling.

## 4. The Gregorian Gap and Historical Edge Cases

`swe_julday` handles the historical switch from the Julian to the Gregorian calendar (October 1582) via `gregflag`.

* **Attention Item:** The plan needs a definitive strategy for handling the "missing days" (October 5 to October 14, 1582). If an invalid calendar date within this gap is passed to `swe_date_conversion`, it must return a properly structured `Result::Err(Error::...)` rather than panicking or producing an undefined Julian Day.

## 5. Automated Golden Data Testing

`CLAUDE.md` states that numerical fidelity is non-negotiable.

* **Attention Item:** Relying solely on a few manual test cases (like J2000 or Unix epoch) is insufficient for a math-heavy module. The implementation plan should include a **differential testing harness**.
* This harness should ingest a generated CSV/JSON file containing thousands of inputs and outputs directly from the C `swetest` utility—covering extreme leap years, deep historical BCE dates, and every single leap second transition transition boundary—and assert bitwise or $epsilon$-exact equality in Rust.

---

Should we formalize the distinct newtypes for `JulianDayTt` and `JulianDayUt1` now, or would you prefer to keep them as a single `JulianDay` wrapper with explicit semantic variable names to start?
