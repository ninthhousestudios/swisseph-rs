use std::path::Path;
use std::process::Command;

#[test]
fn c_smoke_test() {
    let crate_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let c_src = crate_dir.join("tests/c/smoke.c");
    let include_dir = crate_dir.join("include");

    let profile = if cfg!(debug_assertions) {
        "debug"
    } else {
        "release"
    };
    let target_dir = crate_dir.parent().unwrap().join("target").join(profile);

    let out_dir = std::env::var("OUT_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| target_dir.clone());
    let exe_path = out_dir.join("c_smoke_test");

    let cc = std::env::var("CC").unwrap_or_else(|_| "cc".into());

    let lib_name = if cfg!(target_os = "macos") {
        "swisseph_ffi"
    } else {
        "swisseph_ffi"
    };

    let compile = Command::new(&cc)
        .arg(&c_src)
        .arg(format!("-I{}", include_dir.display()))
        .arg(format!("-L{}", target_dir.display()))
        .arg(format!("-l{lib_name}"))
        .arg("-lm")
        .arg("-lpthread")
        .arg("-ldl")
        .arg("-o")
        .arg(&exe_path)
        .output();

    let compile = match compile {
        Ok(o) => o,
        Err(e) => {
            eprintln!("Skipping C smoke test: C compiler `{cc}` not available: {e}");
            return;
        }
    };

    if !compile.status.success() {
        let stderr = String::from_utf8_lossy(&compile.stderr);
        let stdout = String::from_utf8_lossy(&compile.stdout);
        panic!(
            "C compilation failed (exit {}):\nstdout: {stdout}\nstderr: {stderr}",
            compile.status
        );
    }

    let mut run_cmd = Command::new(&exe_path);
    run_cmd.env(
        "LD_LIBRARY_PATH",
        format!(
            "{}:{}",
            target_dir.display(),
            std::env::var("LD_LIBRARY_PATH").unwrap_or_default()
        ),
    );
    if cfg!(target_os = "macos") {
        run_cmd.env(
            "DYLD_LIBRARY_PATH",
            format!(
                "{}:{}",
                target_dir.display(),
                std::env::var("DYLD_LIBRARY_PATH").unwrap_or_default()
            ),
        );
    }

    let run = run_cmd.output().expect("failed to run C smoke test binary");
    let stdout = String::from_utf8_lossy(&run.stdout);
    let stderr = String::from_utf8_lossy(&run.stderr);
    println!("{stdout}");
    if !stderr.is_empty() {
        eprintln!("{stderr}");
    }
    assert!(
        run.status.success(),
        "C smoke test exited with {}",
        run.status
    );
}
