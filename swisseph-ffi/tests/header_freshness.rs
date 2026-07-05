use std::path::Path;

#[test]
fn committed_header_matches_generated() {
    let crate_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let config = cbindgen::Config::from_file(crate_dir.join("cbindgen.toml"))
        .expect("failed to read cbindgen.toml");
    let generated = cbindgen::Builder::new()
        .with_crate(crate_dir)
        .with_config(config)
        .generate()
        .expect("cbindgen generation failed");

    let mut buf = Vec::new();
    generated.write(&mut buf);
    let generated_str = String::from_utf8(buf).expect("generated header not UTF-8");

    let committed = std::fs::read_to_string(crate_dir.join("include/swisseph.h"))
        .expect("include/swisseph.h not found — run cbindgen to generate it");

    if generated_str != committed {
        // Write the generated version for diffing
        let tmp = crate_dir.join("include/swisseph.h.generated");
        std::fs::write(&tmp, &generated_str).ok();
        panic!(
            "include/swisseph.h is stale. Regenerate with:\n  \
             cbindgen --config cbindgen.toml --crate swisseph-ffi --output include/swisseph.h\n\
             Diff written to {:?}",
            tmp
        );
    }
}
