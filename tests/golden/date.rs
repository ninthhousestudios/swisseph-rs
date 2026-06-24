use serde::Deserialize;
use swisseph::CalendarType;
use swisseph::date::{date_conversion, day_of_week, julday, revjul};

#[derive(Deserialize)]
struct JuldayCase {
    y: i32,
    m: i32,
    d: i32,
    h: f64,
    g: i32,
    jd: f64,
}

#[derive(Deserialize)]
struct RevjulCase {
    jd: f64,
    g: i32,
    y: i32,
    m: i32,
    d: i32,
    h: f64,
}

#[derive(Deserialize)]
struct DateConvCase {
    y: i32,
    m: i32,
    d: i32,
    h: f64,
    g: i32,
    valid: bool,
    jd: Option<f64>,
}

#[derive(Deserialize)]
struct DowCase {
    jd: f64,
    dow: u8,
}

#[derive(Deserialize)]
struct DateGolden {
    julday: Vec<JuldayCase>,
    revjul: Vec<RevjulCase>,
    date_conversion: Vec<DateConvCase>,
    day_of_week: Vec<DowCase>,
}

fn load() -> DateGolden {
    let path = super::golden_data_path("date.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

#[test]
fn golden_julday() {
    let data = load();
    for (i, c) in data.julday.iter().enumerate() {
        let cal = CalendarType::try_from(c.g).unwrap();
        let actual = julday(c.y, c.m, c.d, c.h, cal);
        super::assert_f64_exact(
            &format!("julday[{i}]({}, {}, {}, {}, {cal:?})", c.y, c.m, c.d, c.h),
            c.jd,
            actual,
        );
    }
}

#[test]
fn golden_revjul() {
    let data = load();
    for (i, c) in data.revjul.iter().enumerate() {
        let cal = CalendarType::try_from(c.g).unwrap();
        let (y, m, d, h) = revjul(c.jd, cal);
        let label = format!("revjul[{i}](jd={}, {cal:?})", c.jd);
        assert_eq!(y, c.y, "{label}: year");
        assert_eq!(m, c.m, "{label}: month");
        assert_eq!(d, c.d, "{label}: day");
        super::assert_f64_exact(&format!("{label}: hour"), c.h, h);
    }
}

#[test]
fn golden_date_conversion() {
    let data = load();
    for (i, c) in data.date_conversion.iter().enumerate() {
        let cal = CalendarType::try_from(c.g).unwrap();
        let label = format!(
            "date_conversion[{i}]({}, {}, {}, {}, {cal:?})",
            c.y, c.m, c.d, c.h
        );
        let result = date_conversion(c.y, c.m, c.d, c.h, cal);
        if c.valid {
            let jd = result.unwrap_or_else(|e| panic!("{label}: expected Ok, got Err({e:?})"));
            super::assert_f64_exact(&label, c.jd.unwrap(), jd);
        } else {
            assert!(
                result.is_err(),
                "{label}: expected Err, got Ok({:?})",
                result.unwrap()
            );
        }
    }
}

#[test]
fn golden_day_of_week() {
    let data = load();
    for (i, c) in data.day_of_week.iter().enumerate() {
        let actual = day_of_week(c.jd);
        assert_eq!(
            actual, c.dow,
            "day_of_week[{i}](jd={}): expected {}, got {actual}",
            c.jd, c.dow
        );
    }
}
