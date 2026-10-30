//! Optical duplicate tests (`-d`).

use super::{
    Check,
    names::{CoordinateParser, Coordinates},
};

/// Whether two coordinates are within `distance` of each other on both axes, given matching
/// prefixes.
pub(super) fn is_close(original: &Coordinates, x: i64, y: i64, distance: i64) -> bool {
    (original.x - x).abs() <= distance && (original.y - y).abs() <= distance
}

/// `is_optical_duplicate`
pub(super) fn is_optical_duplicate(
    parser: &CoordinateParser,
    original: &[u8],
    duplicate: &[u8],
    distance: i64,
    warnings: &mut u64,
) -> bool {
    let Ok(o) = parser.parse(original, warnings) else {
        return false;
    };

    let Ok(d) = parser.parse(duplicate, warnings) else {
        return false;
    };

    original[o.prefix.clone()] == duplicate[d.prefix] && is_close(&o, d.x, d.y, distance)
}

/// `optical_duplicate_partial`: like [`is_optical_duplicate`] with the original's coordinates
/// already parsed, also recording the duplicate's coordinates in `check`.
pub(super) fn optical_duplicate_partial(
    parser: &CoordinateParser,
    original_name: &[u8],
    original: &Coordinates,
    duplicate: &[u8],
    check: &mut Check,
    distance: i64,
    warnings: &mut u64,
) -> bool {
    let Ok(d) = parser.parse(duplicate, warnings) else {
        return false;
    };

    let is_optical = original_name[original.prefix.clone()] == duplicate[d.prefix.clone()]
        && is_close(original, d.x, d.y, distance);

    check.x = d.x;
    check.y = d.y;
    check.has_coordinates = true;
    check.prefix = duplicate[d.prefix].to_vec();

    is_optical
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_optical_duplicate() {
        let mut warnings = 0;
        let original = b"M:7:FC:1:1101:1000:2000";

        assert!(is_optical_duplicate(
            &CoordinateParser::Colons,
            original,
            b"M:7:FC:1:1101:1100:1900",
            100,
            &mut warnings
        ));
        assert!(!is_optical_duplicate(
            &CoordinateParser::Colons,
            original,
            b"M:7:FC:1:1101:1101:2000",
            100,
            &mut warnings
        ));
        assert!(!is_optical_duplicate(
            &CoordinateParser::Colons,
            original,
            b"M:7:FC:1:1102:1000:2000",
            100,
            &mut warnings
        ));
        assert_eq!(warnings, 0);
    }
}
