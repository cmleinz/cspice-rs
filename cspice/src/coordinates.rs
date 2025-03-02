//! Functions for converting between different types of coordinates.
use crate::spice_unsafe;
use cspice_sys::{azlrec_c, recazl_c, reclat_c, recrad_c, SpiceBoolean, SpiceDouble};
use derive_more::Into;

/// Rectangular coordinates
#[repr(C)]
#[derive(Copy, Clone, Debug, Default, PartialEq, Into)]
pub struct Rectangular {
    pub x: SpiceDouble,
    pub y: SpiceDouble,
    pub z: SpiceDouble,
}

impl From<[SpiceDouble; 3]> for Rectangular {
    fn from(rect: [SpiceDouble; 3]) -> Self {
        Rectangular {
            x: rect[0],
            y: rect[1],
            z: rect[2],
        }
    }
}

impl From<Rectangular> for [SpiceDouble; 3] {
    fn from(rect: Rectangular) -> Self {
        [rect.x, rect.y, rect.z]
    }
}

/// Range, azimuth, and elevation
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct AzEl {
    pub range: SpiceDouble,
    pub az: SpiceDouble,
    pub el: SpiceDouble,
}

impl AzEl {
    /// See [recazl_c](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/recazl_c.html)
    pub fn from_rect(mut rect: Rectangular, azccw: bool, elplsz: bool) -> Self {
        let mut az_el = AzEl::default();
        spice_unsafe!({
            recazl_c(
                &mut rect.x as *mut SpiceDouble,
                azccw as SpiceBoolean,
                elplsz as SpiceBoolean,
                &mut az_el.range,
                &mut az_el.az,
                &mut az_el.el,
            )
        });
        az_el
    }
}

impl Rectangular {
    /// See [azlrec_c](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/azlrec_c.html)
    pub fn from_azel(azel: AzEl, azccw: bool, elplsz: bool) -> Self {
        let mut rect = [0.0f64; 3];
        spice_unsafe!({
            azlrec_c(
                azel.range,
                azel.az,
                azel.el,
                azccw as SpiceBoolean,
                elplsz as SpiceBoolean,
                rect.as_mut_ptr(),
            )
        });
        rect.into()
    }
}

/// Range, right ascension, and declination.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct RaDec {
    pub range: SpiceDouble,
    pub ra: SpiceDouble,
    pub dec: SpiceDouble,
}

impl From<Rectangular> for RaDec {
    /// See [recrad_c](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/recrad_c.html).
    fn from(mut rect: Rectangular) -> Self {
        let mut ra_dec = RaDec::default();
        spice_unsafe!({
            recrad_c(
                &mut rect.x as *mut SpiceDouble,
                &mut ra_dec.range,
                &mut ra_dec.ra,
                &mut ra_dec.dec,
            )
        });
        ra_dec
    }
}

/// Latitudinal coordinates.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Latitudinal {
    pub radius: SpiceDouble,
    pub longitude: SpiceDouble,
    pub latitude: SpiceDouble,
}

impl From<Rectangular> for Latitudinal {
    /// See [reclat_c](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/reclat_c.html).
    fn from(mut rect: Rectangular) -> Self {
        let mut lat = Latitudinal::default();
        spice_unsafe!({
            reclat_c(
                &mut rect.x as *mut SpiceDouble,
                &mut lat.radius,
                &mut lat.longitude,
                &mut lat.latitude,
            )
        });
        lat
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-3;

    // Test data comes from NAIF website https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/recazl_c.html
    const TEST_DATA_F_F: [[SpiceDouble; 6]; 11] = [
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [1.000, 0.000, 0.000, 1.000, 0.000, 0.000],
        [0.000, 1.000, 0.000, 1.000, 270.000, 0.000],
        [0.000, 0.000, 1.000, 1.000, 0.000, -90.000],
        [-1.000, 0.000, 0.000, 1.000, 180.000, 0.000],
        [0.000, -1.000, 0.000, 1.000, 90.000, 0.000],
        [0.000, 0.000, -1.000, 1.000, 0.000, 90.000],
        [1.000, 1.000, 0.000, 1.414, 315.000, 0.000],
        [1.000, 0.000, 1.000, 1.414, 0.000, -45.000],
        [0.000, 1.000, 1.000, 1.414, 270.000, -45.000],
        [1.000, 1.000, 1.000, 1.732, 315.000, -35.264],
    ];

    const TEST_DATA_F_T: [[SpiceDouble; 6]; 11] = [
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [1.000, 0.000, 0.000, 1.000, 0.000, 0.000],
        [0.000, 1.000, 0.000, 1.000, 270.000, 0.000],
        [0.000, 0.000, 1.000, 1.000, 0.000, 90.000],
        [-1.000, 0.000, 0.000, 1.000, 180.000, 0.000],
        [0.000, -1.000, 0.000, 1.000, 90.000, 0.000],
        [0.000, 0.000, -1.000, 1.000, 0.000, -90.000],
        [1.000, 1.000, 0.000, 1.414, 315.000, 0.000],
        [1.000, 0.000, 1.000, 1.414, 0.000, 45.000],
        [0.000, 1.000, 1.000, 1.414, 270.000, 45.000],
        [1.000, 1.000, 1.000, 1.732, 315.000, 35.264],
    ];

    const TEST_DATA_T_F: [[SpiceDouble; 6]; 11] = [
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [1.000, 0.000, 0.000, 1.000, 0.000, 0.000],
        [0.000, 1.000, 0.000, 1.000, 90.000, 0.000],
        [0.000, 0.000, 1.000, 1.000, 0.000, -90.000],
        [-1.000, 0.000, 0.000, 1.000, 180.000, 0.000],
        [0.000, -1.000, 0.000, 1.000, 270.000, 0.000],
        [0.000, 0.000, -1.000, 1.000, 0.000, 90.000],
        [1.000, 1.000, 0.000, 1.414, 45.000, 0.000],
        [1.000, 0.000, 1.000, 1.414, 0.000, -45.000],
        [0.000, 1.000, 1.000, 1.414, 90.000, -45.000],
        [1.000, 1.000, 1.000, 1.732, 45.000, -35.264],
    ];

    const TEST_DATA_T_T: [[SpiceDouble; 6]; 11] = [
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [1.000, 0.000, 0.000, 1.000, 0.000, 0.000],
        [0.000, 1.000, 0.000, 1.000, 90.000, 0.000],
        [0.000, 0.000, 1.000, 1.000, 0.000, 90.000],
        [-1.000, 0.000, 0.000, 1.000, 180.000, 0.000],
        [0.000, -1.000, 0.000, 1.000, 270.000, 0.000],
        [0.000, 0.000, -1.000, 1.000, 0.000, -90.000],
        [1.000, 1.000, 0.000, 1.414, 45.000, 0.000],
        [1.000, 0.000, 1.000, 1.414, 0.000, 45.000],
        [0.000, 1.000, 1.000, 1.414, 90.000, 45.000],
        [1.000, 1.000, 1.000, 1.732, 45.000, 35.264],
    ];

    #[test]
    fn test_azel_rect_conversion() {
        azel_rect_conversion(&TEST_DATA_F_F, false, false);
        azel_rect_conversion(&TEST_DATA_F_T, false, true);
        azel_rect_conversion(&TEST_DATA_T_F, true, false);
        azel_rect_conversion(&TEST_DATA_T_T, true, true);
    }

    fn azel_rect_conversion(test_data: &[[f64; 6]; 11], azccw: bool, elplsz: bool) {
        for test in test_data.iter() {
            let azel = AzEl {
                range: test[3],
                az: test[4].to_radians(),
                el: test[5].to_radians(),
            };
            let rect = Rectangular::from_azel(azel, azccw, elplsz);
            assert!((rect.x - test[0]).abs() < EPSILON);
            assert!((rect.y - test[1]).abs() < EPSILON);
            assert!((rect.z - test[2]).abs() < EPSILON);
            let azel_ = AzEl::from_rect(rect, azccw, elplsz);
            assert!((azel_.range - test[3]).abs() < EPSILON);
            assert!((azel_.az - test[4].to_radians()).abs() < EPSILON);
            assert!((azel_.el - test[5].to_radians()).abs() < EPSILON);
        }
    }
}
