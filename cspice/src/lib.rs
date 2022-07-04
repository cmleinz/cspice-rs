pub mod cell;
pub mod common;
pub mod coordinates;
pub mod data;
pub mod error;
pub mod spk;
pub mod string;
pub mod time;
pub mod vector;

use crate::error::Error;
use crate::string::SpiceString;
use once_cell::sync::OnceCell;
use std::cell::Cell;
use std::fmt::Debug;
use std::thread;
use std::thread::Thread;
use thiserror::Error;

/// Wraps an unsafe SPICE function call.
///
/// First checks that it is safe for the current thread to access SPICE, otherwise panics.
macro_rules! spice_unsafe {
    ($l:block) => {{
        crate::Spice::try_acquire_thread().unwrap();
        unsafe { $l }
    }};
}
pub(crate) use spice_unsafe;

/// Namespace that groups functions that read and write to the SPICE runtime.
pub struct Spice;

impl Spice {
    /// The SPICE library [is not thread safe](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/problems.html#Problem:%20SPICE%20code%20is%20not%20thread%20safe).
    ///
    /// This function checks if it is safe for the current thread to call SPICE functions. SPICE
    /// will be locked to the first thread that calls this function.
    pub fn try_acquire_thread() -> Result<(), SpiceThreadError> {
        static SPICE_THREAD_ID: OnceCell<Thread> = OnceCell::new();
        thread_local! {
            static CACHE: Cell<bool> = Cell::new(false);
        }
        CACHE.with(|cache| {
            if cache.get() {
                return Ok(());
            }
            match SPICE_THREAD_ID.set(thread::current()) {
                Ok(_) => {
                    cache.set(true);
                    Spice::set_error_defaults();
                    Ok(())
                }
                Err(e) => Err(SpiceThreadError(e)),
            }
        })
    }
}

#[derive(Debug, Clone, Error)]
#[error("SPICE is already in use by another thread")]
pub struct SpiceThreadError(pub Thread);

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::sync::Once;

    /// Load test data (once)
    pub fn load_test_data() {
        static SPICE_INIT: Once = Once::new();
        SPICE_INIT.call_once(|| {
            let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_data");
            Spice::furnish(data_dir.join("naif0012.tls").to_string_lossy()).unwrap();
        });
    }

    #[test]
    fn test_acquire_thread() {
        Spice::try_acquire_thread().unwrap();
        Spice::try_acquire_thread().unwrap();
    }

    #[test]
    fn test_acquire_thread_different_thread() {
        Spice::try_acquire_thread().unwrap();
        std::thread::spawn(|| {
            Spice::try_acquire_thread().expect_err("Should be unable to use on another thread")
        })
        .join()
        .unwrap();
    }
}
