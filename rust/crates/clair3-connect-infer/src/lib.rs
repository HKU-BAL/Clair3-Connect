#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

pub mod model_registry;

pub mod torchscript;

pub mod worker;

pub mod batcher;

pub mod device;

#[doc(hidden)]
pub fn _ensure_cuda_linked() -> i64 {
    #[cfg(target_os = "linux")]
    unsafe {
        use std::sync::Once;
        static ONCE: Once = Once::new();
        ONCE.call_once(|| {

            let libname = std::ffi::CString::new("libtorch_cuda.so").unwrap();
            let handle = libc::dlopen(
                libname.as_ptr(),
                libc::RTLD_LAZY | libc::RTLD_GLOBAL,
            );
            if handle.is_null() {
                let err = libc::dlerror();
                if !err.is_null() {
                    let msg = std::ffi::CStr::from_ptr(err);
                    tracing::warn!(
                        "dlopen libtorch_cuda.so failed: {}",
                        msg.to_string_lossy()
                    );
                }
            }
        });
    }
    tch::Cuda::device_count()
}
