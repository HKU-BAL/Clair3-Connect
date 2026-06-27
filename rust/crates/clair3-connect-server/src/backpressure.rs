use std::sync::Arc;
use std::time::Duration;
use tokio::sync::{OwnedSemaphorePermit, Semaphore};
use tonic::Status;

#[derive(Debug, Clone)]
pub struct BackpressureConfig {

    pub max_inflight: usize,

    pub retry_after: Duration,
}

impl Default for BackpressureConfig {
    fn default() -> Self {
        Self {
            max_inflight: 8,
            retry_after: Duration::from_millis(100),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Gate {
    permits: Arc<Semaphore>,
    retry_after_ms: u64,

    max_inflight: usize,
}

impl Gate {

    pub fn new(cfg: BackpressureConfig) -> Self {

        let max = cfg.max_inflight.min(Semaphore::MAX_PERMITS);
        Self {
            permits: Arc::new(Semaphore::new(max)),
            retry_after_ms: cfg.retry_after.as_millis().min(u64::MAX as u128) as u64,
            max_inflight: max,
        }
    }

    #[allow(clippy::result_large_err)]
    pub fn try_acquire(&self) -> Result<OwnedSemaphorePermit, Status> {
        match self.permits.clone().try_acquire_owned() {
            Ok(p) => Ok(p),
            Err(_) => Err(self.exhausted_status()),
        }
    }

    pub fn max_inflight(&self) -> usize {
        self.max_inflight
    }

    pub fn available_permits(&self) -> usize {
        self.permits.available_permits()
    }

    fn exhausted_status(&self) -> Status {
        let mut s = Status::resource_exhausted(format!(
            "server inflight cap {} reached; retry in {} ms",
            self.max_inflight, self.retry_after_ms
        ));

        let v = self
            .retry_after_ms
            .to_string()
            .parse()
            .expect("u64 decimal digits are valid ASCII metadata");
        s.metadata_mut().insert("retry-after-ms", v);
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gate_with_zero_capacity_rejects_everything() {
        let g = Gate::new(BackpressureConfig {
            max_inflight: 0,
            retry_after: Duration::from_millis(50),
        });
        let err = g.try_acquire().unwrap_err();
        assert_eq!(err.code(), tonic::Code::ResourceExhausted);
        let hint = err
            .metadata()
            .get("retry-after-ms")
            .expect("retry-after-ms set")
            .to_str()
            .unwrap()
            .to_string();
        assert_eq!(hint, "50");
    }

    #[test]
    fn gate_acquire_release_acquire() {
        let g = Gate::new(BackpressureConfig {
            max_inflight: 1,
            retry_after: Duration::from_millis(10),
        });
        let p1 = g.try_acquire().expect("first permit");

        let err = g.try_acquire().unwrap_err();
        assert_eq!(err.code(), tonic::Code::ResourceExhausted);
        assert_eq!(g.available_permits(), 0);
        drop(p1);

        assert_eq!(g.available_permits(), 1);
        let _p2 = g.try_acquire().expect("second permit after drop");
    }

    #[test]
    fn gate_concurrent_acquires_up_to_cap() {
        let g = Gate::new(BackpressureConfig {
            max_inflight: 3,
            retry_after: Duration::from_millis(1),
        });
        let p1 = g.try_acquire().unwrap();
        let p2 = g.try_acquire().unwrap();
        let p3 = g.try_acquire().unwrap();

        assert!(g.try_acquire().is_err());

        drop(p2);
        let _p4 = g.try_acquire().expect("permit available after drop");

        assert!(g.try_acquire().is_err());
        drop(p1);
        drop(p3);
    }

    #[test]
    fn gate_default_is_8() {
        let g = Gate::new(BackpressureConfig::default());
        assert_eq!(g.max_inflight(), 8);
        assert_eq!(g.available_permits(), 8);
    }
}
