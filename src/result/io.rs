use std::io::{self, Read, Write};

use crate::compute::ComputeLimits;

use super::{EnvelopeError, check};

pub(super) const CHUNK: usize = 8 * 1024;

fn stopped_io(error: &EnvelopeError) -> io::Error {
    io::Error::other(error.to_string())
}

pub(super) struct BudgetReader<'a> {
    bytes: &'a [u8],
    position: usize,
    limits: &'a ComputeLimits,
    stopped: Option<EnvelopeError>,
}

impl<'a> BudgetReader<'a> {
    pub(super) fn new(bytes: &'a [u8], limits: &'a ComputeLimits) -> Self {
        Self {
            bytes,
            position: 0,
            limits,
            stopped: None,
        }
    }

    pub(super) fn take_stop(&mut self) -> Option<EnvelopeError> {
        self.stopped.take()
    }

    fn poll(&mut self) -> io::Result<()> {
        if let Some(error) = &self.stopped {
            return Err(stopped_io(error));
        }
        match check(self.limits, 0) {
            Ok(()) => Ok(()),
            Err(error) => {
                self.stopped = Some(error.clone());
                Err(stopped_io(&error))
            }
        }
    }
}

impl Read for BudgetReader<'_> {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        self.poll()?;
        if output.is_empty() || self.position == self.bytes.len() {
            return Ok(0);
        }

        let count = output
            .len()
            .min(CHUNK)
            .min(self.bytes.len().saturating_sub(self.position));
        output[..count].copy_from_slice(&self.bytes[self.position..self.position + count]);
        self.position += count;
        self.poll()?;
        Ok(count)
    }
}

pub(super) struct BudgetWriter<'a> {
    bytes: Vec<u8>,
    limits: &'a ComputeLimits,
    base_bytes: usize,
    stopped: Option<EnvelopeError>,
}

impl<'a> BudgetWriter<'a> {
    pub(super) fn new(limits: &'a ComputeLimits, base_bytes: usize) -> Self {
        Self {
            bytes: Vec::new(),
            limits,
            base_bytes,
            stopped: None,
        }
    }

    pub(super) fn capacity(&self) -> usize {
        self.bytes.capacity()
    }

    pub(super) fn take_stop(&mut self) -> Option<EnvelopeError> {
        self.stopped.take()
    }

    pub(super) fn into_inner(self) -> Vec<u8> {
        self.bytes
    }

    fn poll(&mut self) -> io::Result<()> {
        if let Some(error) = &self.stopped {
            return Err(stopped_io(error));
        }
        match check(
            self.limits,
            self.base_bytes.saturating_add(self.bytes.capacity()),
        ) {
            Ok(()) => Ok(()),
            Err(error) => {
                self.stopped = Some(error.clone());
                Err(stopped_io(&error))
            }
        }
    }
}

impl Write for BudgetWriter<'_> {
    fn write(&mut self, input: &[u8]) -> io::Result<usize> {
        self.poll()?;
        if input.is_empty() {
            return Ok(0);
        }

        let count = input.len().min(CHUNK);
        let spare = self.bytes.capacity().saturating_sub(self.bytes.len());
        if count > spare {
            let needed = self.bytes.len().saturating_add(count);
            match check(self.limits, self.base_bytes.saturating_add(needed)) {
                Ok(()) => {}
                Err(error) => {
                    self.stopped = Some(error.clone());
                    return Err(stopped_io(&error));
                }
            }
            if self.bytes.try_reserve_exact(count).is_err() {
                let error = EnvelopeError::MemoryLimitExceeded;
                self.stopped = Some(error.clone());
                return Err(stopped_io(&error));
            }
        }
        self.bytes.extend_from_slice(&input[..count]);
        self.poll()?;
        Ok(count)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.poll()
    }
}

#[cfg(test)]
mod tests {
    use std::io::{Read, Write};

    use super::{BudgetReader, BudgetWriter};
    use crate::{Budget, CancellationToken};

    #[test]
    fn reader_reports_cancellation_between_chunks() {
        let cancellation = CancellationToken::new();
        let budget = Budget::new().cancellation(cancellation.clone());
        let limits = crate::compute::ComputeLimits::of_budget(&budget);
        let mut reader = BudgetReader::new(b"abcdef", &limits);
        let mut output = [0_u8; 1];

        assert_eq!(reader.read(&mut output).expect("first chunk"), 1);
        cancellation.cancel();
        assert!(reader.read(&mut output).is_err());
        assert!(
            reader
                .take_stop()
                .is_some_and(|error| matches!(error, crate::result::EnvelopeError::Timeout))
        );
    }

    #[test]
    fn writer_reports_cancellation_between_chunks() {
        let cancellation = CancellationToken::new();
        let budget = Budget::new().cancellation(cancellation.clone());
        let limits = crate::compute::ComputeLimits::of_budget(&budget);
        let mut writer = BudgetWriter::new(&limits, 0);

        writer.write_all(b"abcdef").expect("first chunk");
        cancellation.cancel();
        assert!(writer.write_all(b"x").is_err());
        assert!(
            writer
                .take_stop()
                .is_some_and(|error| matches!(error, crate::result::EnvelopeError::Timeout))
        );
    }

    #[test]
    fn writer_reserves_a_growing_chunk_before_extending() {
        let budget = Budget::new().memory_limit(25);
        let limits = crate::compute::ComputeLimits::of_budget(&budget);
        let mut writer = BudgetWriter::new(&limits, 0);
        writer.bytes.reserve_exact(16);

        writer.write_all(b"0123456789").expect("first chunk");
        writer.write_all(b"abcdefghij").expect("second chunk");
        assert!(writer.capacity() <= 25);
    }

    #[test]
    fn writer_rejects_a_chunk_before_allocation_at_a_tight_cap() {
        let budget = Budget::new().memory_limit(19);
        let limits = crate::compute::ComputeLimits::of_budget(&budget);
        let mut writer = BudgetWriter::new(&limits, 0);
        writer.bytes.reserve_exact(16);

        writer.write_all(b"0123456789").expect("first chunk");
        let capacity = writer.capacity();
        assert!(writer.write_all(b"abcdefghij").is_err());
        assert_eq!(writer.capacity(), capacity);
        assert!(writer.take_stop().is_some_and(|error| matches!(
            error,
            crate::result::EnvelopeError::MemoryLimitExceeded
        )));
    }
}
