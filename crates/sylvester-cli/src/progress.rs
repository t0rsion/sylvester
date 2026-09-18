//! Optional phase messages for long-running CLI commands.

use clap::Args;

/// Controls phase messages on standard error.
#[derive(Args, Clone, Copy, Debug, Default)]
pub(crate) struct ProgressArgs {
    /// Write phase messages to standard error.
    #[arg(long, conflicts_with = "no_progress")]
    pub progress: bool,
    /// Suppress phase messages.
    #[arg(long, conflicts_with = "progress")]
    pub no_progress: bool,
}

impl ProgressArgs {
    /// Build a reporter from the command-line switches.
    pub(crate) fn reporter(self) -> Reporter {
        Reporter {
            enabled: self.progress && !self.no_progress,
        }
    }
}

/// Emits named phases without estimating work that has no known bound.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Reporter {
    enabled: bool,
}

impl Reporter {
    /// Report a phase when progress was requested.
    pub(crate) fn phase(self, name: &str) {
        if self.enabled {
            eprintln!("phase: {name}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::ProgressArgs;

    #[test]
    fn progress_is_opt_in() {
        assert!(!ProgressArgs::default().reporter().enabled);
        assert!(
            ProgressArgs {
                progress: true,
                no_progress: false,
            }
            .reporter()
            .enabled
        );
    }
}
