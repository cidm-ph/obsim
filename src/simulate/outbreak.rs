use crate::case::History;
use crate::genome::Genome;
use crate::{Count, Time};
use std::io;

/// A simulated outbreak containing a number of cases.
#[derive(Debug)]
pub struct Outbreak<G> {
    pub(super) source: Vec<Option<Count>>,
    pub(super) history: Vec<History>,
    pub(super) genome: Vec<G>,
}

impl<G: Genome> Outbreak<G> {
    /// Get the source (infector) of all cases.
    ///
    /// The values are `None` when a case has no infector. This occurs for the index case in each
    /// simulation, but there can be multiple cases with no infector when outbreaks are combined
    /// with [`extend_with`](Outbreak::extend_with).
    #[inline]
    pub fn sources(&self) -> &[Option<Count>] {
        &self.source
    }

    /// Get the outbreak number of all cases.
    ///
    /// For a single simulation, this will be zero for all cases. When multiple outbreaks are
    /// combined with `extend_with` this will return distinct values for cases that originated in
    /// different simulations.
    #[inline]
    pub fn outbreaks(&self) -> Vec<Count> {
        get_cluster_ids(&self.source)
    }

    /// Get the disease history times of all cases.
    ///
    #[inline]
    pub fn history(&self) -> &[History] {
        &self.history
    }

    /// Get the genomes of all cases.
    #[inline]
    pub fn genomes(&self) -> &[G] {
        &self.genome
    }

    /// Print a FASTA file representing the simulated genomes.
    pub fn write_fasta<W: io::Write>(&self, mut writer: W) -> io::Result<()> {
        let sources = self.outbreaks();

        for (i, genome) in self.genome.iter().enumerate() {
            writeln!(
                writer,
                ">case{:06} day_infected={} day_reported={} outbreak={} parent={}",
                i,
                self.history[i].infected,
                self.history[i].infected,
                sources[i],
                self.source[i]
                    .map(|x| format!("case{:06}", x))
                    .unwrap_or_else(String::new)
            )?;
            genome.write_nucleotides(&mut writer)?;
            writeln!(writer)?;
        }
        Ok(())
    }

    /// Modify all times such that the earliest infection occurs at time zero.
    pub fn rezero_time(&mut self) {
        let start_time = self.history.iter().flat_map(|x| x.iter().min()).min();
        let by_amount = start_time.unwrap_or(0);
        for history in &mut self.history {
            history.time_shift_back(by_amount);
        }
    }

    /// The latest time stored in the outbreak metadata.
    #[inline]
    pub fn end_time(&self) -> Option<Time> {
        self.history.iter().flat_map(|x| x.iter().max()).max()
    }

    /// The number of cases simulated.
    #[inline]
    pub fn n_cases(&self) -> usize {
        self.source.len()
    }

    /// Increase all times in this outbreak by a fixed amount.
    pub fn time_shift(&mut self, by_amount: Time) {
        for history in &mut self.history {
            history.time_shift_forward(by_amount);
        }
    }

    fn id_shift(&mut self, by_amount: Count) {
        for source in self.source.iter_mut().flatten() {
            *source += by_amount;
        }
    }

    /// Append all cases from `other`, shifting their IDs to avoid collisions.
    pub fn extend_with(&mut self, mut other: Outbreak<G>) {
        other.id_shift(self.source.len() as Count);
        self.source.extend(other.source);
        self.history.extend(other.history);
        self.genome.extend(other.genome);
    }
}

/// Convert a vector of sources into a vector of cluster IDs.
pub fn get_cluster_ids(sources: &[Option<Count>]) -> Vec<Count> {
    // NOTE this assumes that all Some(_) entries belong to the same cluster as
    // the previous None. This is true the way the simulation currently works but
    // is not true for arbitrary well-formed source vectors.
    sources
        .iter()
        .scan(0, |state, source| {
            if source.is_none() {
                *state += 1;
            }
            Some(*state - 1)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_cluster_ids() {
        let sources = &[None, Some(0), Some(1), None, Some(3), None];
        assert_eq!(get_cluster_ids(sources), vec![0, 0, 0, 1, 1, 2]);
    }
}
