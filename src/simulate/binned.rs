use rand::Rng;
use rand_distr::{Distribution, Uniform};
use std::collections::{BTreeMap, VecDeque};
use std::fmt;
use thiserror::Error;

use super::{simulate_outbreak, Outbreak};
use crate::disease::DiseaseModel;
use crate::genome::Genome;
use crate::{Count, Time};

/// Configuration for [`binned_outbreaks()`].
#[derive(Debug, Clone)]
pub struct BinnedOutbreakConfig {
    /// Edges of the size bins: provide one more edge than the number of bins.
    pub size_bin_edges: Vec<Count>,

    /// Capacity of each size bin: provide one per bin.
    pub size_counts: Vec<Count>,

    /// Length of the index case window.
    pub latest_importation: Time,

    /// Length of the ancestral divergence time.
    pub time_to_mrca: Time,

    /// Length of the ancestral divergence time for singletons.
    pub time_to_background_mrca: Time,

    /// Number of singletons.
    pub n_background: Count,

    /// Maximum number of simulations to reject.
    pub bad_simulation_cap: usize,
}

#[derive(Error, Debug)]
pub struct BinError {
    discarded: Vec<usize>,
    remaining: Vec<u32>,
    config: BinnedOutbreakConfig,
}

/// Generate and merge many outbreaks.
///
/// Specify a number of size bins. Outbreaks will be generated and if the number of cases fits into
/// one of the bins, that simulation will be retained. New outbreaks are generated until all of the
/// size bins are filled.
///
/// The outbreaks are shifted so that their index cases fall uniformly at random on the index case
/// window. Their genomes are mutated according to the time since an ancestral genome. These
/// parameters are shown in the schematic below.
///
#[doc=include_str!("../../sim_params.svg")]
///
/// See `cargo run --example binned` for an example.
///
/// Panics if `bad_simulation_cap` outbreaks are rejected, since this probably means that the
/// simulation parameters are producing outbreaks that are too big or too small for the size bin
/// configuration.
pub fn binned_outbreaks<D, G, R>(
    ancestral_genome: G,
    disease_model: &D,
    mutation_rate: f64,
    sim_config: &BinnedOutbreakConfig,
    mut rng: R,
) -> Result<Outbreak<G>, BinError>
where
    D: DiseaseModel,
    G: Genome,
    R: Rng,
{
    sim_config.validate();
    let mut size_counts = sim_config.size_counts.clone();

    let total_signal: u32 = size_counts.iter().sum();

    let mut failed = Vec::new();
    let mut outbreak = Outbreak {
        source: Vec::new(),
        history: Vec::new(),
        genome: Vec::new(),
    };

    let importation_dist = Uniform::from(Time::default()..=sim_config.latest_importation);
    let mut importation_times: VecDeque<Time> = importation_dist
        .sample_iter(&mut rng)
        .take(total_signal as usize)
        .collect();

    while size_counts.iter().sum::<u32>() > 0 {
        let generation_time = importation_times[0] + sim_config.time_to_mrca;
        let genome = ancestral_genome.mutate_time(generation_time, mutation_rate, &mut rng);

        match simulate_outbreak(
            genome,
            disease_model,
            mutation_rate,
            sim_config.max_size(),
            &mut rng,
        ) {
            Ok(mut new_ob) => {
                let size_bin = sim_config.size_bin(new_ob.n_cases() as Count);
                if accept(&mut size_counts, size_bin) {
                    new_ob.time_shift(importation_times.pop_front().unwrap());
                    outbreak.extend_with(new_ob);
                } else {
                    failed.push(new_ob.n_cases());
                }
            }
            Err(new_ob) => {
                failed.push(new_ob.outbreak.n_cases());
            }
        }

        if failed.len() >= sim_config.bad_simulation_cap {
            return Err(BinError {
                discarded: failed,
                remaining: size_counts,
                config: sim_config.clone(),
            });
        }
    }

    let last_case = outbreak.end_time().unwrap_or_default();
    let importation_dist = Uniform::from(Time::default()..=last_case);
    for _ in 0..sim_config.n_background {
        let imported_at = importation_dist.sample(&mut rng);
        let generation_time = imported_at + sim_config.time_to_background_mrca;
        let genome = ancestral_genome.mutate_time(generation_time, mutation_rate, &mut rng);

        let (_, history) = disease_model
            .generate_singleton(&mut rng)
            .into_case_history();
        outbreak.source.push(None);
        outbreak.history.push(history);
        outbreak.genome.push(genome);
    }

    outbreak.rezero_time();

    Ok(outbreak)
}

fn accept(size_counts: &mut [Count], size_bin: Option<usize>) -> bool {
    if let Some(bin) = size_bin {
        if size_counts[bin] > 0 {
            size_counts[bin] -= 1;
            return true;
        }
    };
    false
}

impl BinnedOutbreakConfig {
    fn validate(&self) {
        let n = self.size_bin_edges.len();
        assert_eq!(self.size_counts.len(), n - 1);
        for i in 0..(n - 1) {
            assert!(self.size_bin_edges[i + 1] > self.size_bin_edges[i]);
        }
    }

    fn max_size(&self) -> u32 {
        *self.size_bin_edges.last().unwrap()
    }

    fn size_bin(&self, n_cases: Count) -> Option<usize> {
        self.size_bin_edges
            .windows(2)
            .position(|x| x[0] < n_cases && x[1] >= n_cases)
    }

    fn size_bin_labels(&self) -> Vec<String> {
        self.size_bin_edges
            .windows(2)
            .map(|x| format!("({},{}]", x[0], x[1]))
            .collect()
    }

    fn size_bin_with_margins(&self, n_cases: Count) -> usize {
        if self.size_bin_edges.is_empty() || n_cases <= self.size_bin_edges[0] {
            return 0;
        }
        if n_cases > *self.size_bin_edges.last().unwrap() {
            return self.size_bin_edges.len();
        }

        self.size_bin_edges
            .windows(2)
            .position(|x| x[0] < n_cases && x[1] >= n_cases)
            .unwrap()
            + 1
    }
}

impl fmt::Display for BinError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "simulation did not fill all outbreak size bins after {} attempts",
            self.config.bad_simulation_cap
        )?;

        let mut labels = self.config.size_bin_labels();
        labels.insert(0, "underflow".to_owned());
        labels.push("overflow".to_owned());

        let failed_bins = self
            .discarded
            .iter()
            .map(|n| self.config.size_bin_with_margins(*n as u32))
            .fold(BTreeMap::new(), |mut acc, val| {
                acc.entry(val).and_modify(|e| *e += 1).or_insert(1);
                acc
            });
        let failed_bins: Vec<_> = (0..labels.len())
            .map(|x| *failed_bins.get(&x).unwrap_or(&0))
            .collect();

        writeln!(
            f,
            "{:>15}   {:>6}   {:>6}   {:>6}",
            "BIN", "CFG", "TODO", "EXTRA"
        )?;
        for i in 0..labels.len() {
            let (cfg, todo) = if i == 0 || i == labels.len() - 1 {
                ("".to_owned(), "".to_owned())
            } else {
                (
                    self.config.size_counts[i - 1].to_string(),
                    self.remaining[i - 1].to_string(),
                )
            };
            writeln!(
                f,
                "{:>15}   {:>6}   {:>6}   {:>6}",
                labels[i], cfg, todo, failed_bins[i]
            )?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::disease::simple::SimpleDisease;
    use crate::genome::simple::SimpleGenome;
    use crate::simulate::rounded_poisson;
    use rand::SeedableRng;
    use rand_distr::Gamma;
    use rand_xoshiro::Xoshiro256PlusPlus;

    #[test]
    fn test_binned_outbreaks() {
        let dm = SimpleDisease {
            incubation_time: rounded_poisson(2.).unwrap(),
            reporting_time: rounded_poisson(2.).unwrap(),
            reproduction_number: Gamma::new(1.5, 0.75).unwrap(),
            infectiousness: vec![0.4, 0.275, 0.175, 0.1, 0.04, 0.01],
        };
        let mutation_rate = 1.1e-3 / 365. * 30000.;
        let sim_cfg = BinnedOutbreakConfig {
            size_bin_edges: vec![3, 7, 20, 150],
            size_counts: vec![4, 3, 2],
            latest_importation: 30,
            time_to_mrca: 30,
            time_to_background_mrca: 30,
            n_background: 5,
            bad_simulation_cap: 2000,
        };
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(89324 as u64);
        let genome = SimpleGenome::<64>::default();
        let outbreaks = binned_outbreaks(genome, &dm, mutation_rate, &sim_cfg, &mut rng).unwrap();

        assert_eq!(
            outbreaks.sources().len() - outbreaks.sources().iter().flatten().count(),
            4 + 3 + 2 + 5
        );
    }
}
