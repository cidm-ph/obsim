//! Outbreak simulation.

use rand::distributions::{DistMap, WeightedIndex};
use rand::Rng;
use rand_distr::{Distribution, Poisson, PoissonError};
use thiserror::Error;

use crate::case::Case;
use crate::disease::DiseaseModel;
use crate::genome::Genome;
use crate::{Count, Time};

mod binned;
pub(super) mod outbreak;
pub use binned::{binned_outbreaks, BinError, BinnedOutbreakConfig};
use outbreak::Outbreak;

/// See [`rounded_poisson`].
pub type RoundedPoisson = DistMap<Poisson<f64>, fn(f64) -> Time, f64, Time>;

/// Generates a [`Poisson`] distribution with rate `lambda` that returns integer values, suitable
/// for use with disease model configuration.
#[inline]
pub fn rounded_poisson(lambda: f64) -> Result<RoundedPoisson, PoissonError> {
    Poisson::new(lambda).map(|x| x.map(round_time as fn(f64) -> Time))
}

fn round_time(t: f64) -> Time {
    t.round() as Time
}

#[derive(Error, Debug)]
#[error("outbreak exceeded {max_size} cases after time step")]
pub struct GrowthError<G> {
    pub outbreak: Outbreak<G>,
    pub max_size: Count,
}

/// Simulate an outbreak from one index genome.
///
/// The simulation stops when either the size exceeds `max_size` cases or
/// when there are no infectious cases left. These two cases are distinguished
/// by the result: `Ok()` indicates full recovery, and `Err()` indicates
/// termination due to `max_size`.
pub fn simulate_outbreak<D, G, R>(
    index_genome: G,
    disease_model: &D,
    mutation_rate: f64,
    max_size: Count,
    mut rng: R,
) -> Result<Outbreak<G>, GrowthError<G>>
where
    D: DiseaseModel,
    G: Genome,
    R: Rng,
{
    let mut dm_state = D::State::default();

    // start with the index case
    let (index, history) = disease_model
        .generate_case(&mut dm_state, &mut rng)
        .into_case_history();
    let mut t = 0;
    let mut cases = vec![index];
    let mut outbreak = Outbreak {
        source: vec![None],
        history: vec![history],
        genome: vec![index_genome],
    };

    loop {
        let case_infectivity: Vec<f64> = cases.iter_mut().map(Case::step).collect();
        let total_infectivity: f64 = case_infectivity.iter().sum();

        if total_infectivity > 0.0 {
            let case_dist = Poisson::new(total_infectivity).unwrap();
            let new_cases = case_dist.sample(&mut rng) as Count;

            if new_cases > 0 {
                let infector_dist = WeightedIndex::new(case_infectivity).unwrap();

                outbreak.source.reserve(new_cases as usize);
                outbreak.history.reserve(new_cases as usize);
                outbreak.genome.reserve(new_cases as usize);

                for _ in 0..new_cases {
                    let infector = infector_dist.sample(&mut rng) as Count;
                    outbreak.source.push(Some(infector));
                    let (case, mut history) = disease_model
                        .generate_case(&mut dm_state, &mut rng)
                        .into_case_history();
                    cases.push(case);
                    history.time_shift_forward(t);
                    outbreak.history.push(history);

                    let generation_time = t - outbreak.history[infector as usize].infected;
                    let new_genome = if generation_time < 1 {
                        outbreak.genome[infector as usize].clone()
                    } else {
                        outbreak.genome[infector as usize].mutate_time(
                            generation_time,
                            mutation_rate,
                            &mut rng,
                        )
                    };
                    outbreak.genome.push(new_genome);
                }
            }
        }

        if cases.iter().all(Case::is_recovered) {
            return Ok(outbreak);
        }

        if outbreak.n_cases() as Count > max_size {
            return Err(GrowthError { outbreak, max_size });
        }

        t += 1;
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
    fn test_simple_outbreak_simulation() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(893924 as u64);
        let dm = SimpleDisease {
            incubation_time: rounded_poisson(1.).unwrap(),
            reporting_time: rounded_poisson(1.).unwrap(),
            reproduction_number: Gamma::new(1.5, 0.75).unwrap(),
            infectiousness: vec![0.34, 0.33, 0.33],
        };
        let mutation_rate = 2e-4 / 365.;
        let genome = SimpleGenome::default();
        let outbreak = simulate_outbreak(genome, &dm, mutation_rate, 100, &mut rng);

        assert!(outbreak.is_ok());

        let outbreak = outbreak.unwrap();
        assert_eq!(outbreak.n_cases(), 5)
    }
}
