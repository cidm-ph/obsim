use anyhow::Result;
use obsim::simple::{SimpleDisease, SimpleGenome};
use obsim::simulate::{binned_outbreaks, rounded_poisson, BinnedOutbreakConfig};
use rand_distr::Gamma;

// alternative RNG for reproducibility:
// use rand::SeedableRng;
// use rand_xoshiro::Xoshiro256PlusPlus;

fn main() -> Result<()> {
    let disease_model = SimpleDisease {
        incubation_time: rounded_poisson(2.)?,
        reporting_time: rounded_poisson(1.)?,
        reproduction_number: Gamma::new(2.5, 0.3)?,
        infectiousness: vec![0.34, 0.33, 0.33],
    };

    let binned_cfg = BinnedOutbreakConfig {
        size_bin_edges: vec![2, 10, 40, 180],
        size_counts: vec![4, 3, 2],
        latest_importation: 45,
        time_to_mrca: 7,
        time_to_background_mrca: 7,
        n_background: 20,
        bad_simulation_cap: 200,
    };

    // expected mutations per time step
    let mutation_rate = 2e-4 / 365.0 * 30000.0;

    let rng = rand::thread_rng();
    // use this instead for reproducible simulation:
    // let rng = Xoshiro256PlusPlus::seed_from_u64(32189661_u64);

    let genome = SimpleGenome::default();
    let ob = binned_outbreaks(genome, &disease_model, mutation_rate, &binned_cfg, rng)?;

    let stdout = std::io::stdout();
    ob.write_fasta(stdout.lock())?;

    Ok(())
}
