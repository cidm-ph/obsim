use anyhow::Result;
use obsim::simple::{SimpleDisease, SimpleGenome};
use obsim::simulate::{rounded_poisson, simulate_outbreak};
use rand_distr::Gamma;

// alternative RNG for reproducibility:
// use rand::SeedableRng;
// use rand_xoshiro::Xoshiro256PlusPlus;

fn main() -> Result<()> {
    let disease_model = SimpleDisease {
        incubation_time: rounded_poisson(2.)?,
        reporting_time: rounded_poisson(1.)?,
        reproduction_number: Gamma::new(2.8, 0.3)?,
        infectiousness: vec![0.34, 0.33, 0.33],
    };

    // expected mutations per time step
    let mutation_rate = 2e-4 / 365.0 * 30000.0;

    // halt the simulation if there are more than this number of cases
    let max_cases = 200;

    let mut rng = rand::thread_rng();
    // use this instead for reproducible simulation:
    // let mut rng = Xoshiro256PlusPlus::seed_from_u64(9948901_u64);

    let genome = SimpleGenome::<256>::default();
    let ob = simulate_outbreak(genome, &disease_model, mutation_rate, max_cases, &mut rng)?;

    let stdout = std::io::stdout();
    ob.write_fasta(stdout.lock())?;
    Ok(())
}
