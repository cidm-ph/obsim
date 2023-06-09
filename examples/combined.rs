use anyhow::Result;
use obsim::simple::{SimpleDisease, SimpleGenome};
use obsim::simulate::{rounded_poisson, simulate_outbreak};
use obsim::Genome;
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

    // expected mutations per time step
    let mutation_rate = 2e-4 / 365.0 * 30000.0;

    // halt the simulation if there are more than this number of cases
    let max_cases = 200;

    let mut rng = rand::thread_rng();
    // use this instead for reproducible simulation:
    // let mut rng = Xoshiro256PlusPlus::seed_from_u64(32189661_u64);

    // simulate an initial outbreak
    let index1 = SimpleGenome::<64>::default();
    let mut ob = simulate_outbreak(
        index1.clone(),
        &disease_model,
        mutation_rate,
        max_cases,
        &mut rng,
    )?;

    // create a second outbreak from a new introduction of the disease 30 time steps later with 10
    // SNPs separating it from the original index genome
    let index2 = index1.mutate(10, &mut rng);
    let mut ob2 = simulate_outbreak(index2, &disease_model, mutation_rate, max_cases, &mut rng)?;
    ob2.time_shift(30);

    // combine the outbreaks
    ob.extend_with(ob2);

    let stdout = std::io::stdout();
    ob.write_fasta(stdout.lock())?;

    Ok(())
}
