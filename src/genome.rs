use rand::Rng;
use rand_distr::{Distribution, Poisson};
use std::io;

use super::Time;
pub mod simple;

/// Implemented by types that represent a genome sequence and a mutation model.
pub trait Genome: Clone {
    /// Induce the specified number of mutations.
    fn mutate<R: Rng>(&self, n_mutations: usize, rng: R) -> Self;

    /// Induce mutations according to an elapsed time.
    ///
    /// The number of mutations is chosen by a Poisson distribution.
    /// `generation_time` is the interval of time elapsed.
    /// `mutation_rate` is the expected number of mutations per unit time.
    fn mutate_time<R: Rng>(&self, generation_time: Time, mutation_rate: f64, mut rng: R) -> Self
    where
        Self: Sized,
    {
        let n_mutations = mutations_by_time(generation_time, mutation_rate, &mut rng);
        self.mutate(n_mutations, rng)
    }

    /// Measure SNP distance from another genome.
    fn snps(&self, other: &Self) -> u32;

    /// Represent the genome as a nucleotide string (suitable for FASTA).
    fn write_nucleotides<W: io::Write>(&self, writer: W) -> io::Result<()>;
}

fn mutations_by_time<R: Rng>(generation_time: Time, mutation_rate: f64, mut rng: R) -> usize {
    let lambda: f64 = f64::from(generation_time) * mutation_rate;
    if lambda <= 0.0 {
        0
    } else {
        Poisson::new(lambda).unwrap().sample(&mut rng) as usize
    }
}
