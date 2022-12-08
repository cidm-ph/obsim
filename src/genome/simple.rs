use crate::genome::Genome;
use bitvec::prelude::*;
use rand::seq::index;
use rand::Rng;
use std::convert::TryInto;
use std::fmt;
use std::io;

/// number of distinct sites on the genome that can be represented.
pub const GENOME_LENGTH: usize = 64;
pub type GenomeStorage = BitArr!(for GENOME_LENGTH, in u64);

/// Representation of a genome.
///
/// Internally this is represented as a vector of binary bits. Mutation is modelled by
/// selecting bits uniformly at random and flipping them.
///
/// This model enables efficient operations and compact storage.
#[derive(PartialEq, Eq, Clone, Copy, Default)]
pub struct SimpleGenome(GenomeStorage);

impl Genome for SimpleGenome {
    /// Flip exactly `n_mutations` bits chosen at random.
    ///
    /// Panics if the requested number of mutations is greater than the fixed width of the
    /// representation.
    fn mutate<R: Rng>(&self, n_mutations: usize, mut rng: R) -> Self {
        let mut new_genome = self.0;
        assert!(
            n_mutations <= GENOME_LENGTH,
            "Requested number of mutations ({}) exceeds width of genome representation ({})",
            n_mutations,
            GENOME_LENGTH
        );
        for pos in index::sample(&mut rng, GENOME_LENGTH, n_mutations) {
            let val = !new_genome.get(pos).unwrap();
            new_genome.set(pos, val);
        }
        Self(new_genome)
    }

    /// Counts the bitwise differences between the genome representations.
    fn snps(&self, other: &Self) -> u32 {
        (self.0 ^ other.0).count_ones().try_into().unwrap()
    }

    /// Relabels 1 and 0 as A and C respectively.
    fn write_nucleotides<W: io::Write>(&self, mut writer: W) -> io::Result<()> {
        for base in self.0 {
            write!(writer, "{}", if base { 'A' } else { 'C' })?;
        }
        Ok(())
    }
}

impl fmt::Debug for SimpleGenome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SimpleGenome(")?;
        for base in self.0 {
            write!(f, "{}", if base { 'A' } else { 'C' })?;
        }
        write!(f, ")")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    #[test]
    fn test_empty() {
        let genome1 = SimpleGenome::default();
        let genome2 = SimpleGenome::default();
        assert_eq!(genome1.snps(&genome2), 0);
    }

    #[test]
    fn test_mutation() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(89324 as u64);
        let genome = SimpleGenome::default();
        assert_ne!(genome, genome.mutate(4, &mut rng));
    }

    #[test]
    fn test_distance() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(89324 as u64);
        let genome = SimpleGenome::default();
        let child = genome.mutate(5, &mut rng);
        assert_eq!(genome.snps(&child), 5);
        assert_eq!(child.snps(&genome), 5);
    }
}
