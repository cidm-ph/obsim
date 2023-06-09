use crate::genome::Genome;
use bitvec::prelude::*;
use rand::seq::index;
use rand::Rng;
use std::convert::TryInto;
use std::fmt;
use std::io;

type GenomeStorage = BitBox<usize, Lsb0>;

/// Representation of a genome.
///
/// Internally this is represented as an array of binary bits. Mutation is modelled by
/// selecting bits uniformly at random and flipping them.
///
/// This model enables efficient operations and compact storage.
#[derive(PartialEq, Eq, Clone)]
pub struct SimpleGenome<const BP: usize>(GenomeStorage);

impl<const BP: usize> Default for SimpleGenome<BP> {
    fn default() -> Self {
        SimpleGenome(bitbox![usize, Lsb0; 0; BP])
    }
}

impl<const BP: usize> Genome for SimpleGenome<BP> {
    /// Flip exactly `n_mutations` bits chosen at random.
    ///
    /// Panics if the requested number of mutations is greater than the fixed width of the
    /// representation.
    fn mutate<R: Rng>(&self, n_mutations: usize, mut rng: R) -> Self {
        assert!(
            n_mutations <= BP,
            "Requested number of mutations ({}) exceeds width of genome representation ({})",
            n_mutations,
            BP
        );
        let mut new_genome = self.0.clone();
        for pos in index::sample(&mut rng, BP, n_mutations) {
            let val = !new_genome.get(pos).unwrap();
            new_genome.set(pos, val);
        }
        Self(new_genome)
    }

    /// Counts the bitwise differences between the genome representations.
    fn snps(&self, other: &Self) -> u32 {
        (self.0.clone() ^ other.0.clone())
            .count_ones()
            .try_into()
            .unwrap()
    }

    /// Relabels 1 and 0 as A and C respectively.
    fn write_nucleotides<W: io::Write>(&self, mut writer: W) -> io::Result<()> {
        for base in self.0.as_bitslice() {
            write!(writer, "{}", if *base { 'A' } else { 'C' })?;
        }
        Ok(())
    }
}

impl<const BP: usize> fmt::Debug for SimpleGenome<BP> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SimpleGenome(")?;
        for base in self.0.as_bitslice() {
            write!(f, "{}", if *base { 'A' } else { 'C' })?;
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
        let genome1 = SimpleGenome::<64>::default();
        let genome2 = SimpleGenome::<64>::default();
        assert_eq!(genome1.snps(&genome2), 0);
    }

    #[test]
    fn test_mutation() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(89324 as u64);
        let genome = SimpleGenome::<64>::default();
        assert_ne!(genome, genome.mutate(4, &mut rng));
    }

    #[test]
    fn test_distance() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(89324 as u64);
        let genome = SimpleGenome::<64>::default();
        let child = genome.mutate(5, &mut rng);
        assert_eq!(genome.snps(&child), 5);
        assert_eq!(child.snps(&genome), 5);
    }
}
