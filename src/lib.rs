//! Simple simulation of a communicable disease outbreak.
//!
//! # Examples
//! The most basic simulation requires specifying distributions for all relevant quantities:
//!
//! ```
//! use obsim::simulate::{rounded_poisson, simulate_outbreak};
//! use obsim::simple::{SimpleDisease, SimpleGenome};
//! use rand::SeedableRng;
//! use rand_distr::Gamma;
//! use rand_xoshiro::Xoshiro256PlusPlus;
//!
//! let disease_model = SimpleDisease {
//!     // expected incubation and reporting times of 1 time step
//!     incubation_time: rounded_poisson(1.).unwrap(),
//!     reporting_time: rounded_poisson(1.).unwrap(),
//!     // each case is expected to infect about 1 other case
//!     reproduction_number: Gamma::new(1.5, 0.75).unwrap(),
//!     // cases have a constant infectiousness for 3 time steps after incubation
//!     infectiousness: vec![0.34, 0.33, 0.33],
//! };
//!
//! // expected mutations per time step
//! let mutation_rate = 2e-4 / 365.0 * 30000.0;
//!
//! // use a seedable RNG for reproducibility
//! let mut rng = Xoshiro256PlusPlus::seed_from_u64(893924 as u64);
//! // stop the simulation at the end of a time step where the case count exceeds 100
//! let max_cases = 100;
//!
//! let genome = SimpleGenome::<64>::default();
//! let ob = simulate_outbreak(genome, &disease_model, mutation_rate, max_cases, &mut rng);
//! assert!(ob.is_ok());
//!
//! let ob = ob.unwrap();
//! assert_eq!(ob.n_cases(), 5);
//! assert_eq!(ob.sources(), vec![None, Some(0), Some(1), Some(1), Some(1)]);
//! ```
//!
//! See the examples directory for more ways of configuring the simulations, e.g.
//! `cargo run --example combined`.

pub mod case;
mod disease;
mod genome;
pub mod simulate;

pub use disease::DiseaseModel;
pub use genome::Genome;
pub use simulate::outbreak::Outbreak;
pub use simulate::simulate_outbreak;

pub(crate) type Time = u32;
pub(crate) type Count = u32;

pub mod simple {
    //! Simple and efficient models that capture basic features.
    //!
    //! These models make a number of simplifying assumptions, but are a reasonable starting point
    //! and can be configured to match details of a pathogen of interest.

    pub use crate::disease::simple::SimpleDisease;
    pub use crate::genome::simple::SimpleGenome;
}

pub use disease::covid;
