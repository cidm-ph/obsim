pub mod covid;
pub mod simple;

use crate::case::CaseHistory;
use rand::Rng;

/// Implemented by types that model the development of disease.
///
/// The model generates new cases on request, providing the minimal data necessary to reconstruct
/// the case history.
///
/// It can retain state between cases, and is rebuilt from scratch when simulating an independent
/// outbreak.
pub trait DiseaseModel {
    type State: Default;

    // pub trait DiseaseModel<const N: usize = 0> {
    // fn headings() -> [&'static str; N];
    // fn data(offset: Time) -> [Vec<f64>; N];
    fn generate_case<R: Rng>(&self, state: &mut Self::State, rng: R) -> CaseHistory;

    fn generate_singleton<R: Rng>(&self, rng: R) -> CaseHistory {
        let mut state = Self::State::default();
        self.generate_case(&mut state, rng)
    }
}
