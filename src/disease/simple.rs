use super::DiseaseModel;
use crate::case::CaseHistory;
use crate::Time;
use rand::Rng;
use rand_distr::Distribution;

/// Simple outbreak model.
///
/// `SimpleDisease` is a model off disease development. All cases have the same infectiousness
/// profile. There is a random incubation time when the case is not infectious, and then the onset
/// of symptoms and infectiousness occur simultaneously. Reporting occurs randomly after the onset
/// time.
///
/// Each case is assigned a random reproduction number that is used to scale the infectiousness
/// profile.
///
/// See [`rounded_poisson`](crate::simulate::rounded_poisson) and [`rand_distr`] for some useful distributions.
#[derive(Debug)]
pub struct SimpleDisease<DInc, DRep, DR> {
    /// A [`Distribution<Time>`](Distribution) of times between infection and infectiousness/symptom onset.
    pub incubation_time: DInc,

    /// A `Distribution<Time>` of times between symptom onset and case notification.
    pub reporting_time: DRep,

    /// A `Distribution<f64>` of individual reproduction numbers.
    pub reproduction_number: DR,

    /// Normally the infectiousness should curve should have an area of unity, i.e. the values
    /// in the vector should sum to unity. It is multiplied by the case's reproduction number
    /// to give the case's infectivity.
    ///
    /// The latent period at the start is added by the `incubation_time` and therefore the
    /// infectiousness should normally begin with a non-zero number.
    pub infectiousness: Vec<f64>,
}

impl<DInc, DRep, DR> DiseaseModel for SimpleDisease<DInc, DRep, DR>
where
    DInc: Distribution<Time>,
    DRep: Distribution<Time>,
    DR: Distribution<f64>,
{
    type State = ();

    fn generate_case<R: Rng>(&self, _state: &mut Self::State, mut rng: R) -> CaseHistory {
        let onset = self.incubation_time.sample(&mut rng);
        let reported = onset + self.reporting_time.sample(&mut rng);
        let r = self.reproduction_number.sample(&mut rng);
        let infect = std::iter::repeat(0.0)
            .take(onset as usize)
            .chain(self.infectiousness.iter().map(|x| x * r));

        CaseHistory {
            infectivity: infect.collect(),
            symptom_onset: Some(onset),
            reported: Some(reported),
        }
    }

    fn generate_singleton<R: Rng>(&self, mut rng: R) -> CaseHistory {
        let reported = self.incubation_time.sample(&mut rng);

        CaseHistory {
            infectivity: vec![],
            symptom_onset: Some(0),
            reported: Some(reported),
        }
    }
}
